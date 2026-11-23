"""Fixed-model O4ab prediction on portal TOAs; no TEMPO2 parameter fitting."""
from pathlib import Path
import os,json,subprocess,hashlib
from decimal import Decimal
import numpy as np
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'results/timing/portal_arrival_times';OUT=BASE/'validation';OUT.mkdir(exist_ok=True)
CLK=BASE/'clock';EXE=OUT/'evaluate'
subprocess.run(['g++','-O2','-I'+str(ROOT/'.conda/include'),str(ROOT/'scripts/timing/orbital_offset_fit/evaluate_timing.C'),'-L'+str(ROOT/'.conda/lib'),'-Wl,-rpath,'+str(ROOT/'.conda/lib'),'-ltempo2','-o',str(EXE)],check=True)
SWAP=set('RAJ DECJ PMRA PMDEC PX POSEPOCH F0 F1 F2 F3 PEPOCH BINARY PB PBDOT A1 XDOT A1DOT T0 OM ECC OMDOT M2 KIN KOM SINI TASC EPS1 EPS2 EPHEM UNITS PLANET_SHAPIRO'.split())
def lines(p):return [l for l in p.read_text().splitlines() if l.split() and not l.startswith('#')]
def pars(ls):return {l.split()[0]:l.split()[1] for l in ls if len(l.split())>=2}
results={}
for name in ['J0437-4715']:
 wd=OUT/name;wd.mkdir(exist_ok=True);rows=[]
 for p in sorted((BASE/'extracted'/name).rglob('*.tim')):
  for l in lines(p):
   a=l.split()
   if len(a)>4:
    try:t=float(a[2])
    except ValueError:continue
    if 60088.625<=t<=60703+17/24:rows.append(l)
 tim=wd/'o4ab.tim';tim.write_text('FORMAT 1\n'+'\n'.join(rows)+'\n');ids=np.array([l.split()[0] for l in rows])
 release=ROOT/'results/timing/j0437_external_validation/meerkat/partim'/f'{name}.par'
 adopted=ROOT/'gravid/resources/j0437_de440_search.par';published=ROOT/'data/J0437-4715_reardon2024_dr3.par'
 model={k:lines(p) for k,p in [('adopted',adopted),('published',published)]}
 pub=lines(published);pb=next(l.split() for l in pub if l.startswith('PB '));t0=next(l.split() for l in pub if l.startswith('T0 '))
 for key,parts in [('PB',pb),('T0',t0)]:
  model[key+'_plus_1sigma']=[f'{key} {Decimal(parts[1])+Decimal(parts[3])}' if l.split()[0]==key else l for l in pub]
 # Both comparisons use the fixed release DM/FD/instrument prescription; no retuning to portal TOAs.
 nuisance=[l for l in lines(release) if l.split()[0] not in SWAP and l.split()[0] not in ['START','FINISH','NITS','NTOA','CHI2R','TRES','CLK']]
 results[name]={}
 for label,source in model.items():
  par=wd/(label+'.par');par.write_text('\n'.join(nuisance+[l for l in source if l.split()[0] in SWAP and l.split()[0]!='PLANET_SHAPIRO']+['PLANET_SHAPIRO N','CLK TT(BIPM2024)'])+'\n')
  res=wd/(label+'.res')
  with (wd/(label+'.log')).open('w') as log:subprocess.run([str(EXE),str(par),str(tim),str(CLK),str(res)],env=dict(os.environ,TEMPO2=str(ROOT/'.conda/share/tempo2')),stdout=log,stderr=subprocess.STDOUT,check=True)
  d=np.loadtxt(res);assert len(d)==len(rows);p=pars(source);fgw=2*float(p['F0']);output={}
  for span,mask in [('O4a',d[:,0]<60325+16/24),('O4b',d[:,0]>=60410.625),('O4ab',(d[:,0]<60325+16/24)|(d[:,0]>=60410.625))]:
   dd=d[mask];ii=ids[mask];epochs=[]
   for ident in dict.fromkeys(ii):
    e=dd[ii==ident];w=e[:,3]**-2;z=np.average(np.exp(2j*np.pi*fgw*e[:,2]),weights=w)
    # Diagnostic per-epoch constant + nu^-2 fit: removes chromatic dispersion, not an orbital term.
    X=np.column_stack([np.ones(len(e)),(1000/e[:,4])**2]);A=X*np.sqrt(w[:,None]);b=e[:,2]*np.sqrt(w)
    beta=np.linalg.lstsq(A,b,rcond=None)[0];cov=np.linalg.inv(A.T@A);scatter=max(1.,np.sum((b-A@beta)**2)/max(1,len(e)-2))
    epochs.append([np.average(e[:,1],weights=w),z.real,z.imag,beta[0],np.sqrt(cov[0,0]*scatter)])
   ep=np.array(epochs);z=ep[:,1]+1j*ep[:,2];mean=z.mean();r=np.angle(z/mean)/(2*np.pi*fgw)
   row={'n_toas':len(dd),'n_epochs':len(ep),'raw_equal_epoch_power_loss':float(1-abs(mean)**2),'raw_epoch_rms_us':float(np.sqrt(np.mean(r*r))*1e6)}
   # Orbital harmonic test after projecting polynomial spin drift. All fitted terms disclosed.
   t=ep[:,0];tc=(t-t.mean())/365.25
   if 'TASC' in p:ell=2*np.pi*(t-float(p['TASC']))/float(p['PB'])
   else:
    u=(t-float(p['T0']))/float(p['PB']);ell=2*np.pi*(u-.5*float(p.get('PBDOT',0))*u*u)+np.deg2rad(float(p['OM'])+float(p.get('OMDOT',0))*(t-float(p['T0']))/365.25)
   X=np.column_stack([np.ones(len(t)),tc,tc**2,np.sin(ell),np.cos(ell)])
   if len(t)>5:
    for weighting in ['equal','formal']:
     w=np.ones(len(t)) if weighting=='equal' else ep[:,4]**-2;A=X*np.sqrt(w[:,None]);b=ep[:,3]*np.sqrt(w);beta=np.linalg.lstsq(A,b,rcond=None)[0];inv=np.linalg.inv(A.T@A);residual=b-A@beta;h=np.sum((A@inv)*A,axis=1);hc3=inv@(A.T@np.diag((residual/(1-h))**2)@A)@inv
     row[weighting+'_chromatic_projected_orbit']={'sin_delay_us':float(beta[-2]*1e6),'cos_delay_us':float(beta[-1]*1e6),'HC3_sin_sigma_us':float(np.sqrt(hc3[-2,-2])*1e6),'HC3_cos_sigma_us':float(np.sqrt(hc3[-1,-1])*1e6),'post_regression_rms_us':float(np.sqrt(np.mean((ep[:,3]-X@beta)**2))*1e6)}
   np.savetxt(wd/(label+'_'+span+'_epochs.txt'),ep,header='bat_mjd complex_phase_real complex_phase_imag achromatic_intercept_s intercept_error_s');output[span]=row
  results[name][label]=output;print(name,label,json.dumps(output),flush=True)
(OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
