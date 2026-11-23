from pathlib import Path
import sys,json,hashlib
import numpy as np
root=Path.cwd();sys.path.insert(0,str(root))
from scripts.timing import build_j0437_search_model as builder
builder.KEYS.update(["F2","A1DOT","XDOT"])
build=builder.build
from gravfold.j0437_model import SEARCH_PAR,COVARIANCE,KEYS,binary_mean,apply_binary,waveform_metric
from gravfold.j0437_campaign import load_chunk
from gravfold.cwinpy_adapter import _ensure_template_amplitude
from cwinpy.parfile import PulsarParameters
from cwinpy.signal import HeterodynedCWSimulator
import pint.logging
pint.logging.setup(level='ERROR')
wd=root/'results/timing/supporting_checks/refits/J0437-4715';cache=root/'results/timing/j0437_predictive_audit/samples.npz'
if not cache.exists():
 blocks={}
 for path in sorted((root/'results/campaign_het_v2').glob('het_*_g*.hdf5')):
  det,t,d,_=load_chunk(path);ix=np.linspace(0,len(t)-1,32,dtype=int)
  blocks.setdefault(det,[]).append((t[ix],np.full(len(ix),len(t)/len(ix)/np.var(d))))
 np.savez(cache,**{det+'_'+key:np.concatenate([b[i] for b in v]) for det,v in blocks.items() for i,key in enumerate(['times','weights'])})
a=np.load(cache);samples=[(det,a[det+'_times'],a[det+'_weights']) for det in ['H1','L1']]
base=_ensure_template_amplitude(PulsarParameters(str(SEARCH_PAR)));sims=[HeterodynedCWSimulator(base,det,t) for det,t,w in samples]
w=np.concatenate([w*abs(sim.model(newpar=base))**2 for sim,(_,t,w) in zip(sims,samples)]);w/=sum(w)
reference_cov=json.loads(COVARIANCE.read_text());oldraw=np.array(reference_cov['covariance_tcb_internal']);factors=np.sqrt(np.diag(reference_cov['covariance'])/oldraw.diagonal())
summary={}
for case in sorted(wd.glob('*/DONE.json')):
 if json.loads(case.read_text())['returncode']!=0:continue
 d=case.parent;source=d/'precise.par'
 if not source.exists():source=d/'native.par'
 if not source.exists():source=d/'new.par'
 if not source.exists():continue
 artifact=dict(reference_cov);artifact['provenance']=dict(artifact['provenance'],fit_par_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
 covpath=d/'covariance_for_conversion.json';covpath.write_text(json.dumps(artifact))
 text,meta=build(source,covpath);search=d/'search.par';search.write_text(text)
 q=PulsarParameters(str(search));binary=apply_binary(PulsarParameters(str(SEARCH_PAR)),binary_mean(q))
 if q['XDOT'] is not None:binary['XDOT']=q['XDOT']
 row={}
 for label,par,ssb in [('total',q,True),('binary_only',binary,False)]:
  phase=np.concatenate([2*np.pi*sim.model(newpar=par,updateSSB=ssb,updateBSB=True,freqfactor=2.,phase_only=True) for sim in sims]);phase=np.angle(np.exp(1j*phase));phase-=phase@w
  row[label]=dict(phase_span_rad=float(np.ptp(phase)),phase_rms_rad=float(np.sqrt(w@phase**2)),power_loss=float(max(0,1-abs(w@np.exp(1j*phase))**2)))
  np.save(d/(label+'_phase.npy'),phase)
 if (d/'covariance.tsv').exists():
  lines=(d/'covariance.tsv').read_text().splitlines();C=np.loadtxt(lines[lines.index('covariance_internal')+1:]);canon=C*np.outer(factors,factors)
  g=waveform_metric(search,canon,samples)
  if (d/'extended_covariance.tsv').exists():
   linesextra=(d/'extended_covariance.tsv').read_text().splitlines();key,value,err,rawerr=linesextra[0].split();key='XDOT' if key=='A1DOT' else key;extended=np.loadtxt(linesextra[1:])
   scale_tdb=factors[0]/86400
   extra_factor=float(err)/float(rawerr)*(-1/scale_tdb**3 if key=='F2' else 1.)
   extended*=np.outer(np.r_[factors,extra_factor],np.r_[factors,extra_factor])
   L=np.linalg.cholesky(extended);all_derivs=[];all_weights=[]
   qb=_ensure_template_amplitude(PulsarParameters(str(search)));mean=binary_mean(qb)
   for det,t,weight in samples:
    sim=HeterodynedCWSimulator(qb,det,t);all_weights.append(weight*abs(sim.model(newpar=qb))**2);derivs=[]
    for j in range(11):
     phase=[]
     for sign in [-1,1]:
      delta=sign*.1*L[:,j];pp=apply_binary(PulsarParameters(str(search)),mean+delta[:10]);pp[key]=float(qb[key])+delta[-1]
      phase.append(2*np.pi*sim.model(newpar=pp,updateSSB=False,updateBSB=True,phase_only=True))
     derivs.append(np.angle(np.exp(1j*(phase[1]-phase[0])))/.2)
    all_derivs.append(derivs)
   weights=np.concatenate(all_weights);weights/=weights.sum();D=np.concatenate(all_derivs,axis=1);D-=(D@weights)[:,None];g=(D*weights)@D.T
   row['metric_parameter_count']=11;row['extra_parameter']=key;row['extra_parameter_sigma_physical']=float(err);row['extra_internal_to_physical_factor']=extra_factor
  elif any(x in d.name for x in ['fit_F2','fit_A1DOT']):
   row['metric_scope_warning']='Extra fitted parameter excluded from this ten-parameter diagnostic; do not interpret as full phase uncertainty.'
  row['new_metric_max_eigenvalue']=float(np.linalg.eigvalsh(g)[-1]);row['new_metric_trace']=float(np.trace(g));row['pb_omdot_correlation']=float(C[0,6]/np.sqrt(C[0,0]*C[6,6]));row['binary_sigma_over_adopted']=dict(zip(KEYS,np.sqrt(C.diagonal()/oldraw.diagonal()).tolist()))
 if 'fit_F2' in d.name:
  for k in ['new_metric_max_eigenvalue','new_metric_trace']:row.pop(k,None)
  row['metric_scope_warning']='F2 is correlated with F0 and F1, which were not included in this covariance export. No full phase-uncertainty result is claimed for this variant; only its refitted central phase is compared.'
 row['fit_source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();summary[d.name]=row
 print(d.name,json.dumps(row),flush=True)
(wd/'phase_stability.json').write_text(json.dumps(summary,indent=2)+'\n')
