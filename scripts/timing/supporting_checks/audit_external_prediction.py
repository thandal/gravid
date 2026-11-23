"""Audit radio-sampled coherence and project out per-day dispersion.

No spin/binary parameter is fitted. The achromatic intercept diagnostic DOES
fit per-day dispersion; it tests sensitivity to the released DM corrections.
It is not the detector-weighted CW overlap or a confidence bound.
"""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];EXT=ROOT/'results/timing/j0437_external_validation';SRC=ROOT/'results/timing/first_principles';OUT=ROOT/'results/timing/supporting_checks'
p={l.split()[0]:l.split()[1] for l in (ROOT/'data/J0437-4715_reardon2024_dr3.par').read_text().splitlines() if len(l.split())>1}
f=2*float(p['F0']);pb=float(p['PB']);t0=float(p['T0']);out={}
for site in ['meerkat','inpta']:
 names=np.array([l.split()[0] for l in (EXT/site/'flattened.tim').read_text().splitlines() if len(l.split())>4 and l.split()[0] not in ('FORMAT','MODE')]);res={}
 for label in ['published_DR3','control_PB_p1_sigma','control_PB_p1s_OMDOT_partner','control_T0_p1_sigma','control_T0_p1s_OM_partner','old_search_ELL1_(no_OMDOT)']:
  path=SRC/f'{site}_{label}.res';d=np.loadtxt(path);sat,bat,r,err,freq=d.T
  mask=(sat>=60079.59331597222)&(sat<=60325.65554976852);sat,bat,r,err,freq,ids=[a[mask] for a in [sat,bat,r,err,freq,names]]
  epochs=list(dict.fromkeys(ids));z=np.exp(2j*np.pi*f*r)
  ez=np.array([np.average(z[ids==e],weights=err[ids==e]**-2) for e in epochs]);dr=np.angle(ez/ez.mean())/(2*np.pi*f)
  # Independent broadband diagnostic: discard the release's DM correction
  # by marginalizing an arbitrary nu^-2 coefficient every observing day.
  # Residuals of unperturbed models are far smaller than pulse-wrap spacing.
  day=np.floor(sat).astype(int);a=[];aerr=[];times=[];freqranges=[]
  for ep in np.unique(day):
   ix=day==ep;x=(1000/freq[ix])**2;X=np.column_stack([np.ones(sum(ix)),x]);w=1/err[ix]**2
   normal=(X.T*w)@X;cov=np.linalg.inv(normal);coef=cov@(X.T@(w*r[ix]));a.append(coef[0]);aerr.append(np.sqrt(cov[0,0]));times.append(float(np.mean(bat[ix])));freqranges.append([float(min(freq[ix])),float(max(freq[ix]))])
  a=np.array(a);times=np.array(times);phase=np.remainder((times-t0)/pb,1);angle=2*np.pi*phase;X=np.column_stack([np.ones(len(phase)),np.sin(angle),np.cos(angle)]);gram=X.T@X/len(X)
  row=dict(toas=len(r),file_groups=len(epochs),observing_days=len(a),sat_range=[float(min(sat)),float(max(sat))],file_group_coherence_loss=float(1-abs(ez.mean())**2),file_group_rms_ns=float(np.std(dr)*1e9),day_dm_projected_intercept_rms_ns=float(np.std(a)*1e9),day_dm_projected_median_formal_error_ns=float(np.median(aerr)*1e9),day_dm_projected_coherence_loss=float(1-abs(np.exp(2j*np.pi*f*a).mean())**2),orbital_sampling_gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),orbital_cycles=phase.tolist(),frequency_ranges_MHz=freqranges,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
  if 'control' in label or 'old_search' in label:row['dm_projection_warning']='Unwrapped linear intercepts can be ambiguous for large-control residuals; use the direct complex coherence for these controls.'
  official=(d[:,0]>=60088.625)&(d[:,0]<=60325.65554976852)
  oi=names[official];oz=np.exp(2j*np.pi*f*d[official,2]);ow=d[official,3]**-2
  oe=np.array([np.average(oz[oi==e],weights=ow[oi==e]) for e in dict.fromkeys(oi)])
  row['official_O4a']={'file_groups':len(oe),'observing_days':len(np.unique(np.floor(d[official,0]))),'rms_ns':float(np.std(np.angle(oe/oe.mean())/(2*np.pi*f))*1e9),'coherence_loss':float(1-abs(oe.mean())**2)}
  res[label]=row
 out[site]=res
(OUT/'external_prediction_audit.json').write_text(json.dumps(out,indent=2)+'\n')
for site,res in out.items():
 print(site,json.dumps({k:{kk:v for kk,v in r.items() if kk not in ['orbital_cycles','frequency_ranges_MHz','sha256']} for k,r in res.items()},indent=2))
