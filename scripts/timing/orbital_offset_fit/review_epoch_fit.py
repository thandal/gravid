import numpy as np,json,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path('results/timing/j0437_external_validation'));p.add_argument('--output',type=Path,default=Path('scripts/timing/orbital_offset_fit/orbital_fit_review.json'));args=p.parse_args()
b=args.inputs;out={}
for site in ['meerkat','inpta']:
 for sub in ['post','o4a','all']:
  a=np.load(b/f'epoch_{site}_{sub}/epochs.npz');t,r,s,ell=[a[k] for k in ['t','r','s','ell']];tc=(t-np.average(t,weights=s**-2))/365.25
  X=np.column_stack([np.ones(len(t)),tc]+([] if sub=='o4a' else [tc**2])+[np.sin(ell),np.cos(ell)])
  rows={}
  for mode in ['published_weights','equal_epoch_weights']:
   w=s**-2 if mode=='published_weights' else np.ones(len(t));w/=w.mean();A=X*np.sqrt(w[:,None]);y=r*np.sqrt(w)
   beta=np.linalg.lstsq(A,y,rcond=None)[0];inv=np.linalg.inv(A.T@A);e=y-A@beta;h=np.sum((A@inv)*A,axis=1)
   cov=inv*(e@e)/(len(t)-X.shape[1]);hc=inv@(A.T@((e/(1-h))[:,None]**2*A))@inv
   loo=[]
   for i in range(len(t)):
    m=np.arange(len(t))!=i;loo.append(np.linalg.lstsq(A[m],y[m],rcond=None)[0][-1]*1e9)
   # Bartlett HAC estimates are sensitivity checks, not proof of stationarity.
   scores=(A*e[:,None])[np.argsort(t)];hac={}
   for lag in [1,2,4]:
    meat=scores.T@scores
    for k in range(1,min(lag+1,len(t))):
     cross=scores[k:].T@scores[:-k];meat+=(1-k/(lag+1))*(cross+cross.T)
    c=inv@meat@inv*len(t)/(len(t)-X.shape[1]);hac[str(lag)]=float(np.sqrt(max(0,c[-1,-1]))*1e9)
   rows[mode]={'longitude_delay_ns':beta[-1]*1e9,'scaled_sigma_ns':np.sqrt(cov[-1,-1])*1e9,'HC3_sigma_ns':np.sqrt(hc[-1,-1])*1e9,'HAC_sigma_ns_by_epoch_lag':hac,'leave_one_out_range_ns':[min(loo),max(loo)],'dof':len(t)-X.shape[1]}
  out[site+'_'+sub]=rows
print(json.dumps(out,indent=2));args.output.write_text(json.dumps(out,indent=2)+'\n')
