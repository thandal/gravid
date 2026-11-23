"""Forward validation of the exact adopted linear Gaussian problem.
Noise hyperparameters and Fourier basis fixed using full released model:
this is conditional validation, not independently trained noise selection.
"""
from pathlib import Path
import numpy as np,json,argparse
ap=argparse.ArgumentParser();ap.add_argument("--white-scale",type=float,default=1.);ap.add_argument("--prior-width",type=float,default=1.);ap.add_argument("--residuals",type=Path);args=ap.parse_args()
from scipy.linalg import lapack,solve_triangular,cholesky
from scipy.stats import chi2
root=Path.cwd();wd=root/'results/timing/j0437_predictive_audit';pre=root/'results/j0437_timing_precision_audit/native-export/solver'
n,p,nc=map(int,pre.with_suffix('.meta').read_text().splitlines()[0].split())
A=np.fromfile(str(pre)+'.white_design.bin').reshape(n,p);P=np.fromfile(str(pre)+'.constraints.bin').reshape(nc,p);y=np.fromfile(str(pre)+'.white_data.bin');yp=np.fromfile(str(pre)+'.constraint_values.bin')
if args.residuals:
 d=np.loadtxt(args.residuals);assert len(d)==n;y=d[:,2]/d[:,3]
A/=args.white_scale;y/=args.white_scale;P/=args.prior_width;yp/=args.prior_width
raw=[l.split() for l in (root/'results/timing/ppta_dr3_puzzle/csiro/dr3/J0437-4715.tim').read_text().splitlines() if len(l.split())>5 and l.split()[4]=='pks'];t=np.array([float(r[2]) for r in raw]);assert np.all(np.diff(t)>=0)
cut=t.max()-365.25;train=t<=cut;test=~train
B=np.array(np.vstack([A[train],P]),order='F');b=np.r_[y[train],yp];scale=np.linalg.norm(B,axis=0);assert np.all(scale>0),'Unidentifiable training parameter'
B/=scale
work,info=lapack.dgels_lwork(len(b),p,1);qr,fit,info=lapack.dgels(B,b,lwork=int(work),overwrite_a=True);assert info==0
R=np.triu(qr[:p]);coef=fit[:p]/scale;del qr,B
res=y[test]-A[test]@coef
K=solve_triangular(R.T,(A[test]/scale).T,lower=True)
V=np.eye(test.sum())+K.T@K
L=cholesky(V,lower=True);z=solve_triangular(L,res,lower=True)
q=float(z@z);nt=int(test.sum())
train_res=y[train]-A[train]@coef;prior_res=yp-P@coef
train_sse=float(train_res@train_res+prior_res@prior_res)
fullcoef=np.fromfile(str(pre)+'.parameters.bin');fullres=y-A@fullcoef;fullprior=yp-P@fullcoef
full_sse=float(fullres@fullres+fullprior@fullprior)
if args.white_scale==1 and args.prior_width==1 and args.residuals is None:assert abs(full_sse-train_sse-q)<1e-5
if args.white_scale==1 and args.prior_width==1 and args.residuals is None:print('Independent objective-difference check:',full_sse-train_sse,q,flush=True)
groups=np.array([r[r.index('-group')+1] for r in raw])[test]
group_results={}
for g in sorted(set(groups)):
 mask=groups==g;zg=solve_triangular(cholesky(V[np.ix_(mask,mask)],lower=True),res[mask],lower=True)
 group_results[g]=dict(n=int(sum(mask)),marginal_predictive_chi2=float(zg@zg),chi2_per_toa=float(zg@zg/sum(mask)))
orbital_rows=[l.split('\t') for l in (root/'scripts/timing/binary_model/covariance.tsv').read_text().splitlines() if l.startswith('parameter\t')]
orbital_ix=np.array([int(r[2]) for r in orbital_rows]);orbital_cov=np.array(json.loads((root/'gravid/resources/j0437_de440_binary_covariance.json').read_text())['covariance_tcb_internal'])
W=solve_triangular(L,A[test][:,orbital_ix]@np.linalg.cholesky(orbital_cov),lower=True)
u,sv,vh=np.linalg.svd(W,full_matrices=False);rank=int(sum(sv>sv[0]*1e-10));projection=float(np.sum((u[:,:rank].T@z)**2))
result=dict(white_scale=args.white_scale,prior_width=args.prior_width,binary_derivative_subspace=dict(rank=rank,squared_projection=projection,singular_values=sv.tolist(),scope="Diagnostic projection onto binary derivative shapes in held-out data; not a calibrated significance after predictive-noise failure."),training_n=int(train.sum()),heldout_n=nt,cut_mjd=float(cut),conditional_predictive_chi2=q,dof=nt,chi2_per_toa=q/nt,upper_tail_probability=float(chi2.sf(q,nt)),whitened_mean=float(z.mean()),whitened_std=float(z.std()),max_abs_whitened=float(max(abs(z))),scope='Forward prediction within exact linear timing/noise model; full-data noise hyperparameters and basis held fixed, not a fully independent predictive validation. All timing and noise coefficient uncertainty included. Cholesky residual order is chronological.',original_white_residual_rms=float(np.sqrt(np.mean(res**2))),source='native-export exact solver snapshot',objective_difference_check=full_sse-train_sse if args.white_scale==1 and args.prior_width==1 and args.residuals is None else None,training_augmented_sse=train_sse,full_augmented_sse=full_sse if args.white_scale==1 and args.prior_width==1 and args.residuals is None else None,largest_10_whitened_residuals_fraction=float(np.sort(z*z)[-10:].sum()/q),group_diagnostics=group_results)
suffix='' if args.white_scale==1 and args.prior_width==1 and args.residuals is None else f'_white{args.white_scale}_prior{args.prior_width}'
if args.residuals:suffix+='_nonlinear_residuals'
result['residuals_override']=str(args.residuals) if args.residuals else None
np.savez(wd/('heldout_linear_residuals'+suffix+'.npz'),mjd=t[test],residual=res,predictive_whitened=z,predictive_sigma=np.sqrt(V.diagonal()))
(wd/('heldout_linear'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
