"""Epoch-averaged external test of the PPTA orbital prediction.

1. TEMPO2 fits only spin/dispersion nuisance terms (F0, F1, DM...) on the
   external TOA subset with the PPTA DE440 solution otherwise fixed.
2. Post-fit residuals are evaluated with the no-fit evaluator and averaged
   per observation (filename = epoch and band).
3. A weighted linear fit of the epoch means to
      r = c0 + c1 (t - tbar) + c2 (t-tbar)^2 + dA1 sin(l) + A1 dl cos(l)
   with errors scaled by the epoch-level reduced chi-square gives the
   orbital-longitude offset dl and A1 offset with conditional scaled errors.
"""
import json, os, re, subprocess
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
BASE = Path(os.environ.get('EXTERNAL_AUDIT_DIR', ROOT/'results/timing/j0437_external_validation'))
CLK = BASE/'clock'; T2 = ROOT/'.conda/bin/tempo2'; EVAL = BASE/'evaluate_timing'
PPTA_END = 59645.40317295509; O4A = (60079.59331597222, 60325.65554976852)
degyr = np.pi/180/(365.25*86400)
PB_d = 5.7410463602832067856; T0_d = 54530.172210721427433; A1 = 3.3667146504782899253
OM_rad = np.deg2rad(1.3568640571311364396); OMDOT_rpd = 0.015561777438539341396*degyr*86400
PBDOT = 3.7326281719829273242e-12
C = np.array(json.load(open(ROOT/'gravid/resources/j0437_de440_binary_covariance.json'))['covariance_tcb_internal'])
EQ_KOM = -5.7965e-05; EQ_KIN = 2.5096e-06
n_rpd = 2*np.pi/PB_d + OMDOT_rpd
def ppta_sigma_longitude(t_mjd):
    dt = t_mjd - T0_d; u = dt/PB_d; J = np.zeros(10)
    J[0] = -2*np.pi*dt/PB_d**2 + 2*np.pi*PBDOT*u**2/PB_d; J[3] = 1.0; J[4] = -n_rpd; J[5] = -np.pi*u**2
    J[6] = degyr*86400*dt; J[7] = EQ_KIN*degyr*86400*dt; J[8] = EQ_KOM*degyr*86400*dt
    return float(np.sqrt(J@C@J))

def read_toas(tim):
    hdr, rows = [], []
    for line in tim.read_text().splitlines():
        p = line.split()
        if len(p) > 4 and p[0] not in ('FORMAT', 'MODE', 'C', '#'): rows.append((float(p[2]), float(p[3]), p[0], line))
        else: hdr.append(line)
    return hdr, rows
def make_par(src, fit, dest):
    out = []
    for line in src.read_text().splitlines():
        f = line.split()
        if not f or f[0] in ('START','FINISH','TRES','CHI2R','NTOA','NITS'): continue
        if f[0] in ('F0','F1','F2','RAJ','DECJ','PMRA','PMDEC','PX','PB','T0','A1','OM','ECC','PBDOT','OMDOT','M2','KOM','KIN','DM','DM1','DM2','FD1','FD2') and len(f) >= 3:
            flag = '1' if f[0] in fit else '0'
            if f[2] in ('0','1'): f[2] = flag
            else: f.insert(2, flag)
            line = ' '.join(f)
        out.append(line)
    out.append('NITS 3'); dest.write_text('\n'.join(out)+'\n')

def epoch_fit(name, tag, subset, fit, quad=True):
    wd = BASE/f'epoch_{name}_{tag}'; wd.mkdir(exist_ok=True)
    hdr, rows = read_toas(BASE/name/'flattened.tim'); rows = [r for r in rows if subset(r[0])]
    (wd/'input.tim').write_text('\n'.join(hdr+[r[3] for r in rows])+'\n')
    make_par(BASE/name/'replace_all.par', fit, wd/'input.par')
    env = dict(os.environ, TEMPO2=str(ROOT/'.conda/share/tempo2'), OPENBLAS_NUM_THREADS='2')
    log = subprocess.run([str(T2),'-nobs','20000','-clkdir',str(CLK),'-f','input.par','input.tim','-qrfit','-newpar'], cwd=wd, env=env, capture_output=True, text=True, check=True)
    (wd/'fit.log').write_text(log.stdout+log.stderr)
    chi = re.findall(r'Chisqr/nfree\s*=\s*([\d.eE+-]+)/(\d+)\s*=\s*([\d.eE+-]+)', log.stdout)
    evaluated = subprocess.run([str(EVAL), str(wd/'new.par'), str(wd/'input.tim'), str(CLK), str(wd/'resid.txt')], env=env, check=True, capture_output=True, text=True)
    (wd/'evaluate.log').write_text(evaluated.stdout+evaluated.stderr)
    d = np.loadtxt(wd/'resid.txt'); assert len(d) == len(rows) and np.max(abs(d[:,0]-np.array([r[0] for r in rows]))) < 1e-6
    sat, bat, res, err = d[:,0], d[:,1], d[:,2], d[:,3]
    names = np.array([r[2] for r in rows]); w = err**-2
    groups = list(dict.fromkeys(names))
    t_e, r_e, s_e, n_e = [], [], [], []
    for g in groups:
        m = names == g
        t_e.append(np.average(bat[m], weights=w[m])); r_e.append(np.average(res[m], weights=w[m]))
        s_e.append(np.sqrt(1/w[m].sum())); n_e.append(int(m.sum()))
    t_e, r_e, s_e = map(np.array, (t_e, r_e, s_e))
    # orbital longitude from PPTA model at each epoch (TCB days)
    dt = t_e - T0_d; u = dt/PB_d
    ell = 2*np.pi*(u - 0.5*PBDOT*u*u) + OM_rad + OMDOT_rpd*dt
    tc = (t_e - np.average(t_e, weights=s_e**-2))/365.25
    cols = [np.ones_like(tc), tc] + ([tc**2] if quad else []) + [np.sin(ell), A1*np.cos(ell)]
    X = np.column_stack(cols); W = s_e**-2
    def solve(W):
        matrix=X*np.sqrt(W[:,None]); target=r_e*np.sqrt(W)
        beta, _, rank, _ = np.linalg.lstsq(matrix,target,rcond=None)
        if rank != X.shape[1]: raise ValueError('Rank-deficient orbital regression')
        _, triangular = np.linalg.qr(matrix,mode='reduced')
        inv=np.linalg.solve(triangular,np.eye(len(beta)))
        return beta, inv@inv.T
    beta, cov = solve(W)
    chi2 = float(np.sum(W*(r_e-X@beta)**2)); dof = len(r_e)-X.shape[1]; scale = max(1.0, chi2/dof)
    sig = np.sqrt(np.diag(cov)*scale)
    i_sin, i_cos = X.shape[1]-2, X.shape[1]-1
    tmid = float(np.average(t_e, weights=W)); sp = ppta_sigma_longitude(tmid)
    out = dict(n_toas=len(rows), n_epochs=len(groups), weighted_mid_mjd=tmid, tempo2_chi2r=float(chi[-1][2]) if chi else None,
               fit_nuisance=sorted(fit), epoch_chi2r=chi2/dof, epoch_rms_us=float(np.sqrt(np.average((r_e-X@beta)**2, weights=W))*1e6),
               median_epoch_nominal_err_ns=float(np.median(s_e)*1e9),
               dA1_ns=beta[i_sin]*1e9, dA1_sigma_ns=sig[i_sin]*1e9, dA1_sigma_ppta_ns=float(np.sqrt(C[1,1])*1e9),
               dl_rad=beta[i_cos], dl_sigma_ext_rad=sig[i_cos], dl_sigma_ppta_rad=sp, dl_z=beta[i_cos]/np.hypot(sig[i_cos], sp),
               dl_roemer_ns=beta[i_cos]*A1*1e9, dl_sigma_roemer_ns=sig[i_cos]*A1*1e9, gw_phase_amp_rad=2*np.pi*347.3759*A1*beta[i_cos])
    np.savez(wd/'epochs.npz', t=t_e, r=r_e, s=s_e, n=np.array(n_e), ell=ell)
    print(name, tag, json.dumps(out, default=float), flush=True); return out

def main():
    R = {}
    for name, nuis in (('meerkat', {'DM','DM1'}), ('inpta', set())):
        R[name] = {}
        R[name]['post'] = epoch_fit(name, 'post', lambda t: t > PPTA_END, {'F0','F1'}|nuis)
        R[name]['o4a'] = epoch_fit(name, 'o4a', lambda t: O4A[0] <= t <= O4A[1], {'F0'}|({'DM'} if nuis else set()), quad=False)
        R[name]['all'] = epoch_fit(name, 'all', lambda t: True, {'F0','F1'}|nuis|({'DM2'} if nuis else set()))
    (BASE/'epoch_orbit_fit.json').write_text(json.dumps(R, indent=1, default=float)+'\n')

if __name__ == '__main__': main()
