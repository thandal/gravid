#!/usr/bin/env python3
"""Compare uncertainty propagation in independent PINT, TEMPO2 and LAL binaries.

Each implementation applies the canonical perturbations with its own unit and
angle conversion. Shared barycentric epochs/Earth positions isolate the binary
model. No refit is allowed. Deliberately broken covariance directions are
positive controls; agreement on the central waveform alone is insufficient.
"""
import argparse
import hashlib
import itertools
import json
import os
from importlib import metadata
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import astropy.units as u
    from astropy.constants import c
    import lal,lalpulsar,pint.logging
    from pint.models import get_model
    from cwinpy.parfile import PulsarParameters
    from cwinpy.signal import HeterodynedCWSimulator
    from gravid.j0437_model import SEARCH_PAR,COVARIANCE,KEYS,load_covariance,binary_mean,apply_binary
    from gravid.cwinpy_adapter import _ensure_template_amplitude
    from gravid.j0437_campaign import load_chunk
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick',action='store_true')
    ap.add_argument('--work-dir',type=Path,default=ROOT/'results/timing/article1_triple_check')
    ap.add_argument('--out',type=Path,default=ROOT/'scripts/timing/binary_model/uncertainty_crosscheck.json')
    args=ap.parse_args();pint.logging.setup(level='ERROR')
    wd=args.work_dir.resolve();wd.mkdir(exist_ok=True,parents=True)
    par=_ensure_template_amplitude(PulsarParameters(str(SEARCH_PAR)))
    cov,artifact=load_covariance();factor=np.linalg.cholesky(cov);mean=binary_mean(par)
    pintmodel=get_model(str(SEARCH_PAR),allow_T2=True)
    component=pintmodel.components['BinaryDDK'];component.update_binary_object(None)
    binary=component.binary_instance
    # Independent increments in PINT units and its DT96 angular convention.
    units=[u.s,u.lightsecond,u.dimensionless_unscaled,u.rad,u.day,u.dimensionless_unscaled,u.deg/u.yr,u.rad,u.rad,u.Msun]
    scales=np.array([1,1,1,1,1/86400,1e-12,1,-1,-1,1],dtype=np.longdouble)
    pintbase={key:getattr(binary,key).copy() for key in KEYS}
    sampled={}
    for path in sorted((ROOT/'results/campaign_het_v2').glob('het_*_g*.hdf5')):
        det,t,d,_=load_chunk(path)
        ix=np.linspace(0,len(t)-1,min(128,len(t)),dtype=int)[::8]
        sampled.setdefault(det,[]).append((t[ix],np.full(len(ix),len(t)/len(ix)/np.var(d))))
    sims=[];weights=[];epochs=[];earths=[];psrvecs=[];barygps=[];base_lal=[]
    for det,blocks in sorted(sampled.items()):
        times=np.concatenate([b[0] for b in blocks]);w=np.concatenate([b[1] for b in blocks])
        if args.quick:times=times[::8];w=w[::8]
        sim=HeterodynedCWSimulator(par,det,times);sims.append(sim)
        weights.append(w*abs(sim.model(newpar=par))**2)
        gps=times.astype(np.longdouble)+sim.ssbdelay.data.astype(np.longdouble)
        tb=np.longdouble(pintmodel.T0.value)+(gps-np.longdouble(par['T0']))/86400
        pos=[]
        for t in times:
            state=lalpulsar.EarthState();lalpulsar.BarycenterEarth(state,lal.LIGOTimeGPS(float(t)),sim._HeterodynedCWSimulator__edat)
            pos.append(state.posNow.copy())
        epochs.append(tb);earths.append(np.array(pos));psrvecs.append(pintmodel.ssb_to_psb_xyz_ICRS(epoch=tb.astype(float)))
        barygps.append(gps);base_lal.append(sim.bsbdelay.data.copy())
    weights=np.concatenate(weights);weights/=weights.sum();tb=np.concatenate(epochs);earth=np.concatenate(earths)
    psrvec=np.asarray(np.concatenate(psrvecs));bgps=np.concatenate(barygps);base_lal=np.concatenate(base_lal)
    binary.update_input(barycentric_toa=tb*u.day,obs_pos=earth*u.s*c,psr_pos=psrvec)
    base_pint=binary.binary_delay().to_value(u.s).copy()
    # Same full-covariance corners as the prior validation, plus independent draws.
    prior=json.loads((ROOT/'scripts/timing/binary_model/model_validation.json').read_text())
    _,rotation=np.linalg.eigh(prior['metric'])
    rng=np.random.default_rng(20261016)
    cases=[('central',np.zeros(10))]
    for j,key in enumerate(KEYS):
        d=np.zeros(10);d[j]=np.sqrt(cov[j,j]);cases.append(('single_'+key,d))
    for j in range(10):
        for sign in [-1,1]:cases.append((f'correlated_axis_{j}_{sign}',sign*factor[:,j]))
    for i,z in enumerate(rng.normal(size=(8 if args.quick else 256,10))):cases.append((f'gaussian_{i}',factor@z))
    if not args.quick:
        for i,y in enumerate(itertools.product((-3.,3.),repeat=10)):cases.append((f'corner_{i}',factor@rotation@np.array(y)))
    # Remove OMDOT from a correlated draw that shifts PB by exactly one sigma.
    d=factor[:,0].copy();d[6]=0;cases.append(('broken_PB_OMDOT',d))
    d=np.zeros(10);d[[0,6]]=np.sqrt(np.diag(cov))[[0,6]];cases.append(('PB_OMDOT_pair_only',d))
    for eps in [.1,1.]:
        for j in range(10):
            for sign in [-1,1]:cases.append((f'metric_{eps}_{j}_{sign}',sign*eps*factor[:,j]))
    # TEMPO2 inputs: high precision TDB MJDs, Earth light-seconds, unit sky vectors;
    # increments use the same canonical units but C++ performs its own conversion.
    with (wd/'samples.txt').open('w') as f:
        f.write(str(len(tb))+'\n')
        np.savetxt(f,np.column_stack([tb,earth,psrvec]),fmt='%.21g')
    with (wd/'draws.txt').open('w') as f:
        f.write(str(len(cases))+'\n');np.savetxt(f,np.array([d for _,d in cases]),fmt='%.21g')
    source=ROOT/'scripts/timing/check_j0437_tempo2_delays.C';exe=wd/'tempo2_delays'
    subprocess.run(['g++','-O2','-I'+str(ROOT/'.conda/include'),str(source),'-L'+str(ROOT/'.conda/lib'),'-Wl,-rpath,'+str(ROOT/'.conda/lib'),'-ltempo2','-o',str(exe)],check=True)
    with (wd/'tempo2.log').open('w') as log:
        subprocess.run([str(exe),str(SEARCH_PAR),str(wd/'samples.txt'),str(wd/'draws.txt'),str(wd/'delays.bin'),str(wd/'tempo2_parameters.txt')],check=True,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,TEMPO2=str(ROOT/'.conda/share/tempo2')))
    tempo=np.fromfile(wd/'delays.bin').reshape(len(cases),len(tb));base_tempo=tempo[0]
    frequency=2*(np.longdouble(par['F0'])+np.longdouble(par['F1'])*(bgps+base_lal-np.longdouble(par['PEPOCH'])))
    rows=[];metrics={};central_error=float(max(abs(base_pint+base_lal)))
    for i,(name,delta) in enumerate(cases):
        updates={key:pintbase[key]+np.longdouble(delta[j])*scales[j]*units[j] for j,key in enumerate(KEYS)}
        binary.update_input(**updates)
        pd=np.asarray(binary.binary_delay().to_value(u.s)-base_pint,np.longdouble)
        q=apply_binary(PulsarParameters(str(SEARCH_PAR)),mean+delta)
        ld=np.concatenate([lalpulsar.HeterodynedPulsarGetBSBDelay(q.PulsarParameters(),s.gpstimes,s.ssbdelay,s._HeterodynedCWSimulator__edat).data.copy() for s in sims])-base_lal
        td=tempo[i]-base_tempo
        phases={'PINT':-2*np.pi*frequency*pd,'TEMPO2':2*np.pi*frequency*td,'LAL':2*np.pi*frequency*ld}
        losses={};spans={}
        for engine,phase in phases.items():
            phase-=np.sum(weights*phase)
            losses[engine]=float(max(0,1-abs(np.sum(weights*np.exp(1j*phase)))**2))
            spans[engine]=float(np.ptp(phase))
        row=dict(name=name,power_loss=losses,phase_span_rad=spans,
            pint_lal_delta_max_seconds=float(max(abs(pd+ld))),tempo2_lal_delta_max_seconds=float(max(abs(td-ld))),
            pint_tempo2_delta_max_seconds=float(max(abs(pd+td))))
        rows.append(row)
        if name.startswith('metric_'):metrics[name]=phases
        if i%128==0:print(i,len(cases),row,flush=True)
    convergence=[]
    for eps in [.1,1.]:
        for engine in ['PINT','TEMPO2','LAL']:
            deriv=np.array([(metrics[f'metric_{eps}_{j}_1'][engine]-metrics[f'metric_{eps}_{j}_-1'][engine])/(2*eps) for j in range(10)],float)
            g=(deriv*weights)@deriv.T
            convergence.append(dict(engine=engine,eps=eps,eigenvalues=np.linalg.eigvalsh(g).tolist()))
    summary={}
    for group,predicate in [('gaussian',lambda n:n.startswith('gaussian')),('corners',lambda n:n.startswith('corner')),('axes',lambda n:n.startswith('correlated_axis'))]:
        selected=[r for r in rows if predicate(r['name'])]
        if selected:summary[group]=dict(count=len(selected),max_power_loss={e:max(r['power_loss'][e] for r in selected) for e in ['PINT','TEMPO2','LAL']})
    result=dict(scope='Independent binary uncertainty propagation at shared barycentric epochs and Earth positions; no fitting and no independent barycentering.',
        search_par_sha256=sha(SEARCH_PAR),covariance_sha256=sha(COVARIANCE),python_source_sha256=sha(__file__),tempo2_driver_sha256=sha(source),
        tempo2_library_sha256=sha(ROOT/'.conda/lib/libtempo2.so'),samples=len(tb),seed=20261016,
        versions={key:metadata.version(key) for key in ['pint-pulsar','lalsuite','cwinpy','numpy']},
        sample_input_sha256=sha(wd/'samples.txt'),draw_input_sha256=sha(wd/'draws.txt'),
        parameter_order=list(KEYS),canonical_perturbations=[d.tolist() for _,d in cases],
        central_pint_lal_difference_seconds=central_error,central_tempo2_lal_difference_seconds=float(max(abs(base_tempo-base_lal))),
        max_pint_lal_delta_difference_seconds=max(r['pint_lal_delta_max_seconds'] for r in rows),
        max_tempo2_lal_delta_difference_seconds=max(r['tempo2_lal_delta_max_seconds'] for r in rows),
        summary=summary,metric_convergence=convergence,cases=rows)
    checks={
        'central_delays_agree_within_100ns':max(central_error,result['central_tempo2_lal_difference_seconds'])<1e-7,
        'perturbed_delay_changes_agree_within_0p5ns':max(result['max_pint_lal_delta_difference_seconds'],result['max_tempo2_lal_delta_difference_seconds'])<5e-10,
        'correlated_cases_loss_below_1e_minus6':all(max(group['max_power_loss'].values())<1e-6 for group in summary.values()),
        'broken_cancellation_is_detected':all(min(r['power_loss'].values())>.1 for r in rows if r['name'] in ['single_PB','single_OMDOT','broken_PB_OMDOT']),
        'all_metrics_unresolved':max(x['eigenvalues'][-1] for x in convergence)<1e-7,
        'metric_engines_agree_within_one_percent':all(max(x['eigenvalues'][-1] for x in convergence if x['eps']==eps)/min(x['eigenvalues'][-1] for x in convergence if x['eps']==eps)<1.01 for eps in [.1,1.]),
    }
    result['acceptance_checks']=checks
    args.out.parent.mkdir(exist_ok=True,parents=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['cases','metric_convergence','canonical_perturbations']},indent=2))
    if not all(checks.values()):raise RuntimeError(checks)

if __name__=='__main__':main()
