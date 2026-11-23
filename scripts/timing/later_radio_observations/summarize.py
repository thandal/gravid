"""Reassemble the article's test on later radio observations, without the
fitted covariance.

Optionally rerun the native TEMPO2 residual evaluator, which has no fit call.
The radio releases' propagation/instrument terms are retained in the saved
hybrid par files. They are NOT independent predictions of those nuisance terms.
"""
from pathlib import Path
import argparse, hashlib, json, math, os, subprocess
import numpy as np
from scipy.optimize import brentq
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / 'results/timing/first_principles'
EXT = ROOT / 'results/timing/j0437_external_validation'
OUT = Path(__file__).resolve().parent / 'results'
LABELS = {
    'published_parameters': 'S_{site}_published_DR3',
    'search_model': 'S_{site}_search_template_model_(TDB)',
    'release_fit': 'F_{site}_release_own_fit_(uses_O4a_data)',
    'PB_only': '{site}_control_PB_p1_sigma',
    'PB_with_OMDOT': '{site}_control_PB_p1s_OMDOT_partner',
    'T0_only': '{site}_control_T0_p1_sigma',
    'T0_with_OM': '{site}_control_T0_p1s_OM_partner',
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-rerun', action='store_true')
    parser.add_argument('--native-dir', type=Path, default=Path('/tmp/gravid_later_radio_observations'))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    values, errors = {}, {}
    par = ROOT / 'data/J0437-4715_reardon2024_dr3.par'
    for line in par.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] not in ('PB', 'T0', 'OM', 'OMDOT', 'F0', 'A1', 'ECC'):
            continue
        values[fields[0]] = float(fields[1])
        if len(fields) >= 4:
            errors[fields[0]] = float(fields[3])
    pb = values['PB']; fgw = 2 * values['F0']
    ratios = {
        'PB_over_OMDOT_partner': errors['PB'] / (pb**2 * math.radians(errors['OMDOT']) / 365.25 / (2*math.pi)),
        'T0_over_OM_partner': errors['T0'] / (math.radians(errors['OM']) * pb / (2*math.pi)),
    }
    limits = {}
    for loss in (.03, .20):
        amp = brentq(lambda a: 1-j0(a)**2-loss, 0, 2.4)
        limits[str(loss)] = {'GW_phase_amplitude_rad': amp, 'sinusoidal_delay_amplitude_ms': amp/(2*math.pi*fgw)*1000}
    manifest = {str(par.relative_to(ROOT)): sha(par)}
    native_errors = {}
    if args.native_rerun:
        args.native_dir.mkdir(parents=True, exist_ok=True)
        helper = args.native_dir / 'evaluate'
        cpp = ROOT / 'scripts/timing/orbital_offset_fit/evaluate_timing.C'
        subprocess.run(['g++', '-O2', '-I'+str(ROOT/'.conda/include'), str(cpp),
                        '-L'+str(ROOT/'.conda/lib'), '-Wl,-rpath,'+str(ROOT/'.conda/lib'),
                        '-ltempo2', '-o', str(helper)], check=True)
        manifest[str(cpp.relative_to(ROOT))] = sha(cpp)
        manifest[str((ROOT/'.conda/lib/libtempo2.so').relative_to(ROOT))] = sha(ROOT/'.conda/lib/libtempo2.so')
    for clock in sorted((EXT/'clock').glob('*.clk')):
        manifest[str(clock.relative_to(ROOT))] = sha(clock)
    results = {}
    for site in ('meerkat', 'inpta'):
        tim = EXT/site/'flattened.tim'
        manifest[str(tim.relative_to(ROOT))] = sha(tim)
        ids = np.array([line.split()[0] for line in tim.read_text().splitlines()
                        if len(line.split())>4 and line.split()[0] not in ('FORMAT', 'MODE')])
        results[site] = {}
        for label, pattern in LABELS.items():
            stem = pattern.format(site=site)
            p = SRC/(stem+'.par'); saved = SRC/(stem+'.res')
            manifest[str(p.relative_to(ROOT))] = sha(p)
            manifest[str(saved.relative_to(ROOT))] = sha(saved)
            data = np.loadtxt(saved)
            if args.native_rerun:
                fresh = args.native_dir/(stem+'.res')
                completed = subprocess.run([str(helper), str(p), str(tim), str(EXT/'clock'), str(fresh)],
                    env=dict(os.environ, TEMPO2=str(ROOT/'.conda/share/tempo2')), capture_output=True, text=True, check=True)
                fresh.with_suffix('.log').write_text(completed.stdout+completed.stderr)
                new = np.loadtxt(fresh)
                assert new.shape == data.shape
                assert np.array_equal(new[:,0], data[:,0])
                # Compare after removing the unobservable common phase origin.
                delta = new[:,2]-data[:,2]; delta -= delta.mean()
                err = float(np.max(abs(delta)))
                assert err < 1e-9, (site, label, err)
                native_errors[site+'/'+label] = err
                data = new
            assert len(data) == len(ids)
            mask = (data[:,0] >= 60088.625) & (data[:,0] <= 60325.65554976852)
            d = data[mask]; names = ids[mask]
            z = np.exp(2j*np.pi*fgw*d[:,2]); weights = d[:,3]**-2
            groups = list(dict.fromkeys(names))
            group_z = np.array([np.average(z[names==group], weights=weights[names==group]) for group in groups])
            residual = np.angle(group_z/group_z.mean())/(2*np.pi*fgw)
            row = {
                'TOAs': len(d), 'file_groups': len(groups),
                'observing_days': len(np.unique(np.floor(d[:,0]))),
                'MJD_range': [float(d[:,0].min()), float(d[:,0].max())],
                'residual_RMS_ns': float(np.std(residual)*1e9),
                'fractional_power_loss': float(1-abs(group_z.mean())**2),
            }
            results[site][label] = row
    for site in results:
        r = results[site]
        assert r['published_parameters']['fractional_power_loss'] < 1e-4
        assert r['search_model']['fractional_power_loss'] < 1e-4
        assert min(r[key]['fractional_power_loss'] for key in ('PB_only','T0_only')) > .7
        assert max(r[key]['fractional_power_loss'] for key in ('PB_with_OMDOT','T0_with_OM')) < 1e-4
    # Cross-check the independent earlier aggregator, including official-run cuts.
    audit = ROOT/'scripts/timing/supporting_checks/external_prediction_audit.json'
    older = json.loads(audit.read_text())
    mapping = {'published_parameters':'published_DR3','PB_only':'control_PB_p1_sigma',
               'PB_with_OMDOT':'control_PB_p1s_OMDOT_partner','T0_only':'control_T0_p1_sigma','T0_with_OM':'control_T0_p1s_OM_partner'}
    for site in results:
        for newkey, oldkey in mapping.items():
            expected = older[site][oldkey]['official_O4a']
            assert abs(results[site][newkey]['fractional_power_loss']-expected['coherence_loss']) < 1e-10
    manifest[str(audit.relative_to(ROOT))] = sha(audit)
    manifest[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    result = {
        'frequency_Hz': fgw, 'eccentricity':values['ECC'],
        'selection_MJD':[60088.625,60325.65554976852],
        'degeneracy_ratios':ratios,'sinusoidal_tolerances':limits,
        'predictions':results,'native_rerun':args.native_rerun,
        'native_max_centered_residual_difference_s':native_errors,
        'scope':'No spin, astrometric or orbital refit. Retains release propagation/instrument corrections, some estimated with O4a data. Inverse-variance channel averaging then equal file-group weights. Sparse radio-sampled coherence, not detector-weighted mismatch or a confidence bound.',
        'input_sha256':manifest,
    }
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    table = ['| Timing parameters / control | MeerKAT power loss | InPTA power loss |',
             '|---|---:|---:|']
    names = ['Published PPTA parameters','DE440 search model','Release’s own fit (includes O4a)',
             '$P_b+1\\sigma$ alone','Same $P_b$ shift with $\\dot\\omega$ partner',
             '$T_0+1\\sigma$ alone','Same $T_0$ shift with $\\omega$ partner']
    for key, name in zip(LABELS,names):
        a,b = [100*results[site][key]['fractional_power_loss'] for site in ('meerkat','inpta')]
        fmt = lambda v: (f'{v:.2e}%' if v < .001 else f'{v:.3g}%')
        table.append(f'| {name} | {fmt(a)} | {fmt(b)} |')
    (OUT/'prediction_table.md').write_text('\n'.join(table)+'\n')
    print(json.dumps({'ratios':ratios,'tolerances':limits,'native_cases':len(native_errors),
                      'native_max_difference_s':max(native_errors.values(), default=None)},indent=2))
    print('\n'.join(table))

if __name__ == '__main__':
    main()
