"""Table 1: the article's test on later radio observations, without the
fitted covariance.

With ``--native-rerun`` the script starts from the downloaded releases: it
builds the 14 timing models (seven per release), evaluates their residuals
with TEMPO2 (no fit is run), and computes the table.  It needs:

  * the MeerKAT and InPTA files fetched by
    ``scripts/timing/orbital_offset_fit/acquire.py``;
  * the published PPTA parameter file, saved as
    ``data/J0437-4715_reardon2024_dr3.par``;
  * TEMPO2 in a conda environment at ``.conda/`` in the repository root.

Each model keeps the release's own propagation and instrument terms
(dispersion, frequency-dependent terms, jumps) and takes spin, astrometry and
orbit from the model under test.  Those retained terms are NOT independent
predictions.

Without the flag the script reads residuals saved by an earlier run from
``results/timing/first_principles/``.  Where such files exist, a rerun is
also compared against them.
"""
from decimal import Decimal
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

# Parameters taken from the model under test; everything else in a release's
# parameter file is kept.
SWAP = set('RAJ DECJ PMRA PMDEC PX POSEPOCH F0 F1 F2 F3 PEPOCH BINARY PB PBDOT A1 XDOT T0 OM ECC OMDOT M2 KIN KOM SINI '
           'TASC EPS1 EPS2 EPHEM UNITS PLANET_SHAPIRO'.split())
RELEASE_PAR = {'meerkat': 'meerkat/partim/J0437-4715.par', 'inpta': 'inpta/J0437-4715/J0437-4715.DMX.par'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def lines_of(path):
    return {l.split()[0]: l for l in path.read_text().splitlines() if l.split() and not l.startswith('#')}

def build_models(directory, published_par, search_par, release_pars):
    """Write the 14 timing models of Table 1 into ``directory``.

    ``release_pars`` maps 'meerkat' and 'inpta' to each release's parameter
    file.  The four shifted models start from the published parameters and
    use their quoted one-sigma errors; the partner shifts follow the two
    longitude relations in the article.
    """
    pub = lines_of(published_par)
    def val(key): return Decimal(pub[key].split()[1])
    def sig(key): return Decimal(pub[key].split()[3])
    PBd, sPB, sT0 = val('PB'), sig('PB'), sig('T0')
    dOMDOT = 2*Decimal(math.pi)*sPB/PBd**2 * Decimal(180/math.pi) * Decimal('365.25')   # deg/yr
    dOM = 2*Decimal(math.pi)*sT0/PBd * Decimal(180/math.pi)                                # deg
    def perturbed(changes):
        out = dict(pub)
        for k, d in changes.items():
            f = out[k].split(); f[1] = str(Decimal(f[1]) + d); out[k] = ' '.join(f)
        return out
    models = {
        'S_{site}_published_DR3': pub,
        'S_{site}_search_template_model_(TDB)': lines_of(search_par),
        '{site}_control_PB_p1_sigma': perturbed({'PB': sPB}),
        '{site}_control_PB_p1s_OMDOT_partner': perturbed({'PB': sPB, 'OMDOT': dOMDOT}),
        '{site}_control_T0_p1_sigma': perturbed({'T0': sT0}),
        '{site}_control_T0_p1s_OM_partner': perturbed({'T0': sT0, 'OM': dOM}),
    }
    directory.mkdir(parents=True, exist_ok=True)
    for site, release_par in release_pars.items():
        release = [l for l in Path(release_par).read_text().splitlines() if l.split()]
        kept = [l for l in release if l.split()[0] not in SWAP]
        for pattern, model in models.items():
            taken = [model[k] for k in model if k in SWAP and k != 'PLANET_SHAPIRO']
            (directory/(pattern.format(site=site)+'.par')).write_text('\n'.join(kept + taken + ['PLANET_SHAPIRO N']) + '\n')
        own = [l for l in release if l.split()[0] != 'PLANET_SHAPIRO']
        (directory/f'F_{site}_release_own_fit_(uses_O4a_data).par').write_text('\n'.join(own + ['PLANET_SHAPIRO N']) + '\n')

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
    if not args.native_rerun and not SRC.is_dir():
        parser.error(f'no saved residuals in {SRC.relative_to(ROOT)}; run with --native-rerun')
    if args.native_rerun:
        args.native_dir.mkdir(parents=True, exist_ok=True)
        build_models(args.native_dir, par, ROOT/'gravid/resources/j0437_de440_search.par',
                     {site: EXT/name for site, name in RELEASE_PAR.items()})
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
            saved = SRC/(stem+'.res')
            if args.native_rerun:
                p = args.native_dir/(stem+'.par')
                manifest['model built for this run: '+stem+'.par'] = sha(p)
                if (SRC/(stem+'.par')).is_file():
                    assert p.read_bytes() == (SRC/(stem+'.par')).read_bytes(), stem
            else:
                p = SRC/(stem+'.par')
                manifest[str(p.relative_to(ROOT))] = sha(p)
            if saved.is_file():
                manifest[str(saved.relative_to(ROOT))] = sha(saved)
                data = np.loadtxt(saved)
            if args.native_rerun:
                fresh = args.native_dir/(stem+'.res')
                completed = subprocess.run([str(helper), str(p), str(tim), str(EXT/'clock'), str(fresh)],
                    env=dict(os.environ, TEMPO2=str(ROOT/'.conda/share/tempo2')), capture_output=True, text=True, check=True)
                fresh.with_suffix('.log').write_text(completed.stdout+completed.stderr)
                new = np.loadtxt(fresh)
                if saved.is_file():
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
