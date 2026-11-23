#!/usr/bin/env python3
"""F-statistic for PSR J0437-4715 on the reduced O4a strain.

``scripts/j0437_upper_limit.py`` filters with one template of fixed
orientation (cos iota = 0, polarisation angle 0).  A signal whose phase
follows the timing model is, for any inclination, polarisation angle,
amplitude and phase, a complex combination of two templates: the detector
responses at polarisation angles 0 and pi/4.  Maximising over the two complex
coefficients gives the F-statistic (Jaranowski, Krolak & Schutz 1998), here
written for heterodyned data:

    2F = 2 Y^H M^-1 Y,   Y_k = sum_t w m_k^* d,   M_kl = sum_t w m_k^* m_l,

which follows a chi-squared distribution with four degrees of freedom in
Gaussian noise.

Same data, weights (``--noise-weights running-median``), orbit, control
frequencies and injected population as the fixed-orientation search, so the
two can be compared directly; the fixed-orientation statistic is recomputed
here from the same sums as a cross-check.

Memory: one detector per process (about 2.7 GB), then a light combine step.

    python scripts/j0437_f_statistic.py --detector H1 --out-dir DIR
    python scripts/j0437_f_statistic.py --detector L1 --out-dir DIR
    python scripts/j0437_f_statistic.py --combine --out-dir DIR
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PSIS = (0.0, np.pi / 4)
SEED = 20260707                      # as in j0437_upper_limit.py


def _c(z):
    return [float(np.real(z)), float(np.imag(z))]


def _z(pair):
    return complex(pair[0], pair[1])


def _sums(ma, mb, weights, series):
    """Y (2,) and M (2,2) for one template pair against a series."""
    wa = weights * ma
    ya, maa, mab = np.vdot(wa, series), np.vdot(wa, ma).real, np.vdot(wa, mb)
    del wa
    wb = weights * mb
    yb, mbb = np.vdot(wb, series), np.vdot(wb, mb).real
    del wb
    return np.array([ya, yb]), np.array([[maa, mab], [np.conj(mab), mbb]])


def two_f(y, m):
    return float(2.0 * np.real(np.conj(y) @ np.linalg.solve(m, y)))


def fixed(y, m):
    return float(2.0 * abs(y[0]) ** 2 / np.real(m[0, 0]))


def draws(post, n_inj):
    rng = np.random.default_rng(SEED)
    orbits, _ = post.sample_mc(n_inj, rng)
    orientation = [(rng.uniform(-1, 1), rng.uniform(0, np.pi / 2),
                    rng.uniform(0, 2 * np.pi)) for _ in range(n_inj)]
    return np.asarray(orbits), orientation


def run_detector(args):
    import resource
    from scripts import j0437_upper_limit as ul
    from gravid.posterior import OrbitalPosterior
    from cwinpy.parfile import PulsarParameters

    ul.NOISE_WEIGHTS = "running-median"
    files = sorted(glob.glob(os.path.join(ul.HETDIR, "het_*_g*.hdf5")))
    if args.limit_files:                       # smoke test only
        files = [f for det in ("H1", "L1") for f in
                 [g for g in files if f"_{det}_" in os.path.basename(g)][:args.limit_files]]
    identity = ul.verify_selection(files, ul.BASE_PAR, smoke=bool(args.limit_files))
    mine = [f for f in files if f"_{args.detector}_" in os.path.basename(f)]
    demod = ul.MultiDatasetDemodulator(mine, ul.BASE_PAR)
    noise_record = ul._noise_record(demod)
    (chunk,) = demod.chunks
    weights, data = chunk["invvar"], chunk["data"]
    covariance, _ = ul.load_covariance(ul.COVARIANCE_ARTIFACT)
    post = OrbitalPosterior(demod.central_theta, covariance, names=ul.KEYS)
    theta0 = np.asarray(post.mean, float)

    def template(offset, psi):
        par = demod._par_for(theta0, offset)
        par["H0"] = 1.0; par["COSIOTA"] = 0.0
        par["PSI"] = float(psi); par["PHI0"] = 0.0
        return demod._model(chunk, par)

    def pair(offset):
        return [template(offset, psi) for psi in PSIS]

    def peak():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024

    t0 = time.time()
    ma = template(0.0, PSIS[0])
    if args.scan:
        # Dense scan of the fixed-orientation sum against a pure frequency
        # shift, by FFT on the 1-s grid (single precision; a diagnostic).
        from scipy.fft import fft, next_fast_len
        start, stop = args.scan_window
        # Pad to a length with small prime factors: an awkward length makes
        # the transform fall back to a far slower, memory-hungry algorithm.
        n = next_fast_len(int(round(stop - start)) + 1, real=False)
        z = np.zeros(n, np.complex64)
        index = np.round(np.asarray(chunk["times"], float) - start).astype(np.int64)
        product = np.conj(ma)
        product *= data
        product *= weights
        z[index] = product
        del index, product
        spectrum = fft(z, overwrite_x=True)
        del z
        kmax = int(0.3 * n)
        band = np.concatenate([spectrum[n - kmax:], spectrum[:kmax + 1]])
        del spectrum
        np.savez(args.scan, band=band, n_fft=n, start=start, kmax=kmax,
                 norm=float(np.vdot(ma, weights * ma).real))
        del band
        print(f"[{args.detector}] scan saved ({time.time()-t0:.0f}s, "
              f"peak {peak()} MB)", flush=True)
    mb = template(0.0, PSIS[1])
    y_data, m0 = _sums(ma, mb, weights, data)
    weighted0 = [weights * ma, weights * mb]
    del ma, mb
    print(f"[{args.detector}] target sums done ({time.time()-t0:.0f}s, "
          f"peak {peak()} MB)", flush=True)

    offsets = ul._offsets(args.n_offsets)
    control_y, control_m = [], []
    for i, offset in enumerate(offsets):
        ta, tb = pair(offset)
        y, m = _sums(ta, tb, weights, data)
        del ta, tb
        control_y.append(y); control_m.append(m)
        if i % 10 == 9:
            print(f"[{args.detector}] controls {i+1}/{len(offsets)} "
                  f"({time.time()-t0:.0f}s, peak {peak()} MB)", flush=True)

    orbits, orientation = draws(post, args.n_inj)
    y_signal = np.zeros((args.n_inj, 2), complex)
    snr2 = np.zeros(args.n_inj)
    for j, (cosi, psi, phi0) in enumerate(orientation):
        par = demod._apply_orbit(PulsarParameters(demod.search_par), orbits[j])
        par["H0"] = 1.0; par["COSIOTA"] = float(cosi)
        par["PSI"] = float(psi); par["PHI0"] = float(phi0)
        s = demod._model(chunk, par)
        y_signal[j] = [np.vdot(weighted0[0], s), np.vdot(weighted0[1], s)]
        snr2[j] = float(np.vdot(s, weights * s).real)
        del s
        if j % 25 == 24:
            print(f"[{args.detector}] injections {j+1}/{args.n_inj} "
                  f"({time.time()-t0:.0f}s, peak {peak()} MB)", flush=True)

    out = dict(
        detector=args.detector, input_identity=identity,
        noise_weights=noise_record, n_samp=demod.n_samp,
        search_par_sha256=ul.covariance_sha256(ul.SEARCH_PAR),
        covariance_sha256=ul.covariance_sha256(ul.COVARIANCE_ARTIFACT),
        polarisation_angles=list(PSIS), offsets=list(offsets),
        population_seed=SEED, orbital_draws=orbits.tolist(),
        orientation_draws=[list(d) for d in orientation],
        y_data=[_c(v) for v in y_data], m=[[_c(v) for v in row] for row in m0],
        control_y=[[_c(v) for v in y] for y in control_y],
        control_m=[[[_c(v) for v in row] for row in m] for m in control_m],
        y_signal=[[_c(v) for v in y] for y in y_signal], snr2=snr2.tolist(),
        peak_memory_mb=int(peak()),
        invocation=dict(argv=sys.argv[1:],
                        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S")))
    path = os.path.join(args.out_dir, f"f_statistic.{args.detector}.json")
    with open(path + ".part", "w") as handle:
        json.dump(out, handle, indent=1)
    os.replace(path + ".part", path)
    print(f"[{args.detector}] wrote {path} ({time.time()-t0:.0f}s, "
          f"peak {peak()} MB)", flush=True)


def combine(args):
    from scipy.stats import chi2
    from gravid.upper_limit import (heterodyne_statistic_response,
                                      linear_grid_crossing,
                                      paired_grid_bootstrap_crossing)

    parts = [json.load(open(os.path.join(args.out_dir, f"f_statistic.{det}.json")))
             for det in ("H1", "L1")]
    for key in ("offsets", "orbital_draws", "orientation_draws",
                "search_par_sha256", "covariance_sha256", "input_identity"):
        assert parts[0][key] == parts[1][key], key

    def total(key, shape):
        return sum(np.array([_z(v) for v in np.reshape(p[key], (-1, 2))]).reshape(shape)
                   for p in parts)

    offsets = np.array(parts[0]["offsets"])
    n_off, n_inj = len(offsets), len(parts[0]["snr2"])
    y_data, m0 = total("y_data", (2,)), total("m", (2, 2))
    control_y, control_m = total("control_y", (n_off, 2)), total("control_m", (n_off, 2, 2))
    y_signal = total("y_signal", (n_inj, 2))
    snr2 = sum(np.array(p["snr2"]) for p in parts)
    response = np.atleast_1d(heterodyne_statistic_response(offsets))
    controls = np.array([two_f(y, m) for y, m in zip(control_y, control_m)]) / response
    controls_fixed = np.array([fixed(y, m) for y, m in zip(control_y, control_m)]) / response
    captured = np.array([np.real(np.conj(y) @ np.linalg.solve(m0, y))
                         for y in y_signal]) / snr2
    scale = np.sqrt(8.0 / np.median(snr2))
    h0_grid = np.geomspace(0.25 * scale, 12.0 * scale, 49)

    def recovery(statistic, threshold):
        detected = np.array([[statistic(y_data + h0 * ys, m0) >= threshold
                              for ys in y_signal] for h0 in h0_grid])
        efficiency = detected.mean(axis=1)
        boot = paired_grid_bootstrap_crossing(
            h0_grid, detected, target=0.95, confidence=0.90,
            n_bootstrap=20000, seed=20260825)
        return dict(threshold=float(threshold),
                    h0_95=linear_grid_crossing(h0_grid, efficiency, 0.95),
                    h0_95_interval=boot["interval"],
                    efficiency=efficiency.tolist())

    out = dict(
        description=__doc__.strip().splitlines()[0],
        n_samp=sum(p["n_samp"] for p in parts),
        input_identity=parts[0]["input_identity"],
        noise_weights={p["detector"]: p["noise_weights"] for p in parts},
        search_par_sha256=parts[0]["search_par_sha256"],
        covariance_sha256=parts[0]["covariance_sha256"],
        polarisation_angles=list(PSIS),
        target=dict(two_f=two_f(y_data, m0), fixed_orientation=fixed(y_data, m0),
                    two_f_p_chi2_4=float(chi2.sf(two_f(y_data, m0), 4)),
                    controls_at_or_above=int(np.sum(controls >= two_f(y_data, m0)))),
        controls=dict(offsets=offsets.tolist(), response=response.tolist(),
                      two_f=controls.tolist(),
                      fixed_orientation=controls_fixed.tolist()),
        injection=dict(
            n_inj=n_inj, population_seed=SEED, h0_grid=h0_grid.tolist(),
            snr2=snr2.tolist(), captured_power_fraction=captured.tolist(),
            f_statistic=dict(
                loudest_control=recovery(two_f, controls.max()),
                false_alarm_1pct=recovery(two_f, chi2.isf(0.01, 4))),
            fixed_orientation=dict(
                loudest_control=recovery(fixed, controls_fixed.max()),
                false_alarm_1pct=recovery(fixed, chi2.isf(0.01, 2)))),
        products=dict(
            y_data=[_c(v) for v in y_data], m=[[_c(v) for v in r] for r in m0],
            control_y=[[_c(v) for v in y] for y in control_y],
            control_m=[[[_c(v) for v in r] for r in m] for m in control_m],
            y_signal=[[_c(v) for v in y] for y in y_signal]))
    path = os.path.join(args.out_dir, "f_statistic.json")
    with open(path + ".part", "w") as handle:
        json.dump(out, handle, indent=1, default=float)
    os.replace(path + ".part", path)
    print(f"target: 2F = {out['target']['two_f']:.3f} "
          f"(chi-squared(4) probability {out['target']['two_f_p_chi2_4']:.3f}; "
          f"{out['target']['controls_at_or_above']} of {n_off} controls at or above), "
          f"fixed orientation {out['target']['fixed_orientation']:.3f}")
    print(f"controls: 2F mean {controls.mean():.2f} (4 expected), max "
          f"{controls.max():.2f}; fixed orientation mean {controls_fixed.mean():.2f} "
          f"(2 expected), max {controls_fixed.max():.2f}")
    print(f"signal power inside the two-template space: min {captured.min():.7f}")
    for name in ("f_statistic", "fixed_orientation"):
        for rule, r in out["injection"][name].items():
            h0 = "not bracketed" if r["h0_95"] is None else f"{r['h0_95']:.3e}"
            print(f"{name} / {rule}: threshold {r['threshold']:.2f} -> "
                  f"h0_95 = {h0}  {r['h0_95_interval']}")
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--detector", choices=["H1", "L1"])
    parser.add_argument("--combine", action="store_true")
    parser.add_argument("--n-offsets", type=int, default=40)
    parser.add_argument("--n-inj", type=int, default=200)
    parser.add_argument("--limit-files", type=int, default=0,
                        help="smoke test: first N files per detector")
    parser.add_argument("--scan", default=None,
                        help="optional .npy path for a dense frequency scan")
    parser.add_argument("--scan-window", type=float, nargs=2,
                        default=[1368195280.5, 1389455057.5])
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.combine:
        combine(args)
    elif args.detector:
        run_detector(args)
    else:
        parser.error("give --detector or --combine")


if __name__ == "__main__":
    main()
