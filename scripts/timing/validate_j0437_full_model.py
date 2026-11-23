#!/usr/bin/env python3
"""Independent binary-delay and nonlinear detector-weighted coverage checks."""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gravid.j0437_model import SEARCH_PAR, COVARIANCE, KEYS, apply_binary, binary_mean, load_covariance, waveform_metric
from gravid.orbital_metric import covariance_sha256


def main():
    import astropy.units as u
    from astropy.constants import c
    import lal
    import lalpulsar
    import pint.logging
    from pint.models import get_model
    from gravid.j0437_campaign import load_chunk
    from cwinpy.parfile import PulsarParameters
    from cwinpy.signal import HeterodynedCWSimulator
    from gravid.cwinpy_adapter import _ensure_template_amplitude
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/j0437_de440_qr/model_validation.json")
    args = ap.parse_args()
    pint.logging.setup(level="ERROR")
    samples = {}
    for path in sorted((ROOT / "results/campaign_het_v2").glob("het_*_g*.hdf5")):
        detector, times, data, _ = load_chunk(path)
        indices = np.linspace(0, len(times)-1, min(128, len(times)), dtype=int)
        samples.setdefault(detector, []).append((times[indices],
            np.full(len(indices), len(times)/len(indices)/np.var(data))))
    sampled = [(det, np.concatenate([x[0] for x in blocks]),
                np.concatenate([x[1] for x in blocks])) for det, blocks in sorted(samples.items())]
    cov, _ = load_covariance()
    p = _ensure_template_amplitude(PulsarParameters(str(SEARCH_PAR)))
    mean, factor = binary_mean(p), np.linalg.cholesky(cov)
    checks = []
    for stride, eps in ((4, .01), (4, .1), (1, .1), (1, 1.0)):
        sample = [(det, t[::stride], w[::stride]*stride) for det, t, w in sampled]
        metric = waveform_metric(SEARCH_PAR, cov, sample, eps=eps)
        checks.append(dict(stride=stride, eps=eps, eigenvalues=np.linalg.eigvalsh(metric).tolist()))
    # Use the actual placement metric's axes for all corners of its rotated box.
    metric = waveform_metric(SEARCH_PAR, cov, sampled, eps=.1)
    eigenvalues, rotation = np.linalg.eigh(metric)
    rng = np.random.default_rng(20261015)
    zdraws = rng.standard_normal((256, 10))
    corners = np.array(list(itertools.product((-3., 3.), repeat=10))) @ rotation.T
    draws = np.vstack([zdraws, corners])
    sims, weights = [], []
    for det, t, w in sampled:
        sim = HeterodynedCWSimulator(p, det, t[::8])
        sims.append(sim)
        weights.append(w[::8] * abs(sim.model(newpar=p))**2)
    weight = np.concatenate(weights); weight /= weight.sum()
    losses = []
    for i, z in enumerate(draws):
        q = apply_binary(PulsarParameters(str(SEARCH_PAR)), mean + factor @ z)
        phases = [sim.model(newpar=q, updateSSB=False, updateBSB=True,
                            phase_only=True, freqfactor=2.0) for sim in sims]
        overlap = abs(np.sum(weight * np.exp(-2j*np.pi*np.concatenate(phases))))**2
        losses.append(max(0., 1-overlap))
        if (i+1) % 256 == 0:
            print(f"nonlinear coverage {i+1}/{len(draws)}", flush=True)
    # Independent PINT DDK delays at the same barycentric coordinates; PINT
    # automatically converts TEMPO2's KIN/KOM convention for its DDK model.
    times = np.linspace(min(t[0] for _, t, _ in sampled), max(t[-1] for _, t, _ in sampled), 4096)
    sim = HeterodynedCWSimulator(p, "H1", times)
    model = get_model(str(SEARCH_PAR), allow_T2=True)
    binary = model.components["BinaryDDK"]
    binary.update_binary_object(None)
    tb = np.longdouble(model.T0.value) + (times.astype(np.longdouble)
          + sim.ssbdelay.data.astype(np.longdouble) - np.longdouble(p["T0"]))/86400
    positions = []
    for t in times:
        earth = lalpulsar.EarthState()
        lalpulsar.BarycenterEarth(earth, lal.LIGOTimeGPS(float(t)), sim._HeterodynedCWSimulator__edat)
        positions.append(earth.posNow.copy())
    binary.binary_instance.update_input(
        barycentric_toa=tb*u.day, obs_pos=np.array(positions)*u.s*c,
        psr_pos=model.ssb_to_psb_xyz_ICRS(epoch=tb.astype(float)))
    residual = np.asarray(binary.binary_instance.binary_delay().to_value(u.s) + sim.bsbdelay.data, float)
    old = PulsarParameters(str(ROOT / "results/campaign_het_v2/base_cwinpy.par"))
    oldphase = sim.model(newpar=old, updateSSB=True, updateBSB=True, phase_only=True)
    result = dict(
        search_par_sha256=covariance_sha256(SEARCH_PAR), covariance_sha256=covariance_sha256(COVARIANCE),
        parameter_order=list(KEYS), PB_OMDOT_correlation=float(cov[0,6]/np.sqrt(cov[0,0]*cov[6,6])),
        metric_convergence=checks, metric=metric.tolist(), largest_eigenvalue=float(eigenvalues[-1]),
        nonlinear=dict(gaussian_draws=256, rotated_box_corners=1024,
                       max_draw_loss=max(losses[:256]), max_corner_loss=max(losses[256:]), seed=20261015),
        independent_binary_delay=dict(implementation="PINT DDK versus LAL T2; shared barycentric times and Earth positions",
                                      max_absolute_difference_seconds=float(max(abs(residual))),
                                      difference_span_seconds=float(np.ptp(residual))),
        archived_ELL1_vs_refit_T2_uniform_time_power_overlap=float(abs(np.mean(np.exp(-2j*np.pi*oldphase)))**2),
    )
    if max(losses) >= 1e-4 or max(abs(residual)) >= 1e-7 or eigenvalues[-1] >= 1e-6:
        raise RuntimeError(result)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("metric", "metric_convergence")}, indent=2))


if __name__ == "__main__":
    main()
