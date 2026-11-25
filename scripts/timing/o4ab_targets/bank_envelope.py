"""Count-one certificates for explicit sampled uncertainty regions.
Uses phase range, so it bounds any nonnegative weights on the evaluated times.
A simulation-calibrated ellipsoid threshold changes the region, not the covariance.
"""

import json

import numpy as np
from gls import Fit
from paths import BASE, NAMES, OUT
from scipy.stats import binom, chi2

profiles = json.loads((OUT / "noise_profile_results.json").read_text())
boots = json.loads((OUT / "bootstrap_results.json").read_text())
results = {}
for name in NAMES:
    if name not in boots:
        raise RuntimeError("Wait for all bootstrap targets to finish")
    qq = np.sort(np.load(OUT / (name + "_bootstrap.npz"))["q"])
    rank = int(binom.ppf(0.95, len(qq), 0.95))
    critical = max(float(chi2.ppf(0.95, 9)), float(qq[min(rank, len(qq) - 1)]))
    source = BASE / (name + "_ELL1")
    reference = np.fromfile(source / "axes.bin").reshape(19, 8192)[0]
    f0 = float(
        next(
            s.split()[1]
            for s in (source / "native.par").read_text().splitlines()
            if s.startswith("F0 ")
        )
    )
    ph = 2 * np.pi * 2 * f0 * 1.001
    dirs = [source] + (
        sorted((OUT / "physical_grid").iterdir()) if name == "J1231-1411" else []
    )
    worst = None
    total = 0
    for wd in dirs:
        n, p, nc = map(int, (wd / "design.meta").read_text().splitlines()[0].split())
        data = np.fromfile(wd / "toas.bin").reshape(n, 4)
        M = np.fromfile(wd / "design.bin").reshape(n, p)
        meta = np.loadtxt(wd / "design.meta", skiprows=1, dtype=int)
        rows = [
            s.split("\t")
            for s in (wd / "covariance.tsv").read_text().splitlines()
            if s.startswith("parameter\t")
        ]
        fit = Fit(data, M, np.ptp(data[:, 0]))
        idx = [
            next(i for i, (a, c) in enumerate(meta) if a == int(r[2]) and c == 0)
            for r in rows
        ]
        pos = [list(fit.cols).index(i) for i in idx]
        fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
        C0 = fit.C0[np.ix_(pos, pos)] * np.outer(fac, fac)
        sd = np.sqrt(C0.diagonal())
        L = np.linalg.cholesky(C0 / np.outer(sd, sd)) * sd[:, None]
        Li = np.linalg.inv(L)
        d = np.fromfile(wd / "axes.bin").reshape(19, 8192)
        D = (d[1:10] - d[10:]) / 0.2
        D -= D.mean(axis=1)[:, None]
        J = ph * D
        delta = ph * (d[0] - reference)
        delta -= delta.mean()
        maximum = 0
        for row in profiles[name]["cases"]:
            # Inflate the *noise model's* radiometer EFAC to a predeclared sensitivity value.
            # All fitted EFAC values are <1.1; all retained noise profiles are included,
            # even when a physical-grid point would reject that profile.
            theta = np.array(
                [
                    np.log(max(1.1, row["EFAC"])),
                    np.log(row["jitter_us"]),
                    np.log(row["red_rms_us"]),
                    row["gamma"],
                ]
            )
            fit.span = np.ptp(data[:, 0]) * row["span_factor"]
            obj, b, C, K = fit.evaluate(theta, "red", True)
            Cb = C[np.ix_(pos, pos)] * np.outer(fac, fac)
            mu = b[pos] * fac
            W = Li @ Cb @ Li.T
            variance = np.sum(J * (W @ J), axis=0)
            center = delta + (mu @ Li.T) @ J
            # A common constant phase is irrelevant; symmetric enclosure suffices.
            radius = float(
                np.max(abs(center) + np.sqrt(critical * np.maximum(variance, 0)))
            )
            bound = float(np.sin(radius) ** 2) if radius < np.pi / 2 else 1.0
            total += 1
            maximum = max(maximum, bound)
            if worst is None or bound > worst["any_sampled_weight_loss_bound"]:
                worst = {
                    "physical_case": wd.name,
                    "theta": theta.tolist(),
                    "span_factor": row["span_factor"],
                    "any_sampled_weight_loss_bound": bound,
                    "phase_radius_rad": radius,
                    "critical_q": critical,
                }
                np.savez(
                    OUT / (name + "_bank_worst.npz"),
                    binary_cov=Cb,
                    binary_shift=mu,
                    theta=theta,
                    critical_q=critical,
                    physical_dir=wd.name,
                    phase_radius_rad=radius,
                )
        print(name, wd.name, "max any sampled weight loss %", 100 * maximum, flush=True)
    results[name] = {
        "critical_q": critical,
        "bootstrap_upper_quantile_order": rank + 1,
        "bootstrap_replicates": len(qq),
        "cases": total,
        "worst": worst,
        "one_template_sufficient_at_loss": {
            str(tol): bool(worst["any_sampled_weight_loss_bound"] <= tol)
            for tol in [0.03, 0.01, 0.005, 0.001]
        },
        "scope": "Linearized binary phase ellipsoids, calibrated conditionally on fitted red-noise model; union over recorded noise grid (99%-chi4 REML sensitivity cutoff), two Fourier cutoffs, and all 24 physical J1231 cases. Fixed EFAC at least 1.1 is an explicit stress model. Any weights on 8192 O4ab times; continuum requires separate resolution/native checks, not an exact analytic continuum proof.",
    }
    (OUT / "bank_envelope_results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(name, results[name], flush=True)
