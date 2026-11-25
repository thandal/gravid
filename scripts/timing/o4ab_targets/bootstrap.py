"""Conditional parametric calibration with re-estimated noise parameters.
Tests frequentist coverage under the fitted model, not whether that model is true.
"""

import json

import numpy as np
from gls import Fit, load
from paths import NAMES, OUT
from scipy.stats import chi2

# Keep this runnable without importing run.py and executing its work.


def mapping(wd, meta, fit):
    rows = [
        s.split("\t")
        for s in (wd / "covariance.tsv").read_text().splitlines()
        if s.startswith("parameter\t")
    ]
    idx = [
        next(i for i, (k, c) in enumerate(meta) if k == int(r[2]) and c == 0)
        for r in rows
    ]
    pos = [list(fit.cols).index(i) for i in idx]
    fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
    C = np.load(wd / "robust_covariances.npz")["formal_QR"]
    sd = np.sqrt(C.diagonal())
    L = np.linalg.cholesky(C / np.outer(sd, sd)) * sd[:, None]
    d = np.fromfile(wd / "axes.bin").reshape(1 + 2 * len(rows), 8192)
    D = (d[1 : len(rows) + 1] - d[len(rows) + 1 :]) / 0.2
    D -= D.mean(axis=1)[:, None]
    _, v = np.linalg.eigh(D @ D.T)
    a = np.linalg.solve(L.T, v[:, -1])
    return pos, fac, a


rng = np.random.default_rng(2026101701)
result = {}
B = 200
for name in NAMES:
    wd, data, M, meta = load(name)
    fit = Fit(data, M, np.ptp(data[:, 0]))
    kind = "red"
    base = fit.optimize(kind)
    pos, fac, a = mapping(wd, meta, fit)
    noise_root = np.linalg.cholesky(base["K"])
    zs = []
    qs = []
    thetas = []
    for b in range(B):
        y = base["EFAC"] * rng.normal(size=len(data)) + fit.W @ noise_root @ rng.normal(
            size=len(fit.t)
        )
        fit.set_y(y)
        boot = fit.optimize(kind, start=base["theta"], quick=True)
        mu = boot["beta"][pos] * fac
        C = boot["cov"][np.ix_(pos, pos)] * np.outer(fac, fac)
        zs.append(float(a @ mu / np.sqrt(a @ C @ a)))
        sd = np.sqrt(C.diagonal())
        u = mu / sd
        qs.append(float(u @ np.linalg.solve(C / np.outer(sd, sd), u)))
        thetas.append(boot["theta"].tolist())
        if (b + 1) % 25 == 0:
            print(name, b + 1, flush=True)
    z = np.array(zs)
    q = np.array(qs)
    result[name] = {
        "replicates": B,
        "seed": 2026101701,
        "generating_model": "fitted full-data EFAC + epoch jitter + 30-mode power-law red noise",
        "noise_parameters_refitted": True,
        "coverage_68_dominant_phase_direction": float(np.mean(abs(z) < 1)),
        "coverage_95_dominant_phase_direction": float(np.mean(abs(z) < 1.95996398454)),
        "z_rms": float(np.sqrt(np.mean(z * z))),
        "joint_binary_95_coverage": float(np.mean(q < chi2.ppf(0.95, len(pos)))),
        "joint_binary_q95": float(np.quantile(q, 0.95)),
        "scope": "Conditional parametric coverage; not evidence against an unmodelled process; 200 replicates imply roughly 1.5 percentage-point binomial standard error near 95% coverage.",
    }
    np.savez(OUT / (name + "_bootstrap.npz"), z=z, q=q, theta=np.array(thetas))
    (OUT / "bootstrap_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(name, result[name], flush=True)
