"""Disjoint frequency-band fits with shared epoch-noise covariance accounted for."""

import json

import numpy as np
from gls import Fit, load, serialize
from paths import NAMES, OUT
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2

results = {}
for name in NAMES:
    wd, data, M, meta = load(name)
    full = Fit(data, M, np.ptp(data[:, 0]))
    ff = full.optimize("red")
    rows = [
        s.split("\t")
        for s in (wd / "covariance.tsv").read_text().splitlines()
        if s.startswith("parameter\t")
    ]
    idx = [
        next(i for i, (a, c) in enumerate(meta) if a == int(r[2]) and c == 0)
        for r in rows
    ]
    operators = []
    fits = []
    for mask in [data[:, 3] < 1300, data[:, 3] >= 1300]:
        f = Fit(data[mask], M[mask], full.span)
        fit = f.optimize("red")
        K = fit["K"]
        s2 = fit["EFAC"] ** 2
        L = np.linalg.cholesky(K)
        Z = L @ cho_solve(
            cho_factor(np.eye(len(K)) + L.T @ f.WTW @ L / s2, lower=True), L.T
        )
        CiA = (f.A - f.W @ (Z @ (f.W.T @ f.A)) / s2) / s2
        H = fit["cov"] @ CiA.T
        pos = [list(f.cols).index(i) for i in idx]
        fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / f.scale[pos]
        op = np.zeros((len(pos), len(data)))
        op[:, mask] = H[pos] * fac[:, None]
        operators.append(op)
        fits.append(serialize(fit))
    H = operators[0] - operators[1]
    mu = H @ (data[:, 1] / data[:, 2])
    Z = H @ full.W
    C = ff["EFAC"] ** 2 * H @ H.T + Z @ ff["K"] @ Z.T
    sd = np.sqrt(C.diagonal())
    u = mu / sd
    q = float(u @ np.linalg.solve(C / np.outer(sd, sd), u))
    formal = np.load(wd / "robust_covariances.npz")["formal_QR"]
    sd0 = np.sqrt(formal.diagonal())
    L = np.linalg.cholesky(formal / np.outer(sd0, sd0)) * sd0[:, None]
    d = np.fromfile(wd / "axes.bin").reshape(19, 8192)
    D = (d[1:10] - d[10:]) / 0.2
    D -= D.mean(axis=1)[:, None]
    shift = np.linalg.solve(L, mu) @ D
    f0 = float(
        next(
            s.split()[1]
            for s in (wd / "native.par").read_text().splitlines()
            if s.startswith("F0 ")
        )
    )
    loss = float(np.mean((2 * np.pi * 2 * f0 * 1.001 * shift) ** 2))
    results[name] = {
        "bands_MHz": ["below 1300", "1300 and above"],
        "fits": fits,
        "binary_difference_chi2": q,
        "dof": len(mu),
        "nominal_p": float(chi2.sf(q, len(mu))),
        "central_binary_difference_linear_loss": loss,
        "shared_achromatic_noise_included": True,
        "scope": "Local linearized frequency split; statistic conditions on fitted noise parameters. Covariance of the difference includes shared epoch/red noise, not an independence assumption.",
    }
    print(name, results[name], flush=True)
(OUT / "frequency_checks.json").write_text(json.dumps(results, indent=2) + "\n")
