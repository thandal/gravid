"""Fit common noise families, blocked prediction, and orbital covariance propagation."""

import json

import numpy as np
from gls import Fit, epoch_projection, load, serialize
from paths import NAMES, OUT
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2


def binary_map(wd, meta, fit):
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
    physical = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
    cov = np.load(wd / "robust_covariances.npz")["formal_QR"]
    sd = np.sqrt(cov.diagonal())
    L = np.linalg.cholesky(cov / np.outer(sd, sd)) * sd[:, None]
    d = np.fromfile(wd / "axes.bin").reshape(1 + 2 * len(rows), 8192)
    D = (d[1 : len(rows) + 1] - d[len(rows) + 1 :]) / 0.2
    D -= D.mean(axis=1)[:, None]
    Li = np.linalg.inv(L)
    der = Li.T @ D
    f0 = float(
        next(
            s.split()[1]
            for s in (wd / "native.par").read_text().splitlines()
            if s.startswith("F0 ")
        )
    )
    G = (2 * np.pi * 2 * f0 * 1.001) ** 2 * der @ der.T / der.shape[1]
    return rows, pos, physical, G, L


results = {}
for name in NAMES:
    wd, data, M, meta = load(name)
    span = np.ptp(data[:, 0])
    fit = Fit(data, M, span)
    rows, pos, physical, G, L = binary_map(wd, meta, fit)
    res = {"full": {}, "blocked_prediction": {}}
    for kind in ["white", "jitter", "red", "exp"]:
        f = fit.optimize(kind)
        cov = f["cov"][np.ix_(pos, pos)] * np.outer(physical, physical)
        mu = f["beta"][pos] * physical
        rec = serialize(f)
        rec.update(
            {
                "binary_mean_linear_loss": float(np.trace(G @ cov)),
                "binary_3sigma_metric_box_linear_loss": float(9 * np.trace(G @ cov)),
                "central_shift_linear_loss": float(mu @ G @ mu),
            }
        )
        res["full"][kind] = rec
        np.savez(
            OUT / (name + "_" + kind + ".npz"),
            theta=f["theta"],
            binary_cov=cov,
            binary_shift=mu,
            beta=f["beta"],
            cov=f["cov"],
            scale=fit.scale,
            cols=fit.cols,
        )
        print(name, kind, rec, flush=True)
    par = {
        a[0]: a[1]
        for s in (wd / "native.par").read_text().splitlines()
        if len(a := s.split()) > 1
    }
    # Four disjoint validation quarters; the last is forward extrapolation.
    days = np.unique(np.floor(data[:, 0]).astype(int))
    for quarter, testdays in enumerate(np.array_split(days, 4), 1):
        mask = np.isin(np.floor(data[:, 0]).astype(int), testdays)
        test = data[mask]
        E = epoch_projection(test)
        tr = Fit(data[~mask], M[~mask], span)
        bykind = {}
        for kind in ["white", "jitter", "red", "exp"]:
            f = tr.optimize(kind)
            r, V, t = tr.predictive(f, kind, test, M[mask], E)
            cf = cho_factor(V, lower=True)
            q = float(r @ cho_solve(cf, r))
            u = (t - float(par["TASC"])) / float(par["PB"])
            ell = 2 * np.pi * (u - 0.5 * float(par.get("PBDOT", 0)) * u * u)
            tau = (t - t.mean()) / 365.25
            X = np.column_stack(
                [np.ones(len(t)), tau, tau * tau, np.sin(ell), np.cos(ell)]
            )
            ViX = cho_solve(cf, X)
            C = np.linalg.inv(X.T @ ViX)
            b = C @ ViX.T @ r
            qo = float(b[-2:] @ np.linalg.solve(C[-2:, -2:], b[-2:]))
            rec = serialize(f)
            rec.update(
                {
                    "n_epochs": len(t),
                    "chi2": q,
                    "nominal_global_p": float(chi2.sf(q, len(t))),
                    "orbital_coefficients_us": b[-2:].tolist(),
                    "orbital_sigma_us": np.sqrt(C.diagonal()[-2:]).tolist(),
                    "orbital_chi2": qo,
                    "nominal_orbital_p": float(chi2.sf(qo, 2)),
                }
            )
            bykind[kind] = rec
            np.savez(
                OUT / f"{name}_q{quarter}_{kind}.npz",
                residual_us=r,
                cov_us2=V,
                t=t,
                orbital_beta=b[-2:],
                orbital_cov=C[-2:, -2:],
                theta=f["theta"],
            )
        res["blocked_prediction"][f"q{quarter}"] = bykind
        print(
            name,
            "quarter",
            quarter,
            {
                k: (v["chi2"] / v["n_epochs"], v["nominal_orbital_p"])
                for k, v in bykind.items()
            },
            flush=True,
        )
    results[name] = res
    (OUT / "noise_results.json").write_text(json.dumps(results, indent=2) + "\n")
