"""Full training-noise/parameter covariance in the newly acquired O4b test."""

import json
import os
import subprocess
from pathlib import Path

import numpy as np
from gls import Fit, epoch_projection, load
from paths import INPUT, TEMPO, TOOLS, WORK
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2

OUT = WORK / "radio"
BASE = WORK / "fits"
ENV = dict(os.environ, TEMPO2=str((TEMPO / "share/tempo2")))
results = {}
for name in ["J1045-4509", "J1737-0811"]:
    wd, data, M, meta = load(name)
    testdir = OUT / "predictions" / name / "PTA/O4b"
    par = testdir / "earlier_ELL1.par"
    tim = testdir / "input.tim"
    path = testdir / "full_design.res"
    with (testdir / "full_design.log").open("w") as log:
        subprocess.run(
            [
                str(TOOLS / "predict_design"),
                str(par),
                str(tim),
                str((INPUT / "clock")),
                str(path),
            ],
            env=ENV,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    test = np.loadtxt(path)
    n, p = map(int, Path(str(path) + ".meta").read_text().splitlines()[0].split())
    tm = np.loadtxt(str(path) + ".meta", skiprows=1, dtype=int)
    T = np.fromfile(str(path) + ".design").reshape(n, p)
    A = np.zeros((n, len(meta)))
    lookup = {tuple(v): i for i, v in enumerate(tm)}
    for j, key in enumerate(meta):
        if tuple(key) in lookup:
            A[:, j] = T[:, lookup[tuple(key)]]
    td = test[:, [0, 2, 3, 4]]
    E = epoch_projection(td)
    fit = Fit(data, M, np.ptp(data[:, 0]))
    pdict = {
        a[0]: a[1] for s in par.read_text().splitlines() if len(a := s.split()) > 1
    }
    res = {}
    for kind in ["white", "jitter", "red", "exp"]:
        f = fit.optimize(kind)
        r, V, t = fit.predictive(f, kind, td, A, E)
        cf = cho_factor(V, lower=True)
        u = (t - float(pdict["TASC"])) / float(pdict["PB"])
        ell = 2 * np.pi * (u - 0.5 * float(pdict.get("PBDOT", 0)) * u * u)
        tau = (t - t.mean()) / 365.25
        X = np.column_stack([np.ones(len(t)), tau, tau * tau, np.sin(ell), np.cos(ell)])
        ViX = cho_solve(cf, X)
        C = np.linalg.inv(X.T @ ViX)
        b = C @ ViX.T @ r
        q = float(b[-2:] @ np.linalg.solve(C[-2:, -2:], b[-2:]))
        Xp = X[:, :3]
        Cp = np.linalg.inv(Xp.T @ cho_solve(cf, Xp))
        bp = Cp @ Xp.T @ cho_solve(cf, r)
        qp = float(r @ cho_solve(cf, r) - bp @ np.linalg.solve(Cp, bp))
        np.savez(testdir / (kind + "_prediction.npz"), residual_us=r, cov_us2=V, t=t)
        res[kind] = {
            "after_constant_and_spin_chi2": qp,
            "after_constant_and_spin_dof": len(r) - 3,
            "after_constant_and_spin_nominal_p": float(chi2.sf(qp, len(r) - 3)),
            "sin_us": float(b[-2]),
            "cos_us": float(b[-1]),
            "sigma_us": np.sqrt(C.diagonal()[-2:]).tolist(),
            "orbital_chi2": q,
            "nominal_p": float(chi2.sf(q, 2)),
            "global_chi2": float(r @ cho_solve(cf, r)),
            "n_epochs": len(r),
        }
    results[name] = res
    print(name, res, flush=True)
    (OUT / "external_predictive_results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
