"""Cross J1231's physical grid with secular ELL1 and fitted noise.
Full native refits; common reference; no mass-prior weights or posterior claim.
"""

import json
import os
import subprocess

import numpy as np
from gls import Fit, serialize
from paths import BASE, HERE, INPUT, OUT, TEMPO, TOOLS
from scipy.stats import chi2

ENV = dict(
    os.environ,
    TEMPO2=str((TEMPO / "share/tempo2")),
    TEMPO2_PLUG_PATH=str(TOOLS),
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
)
source = BASE / "J1231-1411_ELL1"
out = OUT / "physical_grid"
out.mkdir(exist_ok=True)
reference = np.fromfile(source / "axes.bin").reshape(19, 8192)[0]
f0 = float(
    next(
        s.split()[1]
        for s in (source / "native.par").read_text().splitlines()
        if s.startswith("F0 ")
    )
)
phasefac = 2 * np.pi * 2 * f0 * 1.001


def delays(wd, rows, draws, tag):
    path = wd / (tag + ".txt")
    with path.open("w") as f:
        f.write(
            f"8192 {len(draws)} {len(rows)} 60088.625 {60703 + 17 / 24:.16f}\n"
            + " ".join(r[2] for r in rows)
            + "\n"
        )
        np.savetxt(f, draws, fmt="%.19e")
    with (wd / (tag + ".log")).open("w") as log:
        subprocess.run(
            [
                str((TOOLS / "delays")),
                str(wd / "native.par"),
                str(path),
                str(wd / (tag + ".bin")),
                "unused",
            ],
            env=ENV,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return np.fromfile(wd / (tag + ".bin")).reshape(-1, 8192)


results = {}
grid = json.loads((HERE / "models/j1231_shapiro_grid.json").read_text())
for label, fix in grid.items():
    wd = out / label
    wd.mkdir(exist_ok=True)
    ls = []
    for line in (source / "native.par").read_text().splitlines():
        a = line.split()
        if not a or a[0] in ["START", "FINISH", "NTOA", "TRES", "CHI2R", "NITS"]:
            continue
        if a[0] in fix:
            line = f"{a[0]} {fix[a[0]]:.19g} 0"
        ls.append(line)
    (wd / "input.par").write_text("\n".join(ls + ["NITS 5"]) + "\n")
    (wd / "input.tim").write_bytes((source / "input.tim").read_bytes())
    with (wd / "fit.log").open("w") as log:
        subprocess.run(
            [
                str((TEMPO / "bin/tempo2")),
                "-nobs",
                "20000",
                "-clkdir",
                str((INPUT / "clock")),
                "-f",
                "input.par",
                "input.tim",
                "-qrfit",
                "-norescale",
                "-output",
                "commoncov",
            ],
            cwd=wd,
            env=ENV,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    n, p, nc = map(int, (wd / "design.meta").read_text().splitlines()[0].split())
    assert nc == 0
    data = np.fromfile(wd / "toas.bin").reshape(n, 4)
    M = np.fromfile(wd / "design.bin").reshape(n, p)
    meta = np.loadtxt(wd / "design.meta", skiprows=1, dtype=int)
    fit = Fit(data, M, np.ptp(data[:, 0]))
    rows = [
        s.split("\t")
        for s in (wd / "covariance.tsv").read_text().splitlines()
        if s.startswith("parameter\t")
    ]
    k = len(rows)
    idx = [
        next(i for i, (a, c) in enumerate(meta) if a == int(r[2]) and c == 0)
        for r in rows
    ]
    pos = [list(fit.cols).index(i) for i in idx]
    fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
    C = fit.C0[np.ix_(pos, pos)] * np.outer(fac, fac)
    sd = np.sqrt(C.diagonal())
    L = np.linalg.cholesky(C / np.outer(sd, sd)) * sd[:, None]
    axis = np.vstack([np.zeros(k), 0.1 * L.T, -0.1 * L.T])
    d = delays(wd, rows, axis, "axes")
    D = (d[1 : k + 1] - d[k + 1 :]) / 0.2
    D -= D.mean(axis=1)[:, None]
    delta = d[0] - reference
    delta -= delta.mean()
    central = np.mean((phasefac * delta) ** 2)
    res = {
        "fixed": fix,
        "central_exact_loss": float(1 - abs(np.exp(1j * phasefac * delta).mean()) ** 2),
        "noise": {},
    }
    Li = np.linalg.inv(L)
    q95 = chi2.ppf(0.95, k)
    for kind in ["white", "jitter", "red", "exp"]:
        f = fit.optimize(kind)
        C = f["cov"][np.ix_(pos, pos)] * np.outer(fac, fac)
        mu = f["beta"][pos] * fac
        W = Li @ C @ Li.T
        v, V = np.linalg.eigh((W + W.T) / 2)
        R = V @ np.diag(np.sqrt(np.maximum(v, 0)))
        G = phasefac**2 * D @ D.T / 8192
        ev, evec = np.linalg.eigh(R.T @ G @ R)
        radius = np.sqrt(q95 * max(ev[-1], 0))
        shift = phasefac * (mu @ Li.T @ D)
        # Sign-independent bound includes central and GLS shifts by triangle inequality.
        bound = (np.sqrt(central) + np.sqrt(np.mean(shift**2)) + radius) ** 2
        dirs = np.sqrt(q95) * (L @ R @ evec).T
        draws = np.vstack(
            [np.zeros(k), mu, -mu, mu + dirs, mu - dirs, -mu + dirs, -mu - dirs]
        )
        exact = delays(wd, rows, draws, kind + "_ellipsoid_axes")
        loss = 1 - abs(np.exp(1j * phasefac * (exact - reference)).mean(axis=1)) ** 2
        res["noise"][kind] = {
            **serialize(f),
            "linear_nominal95_region_loss_bound": float(bound),
            "exact_tested_axis_max_loss": float(loss.max()),
            "nominal95_radius": float(radius),
        }
    results[label] = res
    (OUT / "physical_grid_results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(
        label,
        "max bound %",
        100
        * max(v["linear_nominal95_region_loss_bound"] for v in res["noise"].values()),
        flush=True,
    )
