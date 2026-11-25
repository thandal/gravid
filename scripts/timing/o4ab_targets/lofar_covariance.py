"""Common GLS noise sensitivity and native orbital phase bounds for J1400.

An independent dataset/model sensitivity check, not a rescaling of the LVK
published box. Counts here refer to the standalone LOFAR fitted region.
"""

import argparse
import json
import os
import subprocess

import numpy as np
from gls import Fit
from paths import LOFAR, TEMPO, TOOLS
from scipy.stats import chi2

BASE = LOFAR
ap = argparse.ArgumentParser()
ap.add_argument("--branch", default="lofar_fit")
args = ap.parse_args()
wd = BASE / args.branch
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
idx = [
    next(i for i, (a, b) in enumerate(meta) if a == int(r[2]) and b == 0) for r in rows
]
pos = [list(fit.cols).index(i) for i in idx]
k = len(rows)
fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
C0 = fit.C0[np.ix_(pos, pos)] * np.outer(fac, fac)
sd = np.sqrt(C0.diagonal())
L = np.linalg.cholesky(C0 / np.outer(sd, sd)) * sd[:, None]
Li = np.linalg.inv(L)
nn = 8192
draws = np.vstack([np.zeros(k), 0.1 * L.T, -0.1 * L.T])
path = wd / "phase_axes.txt"
with path.open("w") as f:
    f.write(
        f"{nn} {len(draws)} {k} 60088.625 60703.708333333333\n"
        + " ".join(r[2] for r in rows)
        + "\n"
    )
    np.savetxt(f, draws, fmt="%.19e")
with (wd / "phase_axes.log").open("w") as log:
    subprocess.run(
        [
            str((TOOLS / "delays")),
            str(wd / "native.par"),
            str(path),
            str(wd / "phase_axes.bin"),
            "unused",
        ],
        env=dict(os.environ, TEMPO2=str((TEMPO / "share/tempo2"))),
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
    )
d = np.fromfile(wd / "phase_axes.bin").reshape(1 + 2 * k, nn)
D = (d[1 : k + 1] - d[k + 1 :]) / 0.2
D -= D.mean(axis=1)[:, None]
f0 = float(
    next(
        s.split()[1]
        for s in (wd / "native.par").read_text().splitlines()
        if s.startswith("F0 ")
    )
)
J = 2 * np.pi * 2 * f0 * 1.001 * D
q = chi2.ppf(0.95, k)
result = {
    "branch": args.branch,
    "TOAs": n,
    "epochs": len(fit.t),
    "binary_names": [r[1] for r in rows],
    "noise_models": {},
}
lines = (wd / "covariance.tsv").read_text().splitlines()
native = np.loadtxt(lines[lines.index("covariance") + 1 :])
Cqr = fit.C0[np.ix_(pos, pos)] / np.outer(fit.scale[pos], fit.scale[pos])
err = np.max(
    abs(Cqr - native) / np.sqrt(np.outer(native.diagonal(), native.diagonal()))
)
assert err < 1e-5, err
result["native_independent_QR_normalized_error"] = float(err)
for kind in ["white", "jitter", "red", "exp"]:
    f = fit.optimize(kind)
    obj, b, C, K = fit.evaluate(np.array(f["theta"]), kind, True)
    B = C[np.ix_(pos, pos)] * np.outer(fac, fac)
    mu = b[pos] * fac
    W = Li @ B @ Li.T
    var = np.sum(J * (W @ J), axis=0)
    center = (mu @ Li.T) @ J
    radius = float(np.max(abs(center) + np.sqrt(q * np.maximum(var, 0))))
    bound = float(np.sin(radius) ** 2) if radius < np.pi / 2 else 1
    result["noise_models"][kind] = {
        "theta": np.asarray(f["theta"]).tolist(),
        "objective": float(f["objective"]),
        "phase_radius_rad": radius,
        "joint95_any_sampled_weight_loss_bound": bound,
        "one_template_certificate_at_3percent": bool(bound <= 0.03),
        "binary_sigma": np.sqrt(B.diagonal()).tolist(),
    }
    np.savez(wd / (kind + "_fit.npz"), theta=f["theta"], binary_cov=B, binary_shift=mu)
    print(kind, result["noise_models"][kind], flush=True)
(wd / "noise_covariance_results.json").write_text(json.dumps(result, indent=2) + "\n")
