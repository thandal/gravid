"""Use the common 200-replicate conditional calibration on the LOFAR fit."""

import json

import numpy as np
from gls import Fit
from paths import LOFAR
from scipy.stats import binom, chi2, f

wd = LOFAR / "lofar_fit"
n, p, nc = map(int, (wd / "design.meta").read_text().splitlines()[0].split())
data = np.fromfile(wd / "toas.bin").reshape(n, 4)
M = np.fromfile(wd / "design.bin").reshape(n, p)
meta = np.loadtxt(wd / "design.meta", skiprows=1, dtype=int)
fit = Fit(data, M, np.ptp(data[:, 0]))
base = fit.optimize("red")
rows = [
    s.split("\t")
    for s in (wd / "covariance.tsv").read_text().splitlines()
    if s.startswith("parameter\t")
]
idx = [
    next(i for i, (a, b) in enumerate(meta) if a == int(r[2]) and b == 0) for r in rows
]
pos = [list(fit.cols).index(i) for i in idx]
fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
rng = np.random.default_rng(2026101801)
root = np.linalg.cholesky(base["K"])
qs = []
thetas = []
for b in range(200):
    y = base["EFAC"] * rng.normal(size=n) + fit.W @ root @ rng.normal(size=len(fit.t))
    fit.set_y(y)
    boot = fit.optimize("red", start=base["theta"], quick=True)
    mu = boot["beta"][pos] * fac
    C = boot["cov"][np.ix_(pos, pos)] * np.outer(fac, fac)
    sd = np.sqrt(C.diagonal())
    z = mu / sd
    qs.append(float(z @ np.linalg.solve(C / np.outer(sd, sd), z)))
    thetas.append(boot["theta"])
    if (b + 1) % 25 == 0:
        print(b + 1, flush=True)
q = np.array(qs)
critical = max(
    float(chi2.ppf(0.95, len(pos))),
    float(np.sort(q)[int(binom.ppf(0.95, len(q), 0.95))]),
)
result = {
    "replicates": len(q),
    "seed": 2026101801,
    "noise_refitted": True,
    "joint95_nominal_coverage": float(np.mean(q < chi2.ppf(0.95, len(pos)))),
    "joint_q95": float(np.quantile(q, 0.95)),
    "critical_q": critical,
    "upper_quantile_order": int(binom.ppf(0.95, len(q), 0.95)) + 1,
    "white_only_exact_f_radius": float(len(pos) * f.ppf(0.95, len(pos), fit.df)),
    "scope": "Conditional on the fitted generating noise and the retained author-filtered dataset; does not simulate the original RFI/outlier-selection process.",
}
np.savez(wd / "bootstrap.npz", q=q, theta=thetas)
(wd / "bootstrap_results.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
