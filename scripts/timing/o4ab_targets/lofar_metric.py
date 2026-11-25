"""Calculate the uniform-time phase metric of the calibrated LOFAR region."""

import json
import os
import subprocess

import numpy as np
from gls import Fit
from paths import LOFAR, TEMPO, TOOLS

wd = LOFAR / "lofar_fit"
out = LOFAR / "bank_efficiency"
out.mkdir(exist_ok=True)
n, p, _ = map(int, (wd / "design.meta").read_text().splitlines()[0].split())
meta = np.loadtxt(wd / "design.meta", skiprows=1, dtype=int)
data = np.fromfile(wd / "toas.bin").reshape(n, 4)
fit = Fit(data, np.fromfile(wd / "design.bin").reshape(n, p), np.ptp(data[:, 0]))
rows = [
    s.split()
    for s in (wd / "covariance.tsv").read_text().splitlines()
    if s.startswith("parameter\t")
]
indices = [
    next(i for i, (key, count) in enumerate(meta) if key == int(row[2]) and count == 0)
    for row in rows
]
positions = [list(fit.cols).index(i) for i in indices]
factors = (
    np.array([float(row[4]) / float(row[5]) for row in rows]) / fit.scale[positions]
)
cov = fit.C0[np.ix_(positions, positions)] * np.outer(factors, factors)
sd = np.sqrt(cov.diagonal())
L = np.linalg.cholesky(cov / np.outer(sd, sd)) * sd[:, None]
red = np.load(wd / "red_fit.npz")
whitened = np.linalg.solve(L, np.linalg.solve(L, red["binary_cov"]).T).T
root = np.linalg.cholesky((whitened + whitened.T) / 2)
spin = float(
    next(
        s.split()[1]
        for s in (wd / "native.par").read_text().splitlines()
        if s.startswith("F0 ")
    )
)
phase = 2 * np.pi * 2 * spin * 1.001


def phase_derivatives(delays):
    derivatives = phase * (delays[1:10] - delays[10:]) / 0.2
    return derivatives - derivatives.mean(axis=1, keepdims=True)


coarse = phase_derivatives(np.fromfile(wd / "phase_axes.bin").reshape(19, 8192))
coarse = root.T @ coarse
_, directions = np.linalg.eigh(coarse @ coarse.T)
axes = L @ root @ directions[:, ::-1]
np.savez(wd / "bank.npz", physical_axes=axes)
offsets = np.vstack([np.zeros(9), 0.1 * axes.T, -0.1 * axes.T])
with (wd / "dense.txt").open("w") as stream:
    stream.write(
        "131072 19 9 60088.625 60703.708333333333\n"
        + " ".join(row[2] for row in rows)
        + "\n"
    )
    np.savetxt(stream, offsets, fmt="%.19e")
subprocess.run(
    [TOOLS / "delays", wd / "native.par", wd / "dense.txt", wd / "dense.bin", "unused"],
    check=True,
    env=dict(os.environ, TEMPO2=str(TEMPO / "share/tempo2")),
)
dense = phase_derivatives(np.fromfile(wd / "dense.bin").reshape(19, 131072))
values, vecs = np.linalg.eigh(dense @ dense.T / dense.shape[1])
q = json.loads((wd / "bootstrap_results.json").read_text())["critical_q"]
np.savez(
    out / "uniform_bank.npz",
    eigenvalues=np.maximum(values[::-1], 0),
    physical_axes=axes @ vecs[:, ::-1],
    critical_q=q,
    mean_shift=red["binary_shift"],
)
