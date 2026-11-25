"""Evaluate orbital-delay derivatives about the refitted MeerKAT reference orbits."""

import os
import subprocess

import numpy as np
from gls import Fit, load
from paths import NAMES, TEMPO, TOOLS


def orbital_axes(directory, data, design, parameters):
    fit = Fit(data, design, np.ptp(data[:, 0]))
    rows = [
        line.split()
        for line in (directory / "covariance.tsv").read_text().splitlines()
        if line.startswith("parameter\t")
    ]
    indices = [
        next(
            i
            for i, (key, counter) in enumerate(parameters)
            if key == int(row[2]) and counter == 0
        )
        for row in rows
    ]
    positions = [list(fit.cols).index(i) for i in indices]
    conversion = (
        np.array([float(row[4]) / float(row[5]) for row in rows]) / fit.scale[positions]
    )
    covariance = fit.C0[np.ix_(positions, positions)] * np.outer(conversion, conversion)
    sd = np.sqrt(covariance.diagonal())
    root = np.linalg.cholesky(covariance / np.outer(sd, sd)) * sd[:, None]
    offsets = np.vstack([np.zeros(len(rows)), 0.1 * root.T, -0.1 * root.T])
    np.savez(directory / "robust_covariances.npz", formal_QR=covariance)
    with (directory / "axes.txt").open("w") as output:
        output.write(f"8192 {len(offsets)} {len(rows)} 60088.625 60703.708333333336\n")
        output.write(" ".join(row[2] for row in rows) + "\n")
        np.savetxt(output, offsets, fmt="%.19e")
    subprocess.run(
        [
            TOOLS / "delays",
            directory / "native.par",
            directory / "axes.txt",
            directory / "axes.bin",
            "unused",
        ],
        env=dict(os.environ, TEMPO2=str(TEMPO / "share/tempo2")),
        check=True,
    )


if __name__ == "__main__":
    for name in NAMES:
        orbital_axes(*load(name))
