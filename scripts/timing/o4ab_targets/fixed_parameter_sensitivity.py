"""Refit with PX and F2 free: a shared deterministic-model sensitivity test."""

import json
import os
import subprocess

import numpy as np
from gls import Fit
from paths import BASE, INPUT, NAMES, OUT, TEMPO, TOOLS

out = OUT / "free_px_f2"
out.mkdir(exist_ok=True)
ENV = dict(
    os.environ,
    TEMPO2=str((TEMPO / "share/tempo2")),
    TEMPO2_PLUG_PATH=str(TOOLS),
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
)
profiles = json.loads((OUT / "noise_profile_results.json").read_text())
banks = json.loads((OUT / "bank_envelope_results.json").read_text())
result = {}
for name in NAMES:
    base = BASE / (name + "_ELL1")
    wd = out / name
    wd.mkdir(exist_ok=True)
    ls = []
    values = {
        s.split()[0]: s.split()[1]
        for s in (base / "native.par").read_text().splitlines()
        if len(s.split()) > 1
    }
    for s in (base / "native.par").read_text().splitlines():
        a = s.split()
        if not a or a[0] in [
            "START",
            "FINISH",
            "NTOA",
            "TRES",
            "CHI2R",
            "NITS",
            "PX",
            "F2",
        ]:
            continue
        ls.append(s)
    ls += [
        "PX " + values.get("PX", "0") + " 1",
        "F2 " + values.get("F2", "0") + " 1",
        "NITS 5",
    ]
    (wd / "input.par").write_text("\n".join(ls) + "\n")
    (wd / "input.tim").write_bytes((base / "input.tim").read_bytes())
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
    C = fit.C0[np.ix_(pos, pos)] * np.outer(fac, fac)
    sd = np.sqrt(C.diagonal())
    L = np.linalg.cholesky(C / np.outer(sd, sd)) * sd[:, None]
    Li = np.linalg.inv(L)
    axis = np.vstack([np.zeros(9), 0.1 * L.T, -0.1 * L.T])
    path = wd / "axes.txt"
    with path.open("w") as f:
        f.write(
            f"8192 19 9 60088.625 {60703 + 17 / 24:.16f}\n"
            + " ".join(r[2] for r in rows)
            + "\n"
        )
        np.savetxt(f, axis, fmt="%.19e")
    with (wd / "axes.log").open("w") as log:
        subprocess.run(
            [
                str((TOOLS / "delays")),
                str(wd / "native.par"),
                str(path),
                str(wd / "axes.bin"),
                "unused",
            ],
            env=ENV,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    d = np.fromfile(wd / "axes.bin").reshape(19, 8192)
    ref = np.fromfile(base / "axes.bin").reshape(19, 8192)[0]
    ph = 2 * np.pi * 2 * float(values["F0"]) * 1.001
    D = ph * (d[1:10] - d[10:]) / 0.2
    D -= D.mean(axis=1)[:, None]
    central = ph * (d[0] - ref)
    central -= central.mean()
    critical = banks[name]["critical_q"]
    worst = 0
    for row in profiles[name]["cases"]:
        fit.span = np.ptp(data[:, 0]) * row["span_factor"]
        th = np.array(
            [
                np.log(max(1.1, row["EFAC"])),
                np.log(row["jitter_us"]),
                np.log(row["red_rms_us"]),
                row["gamma"],
            ]
        )
        obj, b, C, K = fit.evaluate(th, "red", True)
        Cb = C[np.ix_(pos, pos)] * np.outer(fac, fac)
        mu = b[pos] * fac
        W = Li @ Cb @ Li.T
        center = central + (mu @ Li.T) @ D
        v = np.sum(D * (W @ D), axis=0)
        radius = np.max(abs(center) + np.sqrt(critical * np.maximum(v, 0)))
        bound = float(np.sin(radius) ** 2) if radius < np.pi / 2 else 1.0
        worst = max(worst, bound)
    result[name] = {
        "free_extra_parameters": ["PX", "F2"],
        "central_binary_uniform_loss": float(1 - abs(np.exp(1j * central).mean()) ** 2),
        "largest_profile_any_sampled_weight_bound": worst,
        "scope": "Model sensitivity at the reference physical case; same bootstrap critical radius and reference noise-profile grid, not a new coverage calibration or a cross with every J1231 mass/inclination case.",
    }
    print(name, result[name], flush=True)
    (OUT / "fixed_parameter_sensitivity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
