"""Native nonlinear checks of new banks; point samples do not certify coverage."""

import json
import os
import subprocess

import numpy as np
from paths import LOFAR, TEMPO, TOOLS

WD = LOFAR / "lofar_fit"
OUT = WD.parent / "bank_efficiency"
label = "uniform_slabs"
bank = np.load(OUT / (label + "_bank.npz"))
A = bank["physical_axes"]
mu = bank["mean_shift"]
R = np.sqrt(float(bank["critical_q"]))
centers = bank["coordinates"]
rows = [
    r.split("\t")
    for r in (WD / "covariance.tsv").read_text().splitlines()
    if r.startswith("parameter\t")
]
f0 = float(
    next(
        s.split()[1]
        for s in (WD / "native.par").read_text().splitlines()
        if s.startswith("F0 ")
    )
)
ph = 2 * np.pi * 2 * f0 * 1.001
rng = np.random.default_rng(2026101804)
v = rng.normal(size=(4096, 9))
v *= R / np.linalg.norm(v, axis=1)[:, None]
v[2048:] *= rng.random(2048)[:, None] ** (1 / 9)
v = np.vstack([v, np.eye(9) * R, -np.eye(9) * R, np.zeros((1, 9))])
extra = []
for endpoint in bank["edges"]:
    rem = max(R * R - endpoint**2, 0)
    for axis in [1, 2]:
        for sign in [-1, 1]:
            p = np.zeros(9)
            p[0] = endpoint
            p[axis] = sign * np.sqrt(rem)
            extra.append(p)
v = np.vstack([v, extra])
slab = np.clip(
    np.searchsorted(bank["edges"], v[:, 0], side="right") - 1, 0, len(bank["edges"]) - 2
)
assigned = bank["slab_assignment"][slab, (v[:, 1] >= 0).astype(int)]
truth = mu + v @ A.T
templates = bank["binary_offsets"][assigned]
# Project the new bank's physical offsets into the preserved dense native map.
old = np.load(WD / "bank.npz")
nn = 8192


# Re-evaluate Jacobian exactly on the validation times, rather than interpolating
# it across a different grid when judging tiny nonlinear errors.
def native(draws, nt, name, retain=False):
    path = OUT / (name + ".txt")
    binary = OUT / (name + ".bin")
    with path.open("w") as f:
        f.write(
            f"{nt} {len(draws)} 9 60088.625 60703.708333333333\n"
            + " ".join(r[2] for r in rows)
            + "\n"
        )
        np.savetxt(f, draws, fmt="%.19e")
    with (OUT / (name + ".log")).open("w") as log:
        subprocess.run(
            [
                str((TOOLS / "delays")),
                str(WD / "native.par"),
                str(path),
                str(binary),
                "unused",
            ],
            env=dict(os.environ, TEMPO2=str((TEMPO / "share/tempo2"))),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    d = np.fromfile(binary).reshape(len(draws), nt)
    if not retain:
        binary.unlink()
        path.unlink()
    return d


axes = native(
    np.vstack(
        [np.zeros(9), 0.1 * old["physical_axes"].T, -0.1 * old["physical_axes"].T]
    ),
    nn,
    label + "_validation_axes",
)
J = ph * (axes[1:10] - axes[10:]) / 0.2
J -= J.mean(axis=1, keepdims=True)
H = np.linalg.solve(old["physical_axes"].T, J)
losses = []
bounds = []
errors = []
for start in range(0, len(v), 256):
    stop = min(start + 256, len(v))
    t = truth[start:stop]
    b = templates[start:stop]
    d = native(np.vstack([t, b]), nn, label + "_validation_batch")
    phase = ph * (d[: len(t)] - d[len(t) :])
    phase -= phase.mean(axis=1, keepdims=True)
    linear = (t - b) @ H
    error = np.max(abs(phase - linear), axis=1)
    loss = 1 - abs(np.exp(1j * phase).mean(axis=1)) ** 2
    radius = np.max(abs(phase), axis=1)
    losses.extend(loss)
    bounds.extend(np.sin(np.minimum(radius, np.pi / 2)) ** 2)
    errors.extend(error)
    if stop % 1024 == 0:
        print(label, stop, "of", len(v), flush=True)
np.savez(
    OUT / (label + "_native_checks.npz"),
    truth_offsets=truth,
    template_indices=assigned,
    uniform_losses=losses,
    any_weight_bounds=bounds,
    linearization_errors_rad=errors,
)
result = {
    "templates": len(centers),
    "tested_truths": len(v),
    "times_per_truth": nn,
    "max_exact_native_uniform_loss": float(np.max(losses)),
    "max_exact_native_any_sampled_weight_loss_bound": float(np.max(bounds)),
    "max_centered_native_linearization_error_rad": float(np.max(errors)),
    "seed": 2026101804,
    "scope": "Native boundary, interior, axis and selected cell-extremum checks; complete region coverage comes from the separate linear certificate, with sampled support for the nonlinear reserve.",
}
assert result["max_centered_native_linearization_error_rad"] < 1e-4
assert result["max_exact_native_uniform_loss"] < 0.03
(OUT / (label + "_native_validation.json")).write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
