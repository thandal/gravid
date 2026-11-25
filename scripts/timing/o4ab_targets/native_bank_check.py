"""Check the limiting orbital regions on dense grids with TEMPO2 binary delays."""

import argparse
import json
import os
import subprocess

import numpy as np
from paths import BASE, HERE, NAMES, OUT, TEMPO, TOOLS

SAVED = HERE / "results"


def native_delays(model, parameter_ids, offsets, count, directory, label):
    spec = directory / (label + ".txt")
    binary = directory / (label + ".bin")
    with spec.open("w") as output:
        output.write(
            f"{count} {len(offsets)} {len(parameter_ids)} 60088.625 {60703 + 17 / 24:.16f}\n"
        )
        output.write(" ".join(map(str, parameter_ids)) + "\n")
        np.savetxt(output, offsets, fmt="%.19e")
    with (directory / (label + ".log")).open("w") as log:
        subprocess.run(
            [TOOLS / "delays", model, spec, binary, "unused"],
            env=dict(os.environ, TEMPO2=str(TEMPO / "share/tempo2")),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return np.fromfile(binary).reshape(len(offsets), count)


def limiting_fit(name):
    data = np.load(OUT / (name + "_bank_worst.npz"))
    reference = BASE / (name + "_ELL1")
    label = str(data["physical_dir"])
    fit = reference if label == reference.name else OUT / "physical_grid" / label
    parameters = [
        line.split()
        for line in (fit / "covariance.tsv").read_text().splitlines()
        if line.startswith("parameter\t")
    ]
    return {
        "reference": reference / "native.par",
        "model": fit / "native.par",
        "parameters": [r[1] for r in parameters],
        "parameter_ids": [int(r[2]) for r in parameters],
        "covariance": data["binary_cov"].tolist(),
        "shift": data["binary_shift"].tolist(),
        "critical_q": float(data["critical_q"]),
    }


def check(fit, count, directory, rng, nonlinear):
    sd = np.sqrt(np.diag(fit["covariance"]))
    root = np.linalg.cholesky(fit["covariance"] / np.outer(sd, sd)) * sd[:, None]
    k = len(sd)
    scale = (
        2
        * np.pi
        * 2
        * float(
            next(
                line.split()[1]
                for line in fit["model"].read_text().splitlines()
                if line.startswith("F0 ")
            )
        )
        * 1.001
    )

    def delays(model, offsets, label):
        return native_delays(
            model, fit["parameter_ids"], offsets, count, directory, label
        )

    ref = delays(fit["reference"], np.zeros((1, k)), "reference")[0]
    axes = delays(
        fit["model"], np.vstack([np.zeros(k), 0.1 * root.T, -0.1 * root.T]), "axes"
    )
    jacobian = scale * (axes[1 : k + 1] - axes[k + 1 :]) / 0.2
    jacobian -= jacobian.mean(axis=1, keepdims=True)
    central = scale * (axes[0] - ref)
    central -= central.mean()
    shift = np.linalg.solve(root, fit["shift"])
    center = central + shift @ jacobian
    q = fit["critical_q"]
    envelope = abs(center) + np.sqrt(q * np.sum(jacobian**2, axis=0))
    radius = float(envelope.max())
    result = {
        "phase_radius_rad": radius,
        "any_sampled_weight_loss_bound": float(np.sin(min(radius, np.pi / 2)) ** 2),
    }
    if not nonlinear:
        return result, None
    extremal = jacobian[:, np.argmax(envelope)]
    extremal = extremal / np.linalg.norm(extremal)
    directions = rng.normal(size=(256, k))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    z = shift + np.sqrt(q) * np.vstack([directions, extremal, -extremal])
    phase = scale * (delays(fit["model"], z @ root.T, "boundary") - ref)
    phase -= phase.mean(axis=1, keepdims=True)
    return result, {
        "boundary_draws": len(z),
        "max_phase_linearization_error_rad": float(
            np.max(abs(phase - (central + z @ jacobian)))
        ),
        "largest_exact_uniform_loss": float(
            np.max(1 - abs(np.exp(1j * phase).mean(axis=1)) ** 2)
        ),
        "max_exact_phase_radius_rad": float(np.max(abs(phase))),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--high-only", action="store_true")
    parser.add_argument(
        "--from-saved",
        action="store_true",
        help="Read the published fitted regions and timing models",
    )
    args = parser.parse_args()
    out = OUT / ("saved_replay" if args.from_saved else "native_bank")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(72026)
    results = {}
    for name in NAMES:
        if args.from_saved:
            fit = json.loads((SAVED / "limiting_orbits.json").read_text())[name]
            for key in ("model", "reference"):
                fit[key] = SAVED / fit[key]
        else:
            fit = limiting_fit(name)
        per_target = {}
        for n in [131072] if args.high_only else [8192, 32768]:
            directory = out / name / str(n)
            directory.mkdir(parents=True, exist_ok=True)
            per_target[str(n)], boundary = check(fit, n, directory, rng, n == 32768)
            if boundary:
                per_target["nonlinear"] = boundary
        results[name] = per_target
        print(name, per_target, flush=True)
    filename = (
        "native_bank_high_resolution.json"
        if args.high_only
        else "native_bank_checks.json"
    )
    destination = out if args.from_saved else OUT
    (destination / filename).write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
