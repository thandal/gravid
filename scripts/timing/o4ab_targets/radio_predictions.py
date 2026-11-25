"""Predict the O4b arrival times with the earlier fits, without a timing refit."""

import json
import os
import shutil
import subprocess

import numpy as np
from paths import BASE, INPUT, NAMES, OUT, TEMPO, TOOLS, WORK


def prediction(directory, label, model, observations, group_ids, change=None):
    lines = []
    for line in model.read_text().splitlines():
        words = line.split()
        if (
            not words
            or words[0] in {"START", "FINISH", "NITS", "NTOA", "TRES", "CHI2R"}
            or words[0].startswith("DMX")
        ):
            continue
        if change and words[0] in change:
            line = f"{words[0]} {float(words[1]) + change[words[0]]:.18g}"
        lines.append(line)
    spin = float(next(line.split()[1] for line in lines if line.startswith("F0 ")))
    par = directory / (label + ".par")
    par.write_text("\n".join(lines) + "\n")
    residuals = directory / (label + ".res")
    with (directory / (label + ".log")).open("w") as log:
        subprocess.run(
            [TOOLS / "evaluate", par, observations, INPUT / "clock", residuals],
            env=dict(os.environ, TEMPO2=str(TEMPO / "share/tempo2")),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    data = np.loadtxt(residuals)
    if len(data) != len(group_ids):
        raise ValueError("Different numbers of TOAs and residuals")
    epoch_phasors = []
    for epoch in dict.fromkeys(group_ids):
        group = data[group_ids == epoch]
        epoch_phasors.append(
            np.average(
                np.exp(2j * np.pi * 2 * spin * group[:, 2]), weights=group[:, 3] ** -2
            )
        )
    return {
        "TOAs": len(data),
        "epochs": len(epoch_phasors),
        "raw_sampled_loss": float(1 - abs(np.mean(epoch_phasors)) ** 2),
        "epoch_phasors": [[float(z.real), float(z.imag)] for z in epoch_phasors],
    }


def main():
    results = {}
    for name in NAMES:
        stream = "RelBin" if name == "J1231-1411" else "PTA"
        directory = WORK / "radio" / "predictions" / name / stream / "O4b"
        directory.mkdir(parents=True, exist_ok=True)
        observations = directory / "input.tim"
        shutil.copyfile(INPUT / (name + "_o4b.tim"), observations)
        ids = np.array(
            [
                s.split()[0]
                for s in observations.read_text().splitlines()
                if len(s.split()) > 4
            ]
        )
        ref = BASE / (name + "_ELL1") / "native.par"
        cases = {
            "earlier_ELL1": prediction(
                directory, "earlier_ELL1", ref, observations, ids
            )
        }
        if name == "J1231-1411":
            grid = sorted((OUT / "physical_grid").glob("*/native.par"))
            if len(grid) != 24:
                raise ValueError("Expected 24 fitted mass and inclination combinations")
            for par in grid:
                cases[par.parent.name] = prediction(
                    directory, par.parent.name, par, observations, ids
                )
        else:
            params = {
                s.split()[0]: s.split()[1]
                for s in ref.read_text().splitlines()
                if len(s.split()) >= 2
            }
            for key, delta in [
                ("A1", 0.0003),
                (
                    "TASC",
                    0.0003 * float(params["PB"]) / (2 * np.pi * float(params["A1"])),
                ),
            ]:
                label = key + "_300us_control"
                cases[label] = prediction(
                    directory, label, ref, observations, ids, {key: delta}
                )
        results[name] = cases
    (WORK / "radio/predictions.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
