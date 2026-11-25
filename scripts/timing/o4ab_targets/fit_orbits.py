"""Refit the three MeerKAT orbits and the free/fixed-orbit LOFAR models."""

import argparse
import os
import shutil
import subprocess

from paths import BASE, HERE, INPUT, LOFAR, NAMES, TEMPO, TOOLS


def fit(seed, toas, clock, directory):
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(seed, directory / "input.par")
    shutil.copyfile(toas, directory / "input.tim")
    env = dict(
        os.environ, TEMPO2=str(TEMPO / "share/tempo2"), TEMPO2_PLUG_PATH=str(TOOLS)
    )
    cmd = [
        TEMPO / "bin/tempo2",
        "-nobs",
        "20000",
        "-clkdir",
        clock,
        "-f",
        "input.par",
        "input.tim",
        "-qrfit",
        "-norescale",
        "-output",
        "commoncov",
    ]
    with (directory / "fit.log").open("w") as log:
        subprocess.run(
            cmd,
            env=env,
            cwd=directory,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    print(directory.name, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", choices=["meerkat", "lofar"])
    args = parser.parse_args()
    if args.sample == "meerkat":
        for name in NAMES:
            fit(
                HERE / "models" / (name + "_seed.par"),
                INPUT / (name + "_training.tim"),
                INPUT / "clock",
                BASE / (name + "_ELL1"),
            )
    else:
        for branch in ["lofar_fit", "lofar_fit_fixed_lvk", "lofar_fit_no_secular"]:
            fit(
                HERE / "models" / (branch + "_seed.par"),
                INPUT / "J1400-1431_training.tim",
                INPUT / "low_frequency_clock",
                LOFAR / branch,
            )
