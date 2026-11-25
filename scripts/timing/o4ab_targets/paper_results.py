"""Write compact publication data, verify a full rerun, or recompute Table 1."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from native_bank_check import SAVED, limiting_fit
from paths import HERE, INPUT, LOFAR, NAMES, OUT, WORK


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def table(bounds, predictions):
    result = {}
    for name in NAMES:
        cases = {
            k: v for k, v in predictions[name].items() if not k.endswith("_control")
        }
        losses = []
        for case in cases.values():
            phasors = np.array(case["epoch_phasors"])
            if phasors.shape != (case["epochs"], 2):
                raise ValueError("Epoch count disagrees with the phase data")
            mean = phasors.mean(axis=0)
            losses.append(float(1 - mean @ mean))
        ref = cases["earlier_ELL1"]
        result[name] = {
            "orbital_loss_bound": bounds[name]["131072"][
                "any_sampled_weight_loss_bound"
            ],
            "epochs": ref["epochs"],
            "TOAs": ref["TOAs"],
            "radio_models": len(cases),
            "min_radio_loss": min(losses),
            "max_radio_loss": max(losses),
        }
    return result


def epoch_predictions():
    predictions = read(WORK / "radio/predictions.json")
    for name, cases in predictions.items():
        epochs = {}
        for line in (INPUT / (name + "_o4b.tim")).read_text().splitlines():
            words = line.split()
            if len(words) >= 5:
                epochs.setdefault(words[0], []).append(float(words[2]))
        ref = cases["earlier_ELL1"]
        if (
            len(epochs) != ref["epochs"]
            or sum(map(len, epochs.values())) != ref["TOAs"]
        ):
            raise ValueError(f"Inconsistent observation list for {name}")
        ref["epoch_mean_mjd_utc"] = [float(np.mean(times)) for times in epochs.values()]
        ref["channels_per_epoch"] = [len(times) for times in epochs.values()]
    return predictions


def full_results():
    bank = LOFAR / "bank_efficiency"
    fits = read(OUT / "noise_results.json")
    bands = read(OUT / "frequency_checks.json")
    return {
        "timing_splits": {
            name: {
                "frequency_halves_p": bands[name]["nominal_p"],
                "quarter_orbital_p": [
                    fits[name]["blocked_prediction"][f"q{i}"]["red"][
                        "nominal_orbital_p"
                    ]
                    for i in range(1, 5)
                ],
            }
            for name in NAMES
        },
        "table_1": table(
            read(OUT / "native_bank_high_resolution.json"),
            read(WORK / "radio/predictions.json"),
        ),
        "bootstrap": read(OUT / "bootstrap_results.json"),
        "noise_grid": read(OUT / "bank_envelope_results.json"),
        "grid_and_nonlinear_checks": read(OUT / "native_bank_checks.json"),
        "extra_parameters": read(OUT / "fixed_parameter_sensitivity.json"),
        "radio_prediction_statistics": read(
            WORK / "radio/external_predictive_results.json"
        ),
        "J0437": read(WORK / "j0437_comparison.json"),
        "historical": read(WORK / "historical/results.json"),
        "LOFAR_placement": read(bank / "uniform_slabs_bank_results.json"),
        "LOFAR_native_checks": read(bank / "uniform_slabs_native_validation.json"),
    }


def package():
    SAVED.mkdir(parents=True, exist_ok=True)
    orbits = SAVED / "orbits"
    orbits.mkdir(exist_ok=True)
    shutil.copyfile(INPUT / "J1400-1431_lwa.tim", SAVED / "J1400-1431_lwa.tim")
    fits = {}
    for name in NAMES:
        fit = limiting_fit(name)
        reference = fit["reference"]
        for key in ("reference", "model"):
            filename = name + (
                "_reference.par" if fit[key] == reference else "_limiting.par"
            )
            shutil.copyfile(fit[key], orbits / filename)
            fit[key] = str(Path("orbits") / filename)
        fits[name] = fit
    write(SAVED / "limiting_orbits.json", fits)
    write(SAVED / "epoch_predictions.json", epoch_predictions())
    bank = LOFAR / "bank_efficiency"
    shutil.copyfile(
        bank / "uniform_slabs_bank_offsets.csv", SAVED / "lofar_offsets.csv"
    )
    shutil.copyfile(LOFAR / "lofar_fit/native.par", orbits / "J1400-1431_lofar.par")
    write(SAVED / "paper_results.json", full_results())


def compare(actual, expected, key=""):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f"Result fields differ at {key}")
        for k, value in expected.items():
            compare(actual[k], value, f"{key}/{k}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"Result lengths differ at {key}")
        for i, value in enumerate(expected):
            compare(actual[i], value, f"{key}/{i}")
    elif isinstance(expected, float):
        np.testing.assert_allclose(actual, expected, rtol=5e-5, atol=1e-9, err_msg=key)
    elif actual != expected:
        raise ValueError(f"Result differs at {key}: {actual} vs {expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--package",
        action="store_true",
        help="Write publication files from a completed timing run",
    )
    mode.add_argument(
        "--verify-run",
        action="store_true",
        help="Compare a completed full run with the publication files",
    )
    args = parser.parse_args()
    if args.package:
        package()
        return
    expected = read(SAVED / "paper_results.json")
    if args.verify_run:
        result = full_results()
        compare(result, expected)
        compare(epoch_predictions(), read(SAVED / "epoch_predictions.json"))
        print("Full timing run agrees with the saved results.")
        rows = result["table_1"]
    else:
        subprocess.run(
            [
                sys.executable,
                HERE / "native_bank_check.py",
                "--from-saved",
                "--high-only",
            ],
            check=True,
        )
        rows = table(
            read(OUT / "saved_replay/native_bank_high_resolution.json"),
            read(SAVED / "epoch_predictions.json"),
        )
        compare(rows, expected["table_1"])
    for name, row in rows.items():
        print(
            f"{name}   orbital bound {100 * row['orbital_loss_bound']:.4f}%   "
            f"{row['epochs']} epochs ({row['TOAs']} TOAs)   "
            f"radio loss {100 * row['min_radio_loss']:.5f}% to {100 * row['max_radio_loss']:.5f}%"
        )


if __name__ == "__main__":
    main()
