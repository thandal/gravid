"""Check the published five-parameter orbits and the independent J1400/J1745 data."""

import io
import json
import os
import subprocess

import numpy as np
from gls import Fit
from paths import HERE, INPUT, LOFAR, TEMPO, TOOLS, WORK
from scipy.stats import f

MODELS = HERE / "models"
OUTPUT = WORK / "historical"
ORBIT_KEYS = ("PB", "A1", "TASC", "EPS1", "EPS2")


def fields(model):
    return {
        words[0]: words[1:]
        for line in model.splitlines()
        if len(words := line.split()) >= 2
    }


def ell1_box(path, sigma):
    """Bound the varying part of the first-order Roemer delay on the whole time interval."""
    parameters = fields(path.read_text())
    p = {key: float(parameters[key][0]) for key in ORBIT_KEYS}
    errors = {key: sigma * float(parameters[key][2]) for key in ORBIT_KEYS}
    min_period = p["PB"] - errors["PB"]
    span = max(abs(t - p["TASC"]) for t in (60088.625, 60703 + 17 / 24))
    angle = (
        2
        * np.pi
        * (span * errors["PB"] / (p["PB"] * min_period) + errors["TASC"] / min_period)
    )
    ecc = abs(p["EPS1"]) + abs(p["EPS2"])
    delta_ecc = errors["EPS1"] + errors["EPS2"]
    delay = errors["A1"] * (1 + (ecc + delta_ecc) / 2) + p["A1"] * (
        (1 + ecc) * angle + delta_ecc / 2
    )
    radius = 2 * np.pi * 2 * float(parameters["F0"][0]) * 1.001 * delay
    return {
        "delay_seconds": delay,
        "loss_bound": float(np.sin(min(radius, np.pi / 2)) ** 2),
    }


def residuals(label, model, tim):
    par = OUTPUT / (label + ".par")
    par.write_text(model)
    output = OUTPUT / (label + ".res")
    with (OUTPUT / (label + ".log")).open("w") as log:
        subprocess.run(
            [
                TOOLS / "evaluate",
                par,
                INPUT / tim,
                INPUT / "low_frequency_clock",
                output,
            ],
            env=dict(os.environ, TEMPO2=str(TEMPO / "share/tempo2")),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return np.loadtxt(output)


def historical_predictions():
    from pint import logging
    from pint.models import get_model

    text = (MODELS / "J1745-0952_folding.par").read_text()
    tcb = text.replace("UNITS TDB", "UNITS TCB")
    logging.setup(level="ERROR")
    converted = get_model(io.StringIO(tcb), allow_tcb=True, allow_T2=True).as_parfile()
    models = {"labelled_TDB": text, "numeric_TCB": tcb, "converted_TDB": converted}
    result = {}
    times = {}
    for label, model in models.items():
        data = residuals(label, model, "J1745-0952_historical.tim")
        times[label] = data[:, 2]
        spin = float(fields(model)["F0"][0])
        days = np.floor(data[:, 0])
        z = np.exp(2j * np.pi * 2 * spin * data[:, 2])
        epochs = [
            np.average(z[days == day], weights=data[days == day, 3] ** -2)
            for day in np.unique(days)
        ]
        result[label] = {
            "epochs": len(epochs),
            "sampled_loss": float(1 - abs(np.mean(epochs)) ** 2),
        }
    diff = times["converted_TDB"] - times["numeric_TCB"]
    result["conversion_difference_ns"] = float(1e9 * np.max(abs(diff - diff.mean())))
    return result


def lofar_diagnostics():
    def load_fit(directory):
        n, p, _ = map(
            int, (directory / "design.meta").read_text().splitlines()[0].split()
        )
        data = np.fromfile(directory / "toas.bin").reshape(n, 4)
        return Fit(
            data,
            np.fromfile(directory / "design.bin").reshape(n, p),
            np.ptp(data[:, 0]),
        )

    wd = LOFAR / "lofar_fit"
    free = load_fit(wd)
    result = {}
    for label in ("lofar_fit_fixed_lvk", "lofar_fit_no_secular"):
        reduced = load_fit(LOFAR / label)
        count = reduced.df - free.df
        stat = (reduced.q - free.q) / count / (free.q / free.df)
        result[label] = {
            "F": float(stat),
            "dof": [int(count), int(free.df)],
            "nominal_p": float(f.sf(stat, count, free.df)),
        }
    rows = [
        line.split()
        for line in (wd / "covariance.tsv").read_text().splitlines()
        if line.startswith("parameter\t")
    ]
    i = next(i for i, row in enumerate(rows) if row[1] == "PBDOT")
    noise = np.load(wd / "red_fit.npz")
    scale = 365.25 * 86400 * 1000
    result["period_drift_ms_per_year"] = {
        "mean": float((float(rows[i][3]) + noise["binary_shift"][i]) * scale),
        "sigma": float(np.sqrt(noise["binary_cov"][i, i]) * scale),
    }
    return result


def main():
    OUTPUT.mkdir(exist_ok=True, parents=True)
    boxes = {
        name: {str(sigma): ell1_box(MODELS / par, sigma) for sigma in ns}
        for name, par, ns in [
            ("J1400-1431", "J1400-1431.par", [3, 8]),
            ("J1745-0952", "J1745-0952_historical_only.par", [3]),
        ]
    }
    model = (MODELS / "J1400-1431_lwa.par").read_text()
    data = residuals("lwa", model, "J1400-1431_lwa.tim")
    frequency = 2 * float(fields(model)["F0"][0])
    phases = np.exp(2j * np.pi * frequency * data[:, 2])
    contrasts = np.angle(phases[1:] * phases[0].conj()) / (2 * np.pi * frequency) * 1e6
    result = {
        "ell1_boxes": boxes,
        "J1745_historical": historical_predictions(),
        "J1400_lofar": lofar_diagnostics(),
        "J1400_lwa_contrasts_us": contrasts.tolist(),
    }
    (OUTPUT / "results.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
