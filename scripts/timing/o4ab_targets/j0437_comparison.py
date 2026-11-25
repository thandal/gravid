"""Extend the companion article's correlated-orbit test to the O4ab grid."""

import itertools
import json

import numpy as np
from astropy.time import Time
from cwinpy.parfile import PulsarParameters
from cwinpy.signal import HeterodynedCWSimulator
from paths import WORK

from gravid.cwinpy_adapter import _ensure_template_amplitude
from gravid.j0437_model import SEARCH_PAR, apply_binary, binary_mean, load_covariance

times = Time(
    np.linspace(60088.625, 60703 + 17 / 24, 16384), format="mjd", scale="utc"
).gps
base = _ensure_template_amplitude(PulsarParameters(str(SEARCH_PAR)))
simulator = HeterodynedCWSimulator(base, "H1", times)
mean = binary_mean(base)
covariance, _ = load_covariance()
root = np.linalg.cholesky(covariance)


def phase(direction):
    par = apply_binary(PulsarParameters(str(SEARCH_PAR)), mean + root @ direction)
    cycles = simulator.model(
        newpar=par, updateSSB=False, updateBSB=True, phase_only=True, freqfactor=2.0
    )
    return np.unwrap(2 * np.pi * np.asarray(cycles, float))


axes = np.array(
    [(phase(axis) - phase(-axis)) / 0.2 for axis in 0.1 * np.eye(len(mean))]
)
axes -= axes.mean(axis=1)[:, None]
_, directions = np.linalg.eigh(axes @ axes.T / len(times))
samples = {
    "gaussian": np.random.default_rng(20261017).normal(size=(256, len(mean))),
    "corners": np.array(list(itertools.product([-3.0, 3.0], repeat=len(mean))))
    @ directions.T,
}
result = {}
for group, offsets in samples.items():
    result[group] = {
        "orbits": len(offsets),
        "largest_uniform_loss": float(
            max(
                1 - abs(np.exp(1j * 1.001 * phase(offset)).mean()) ** 2
                for offset in offsets
            )
        ),
    }
(WORK / "j0437_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
