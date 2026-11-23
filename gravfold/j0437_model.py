"""Full binary-model coordinates for the DE440 J0437 search.

PBDOT is expressed in 1e-12 and OMDOT in degrees/year to avoid artificial
variance floors when placing nodes. LAL conversions are explicit at the edge.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEARCH_PAR = ROOT / "gravfold/resources/j0437_de440_search.par"
COVARIANCE = ROOT / "gravfold/resources/j0437_de440_binary_covariance.json"
KEYS = ("PB", "A1", "ECC", "OM", "T0", "PBDOT", "OMDOT", "KIN", "KOM", "M2")
UNITS = ("s", "light-s", "1", "rad", "s", "1e-12", "deg/yr", "rad", "rad", "solar_mass")


def lal_factors():
    import lal
    return np.array([1, 1, 1, 1, 1, 1e-12, np.pi / (180 * 365.25 * 86400),
                     1, 1, lal.MSUN_SI])


def binary_mean(par):
    return np.array([float(par[key]) for key in KEYS]) / lal_factors()


def apply_binary(par, theta):
    theta = np.asarray(theta, float)
    if theta.shape != (len(KEYS),):
        raise ValueError("full binary model requires all ten correlated parameters")
    for key, value in zip(KEYS, theta * lal_factors()):
        par[key] = float(value)
    par["SINI"] = float(np.sin(par["KIN"]))
    return par


def load_covariance(path=COVARIANCE):
    artifact = json.loads(Path(path).read_text())
    if artifact["parameter_order"] != list(KEYS) or artifact["canonical_units"] != list(UNITS):
        raise ValueError("wrong full-binary parameter order or units")
    cov = np.asarray(artifact["covariance"], float)
    if cov.shape != (10, 10) or not np.all(np.isfinite(cov)):
        raise ValueError("invalid full binary covariance")
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=1e-30):
        raise ValueError("asymmetric full binary covariance")
    np.linalg.cholesky(cov)
    return cov, artifact


def waveform_metric(parfile, cov, samples, eps=0.1):
    """Fixed-frequency phase metric with sampled detector noise/amplitude weights.

    samples contains (detector, GPS times, inverse-variance quadrature weights).
    Only constant phase is projected out. SSB parameters remain at the fitted
    values; all fitted binary parameters are varied together.
    """
    from cwinpy.parfile import PulsarParameters
    from cwinpy.signal import HeterodynedCWSimulator
    from .cwinpy_adapter import _ensure_template_amplitude
    base = _ensure_template_amplitude(PulsarParameters(str(parfile)))
    mean = binary_mean(base)
    factor = np.linalg.cholesky(cov)
    derivative_blocks, weight_blocks = [], []
    for detector, times, inverse_variance in samples:
        sim = HeterodynedCWSimulator(base, detector, np.asarray(times, float))
        amplitude = np.asarray(sim.model(newpar=base, freqfactor=2.0))
        weight_blocks.append(np.asarray(inverse_variance) * abs(amplitude)**2)
        derivatives = []
        for j in range(len(KEYS)):
            phases = []
            for sign in (1, -1):
                par = apply_binary(PulsarParameters(str(parfile)), mean + sign*eps*factor[:, j])
                phase = sim.model(newpar=par, updateSSB=False, updateBSB=True,
                                  freqfactor=2.0, phase_only=True)
                phases.append(np.exp(-2j*np.pi*phase))
            derivatives.append(np.angle(phases[0] * phases[1].conj()) / (2*eps))
        derivative_blocks.append(np.asarray(derivatives))
    derivatives = np.concatenate(derivative_blocks, axis=1)
    weights = np.concatenate(weight_blocks)
    weights /= weights.sum()
    derivatives -= (derivatives @ weights)[:, None]
    return (derivatives * weights) @ derivatives.T
