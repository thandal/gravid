"""Real-data demodulator built on cwinpy (optional dependency).

Replaces SyntheticDemodulator once a heterodyned series is in hand. cwinpy does
the expensive broadband->narrowband heterodyne ONCE per pulsar-detector; the
`HeterodynedCWSimulator.model(newpar, updateSSB=True, updateBSB=True)` re-phases
that cached series for ANY orbit with no raw-data round trip -- which is what
makes per-node marginalisation cheap (REPORT.md sec. 7).

Validated against cwinpy 1.5.0: `HeterodynedCWSimulator` lives in `cwinpy.signal`;
`.model(newpar=..., updateSSB=True, updateBSB=True, freqfactor=2.0)` returns the
complex heterodyned signal template. cwinpy is NOT needed to import gravfold.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .demodulate import DemodResult

_IMPORT_HINT = (
    "CwinpyDemodulator requires cwinpy + lalsuite (`conda install -c conda-forge "
    "cwinpy` into ./.conda). Use SyntheticDemodulator for sandbox/logic work."
)


def _require_cwinpy():
    try:
        import cwinpy  # noqa: F401
        from cwinpy.signal import HeterodynedCWSimulator  # noqa: F401
    except Exception as exc:  # pragma: no cover - only without cwinpy
        raise ImportError(_IMPORT_HINT) from exc


def theta_from_par(par):
    """Extract the canonical orbit theta=(pb[s], x[ls], e, om[rad], t0[s]) from a
    cwinpy PulsarParameters, using cwinpy's OWN (exact) epochs and units.

    Critically this does NOT go through gravfold.ephemeris.parse_par: that uses an
    approximate leap-second-free MJD->GPS which is offset ~tens of seconds from
    cwinpy's exact TASC/T0, enough to decorrelate the matched filter at ~350 Hz.
    Always seed the cwinpy-path posterior from this, not from parse_par.
    """
    pb = float(par["PB"]); x = float(par["A1"])
    if par["EPS1"] is not None:                       # ELL1
        eps1, eps2 = float(par["EPS1"]), float(par["EPS2"])
        e = float(np.hypot(eps1, eps2))
        om = float(np.arctan2(eps1, eps2)) if e > 0 else 0.0
        t0 = float(par["TASC"]) + om / (2.0 * np.pi) * pb
    else:                                             # DD / BT / T2 (cwinpy OM in rad)
        e = float(par["ECC"] or 0.0)
        om = float(par["OM"] or 0.0)
        t0 = float(par["T0"])
    return np.array([pb, x, e, om, t0], dtype=float)


def _ensure_template_amplitude(par):
    """Force a unit-amplitude, fiducial-orientation template. A par without H0
    makes HeterodynedCWSimulator emit an ALL-ZERO model -> 2F silently 0.
    2F = 2|<d,m>|^2/<m,m> is amplitude-invariant, so H0=1 is safe."""
    if par["H0"] is None or float(par["H0"]) == 0.0:
        par["H0"] = 1.0
    for k, v in (("COSIOTA", 0.0), ("PSI", 0.0), ("PHI0", 0.0)):
        if par[k] is None:
            par[k] = v
    return par


def _apply_orbit(par, theta):
    """Write the canonical orbit theta=(pb[s], x[ls], e, om[rad], t0[s]) into a
    cwinpy PulsarParameters, respecting its binary model (ELL1 or DD/BT).

    NB cwinpy stores OM in RADIANS (not degrees) internally."""
    pb, x, e, om, t0 = [float(v) for v in theta]
    par["PB"] = pb
    par["A1"] = x
    if par["EPS1"] is not None or (par["BINARY"] is not None
                                   and str(par["BINARY"]).upper().startswith("ELL1")):
        par["EPS1"] = e * np.sin(om)
        par["EPS2"] = e * np.cos(om)
        par["TASC"] = t0 - om / (2.0 * np.pi) * pb    # inverse of theta_from_par
    else:  # DD / BT / T2
        par["ECC"] = e
        par["OM"] = om                                # radians, matching cwinpy
        par["T0"] = t0
    return par


@dataclass
class CwinpyDemodulator:
    """Wrap a cached heterodyned series as an orbit -> 2F callable.

    Parameters
    ----------
    het : cwinpy.HeterodynedData   (cached coarse series for one pulsar+detector)
    base_par : str                 (path to the pulsar .par; orbit overwritten per call)
    """
    het: object
    base_par: str

    def __post_init__(self):
        _require_cwinpy()
        from cwinpy.parfile import PulsarParameters
        from cwinpy.signal import HeterodynedCWSimulator
        self._PulsarParameters = PulsarParameters
        base = _ensure_template_amplitude(PulsarParameters(self.base_par))
        self._base = base
        self._sim = HeterodynedCWSimulator(par=base, det=self.het.detector,
                                           times=self.het.times)
        self._data = np.asarray(self.het.data, dtype=complex)
        # per-sample noise variance for the noise-weighted matched filter; cwinpy
        # estimates it (het.vars). Fall back to a robust constant if absent.
        vars_ = getattr(self.het, "vars", None)
        if vars_ is None:
            mad = np.median(np.abs(self._data - np.median(self._data)))
            v = (1.4826 * mad) ** 2 if mad > 0 else np.var(self._data)
            self._invvar = np.full(self._data.size, 1.0 / max(v, 1e-300))
        else:
            self._invvar = 1.0 / np.clip(np.asarray(vars_, dtype=float), 1e-300, None)

    @property
    def central_theta(self):
        """Canonical orbit (pb, x, e, om, t0) from cwinpy's exact par -- seed the
        posterior from THIS (not gravfold.ephemeris.parse_par; see theta_from_par)."""
        return theta_from_par(self._base)

    def __call__(self, theta) -> DemodResult:
        par = _ensure_template_amplitude(
            _apply_orbit(self._PulsarParameters(self.base_par), theta))
        model = np.asarray(self._sim.model(newpar=par, updateSSB=True,
                                           updateBSB=True, freqfactor=2.0),
                           dtype=complex)
        if not np.any(model):
            raise RuntimeError("zero template from HeterodynedCWSimulator — "
                               "check par amplitude/orientation fields")
        # noise-weighted matched filter: 2F = 2|<d,m>_w|^2 / <m,m>_w, weights 1/sigma^2
        wm = np.conj(model) * self._invvar
        y = np.sum(wm * self._data)
        norm = float(np.sum((np.abs(model) ** 2) * self._invvar))
        twoF = 2.0 * (np.abs(y) ** 2) / norm if norm > 0 else 0.0
        return DemodResult(twoF=float(twoF), complex_sum=complex(y),
                           n=self._data.size)
