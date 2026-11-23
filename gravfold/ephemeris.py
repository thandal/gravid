"""Pulsar ephemeris container and a minimal TEMPO2 .par reader.

The orbital covariance returned here is a FIRST-CUT DIAGONAL built from the
.par 1-sigma error columns. A correct *correlated* 5x5 covariance (which is what
ultimately matters for marginalisation) must come from the TEMPO2 fit covariance
via libstempo (`libstempo.tempopulsar(...).fitterrs`/covariance) or a PINT
posterior; this module does not provide it.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from .orbit import PARAM_NAMES

# Approximate MJD->GPS seconds. Ignores the (slowly varying) leap-second offset
# (~ -18 s in 2024); irrelevant for *relative* orbital phase, which is all the
# synthetic demodulator uses. The real cwinpy path uses exact barycentring.
_MJD_GPS_OFFSET_DAYS = 44244.0


def mjd_to_gps(mjd: float) -> float:
    return (float(mjd) - _MJD_GPS_OFFSET_DAYS) * 86400.0


@dataclass
class OrbitalEphemeris:
    name: str
    f0: float                       # spin frequency [Hz]
    pb: float                       # orbital period [s]
    x: float                        # a_p sin i / c [light-s]
    e: float                        # eccentricity
    om: float                       # argument of periapsis [rad]
    t0: float                       # epoch of periastron [s]
    pb_err: float = 0.0
    x_err: float = 0.0
    e_err: float = 0.0
    om_err: float = 0.0
    t0_err: float = 0.0
    f0_err: float = 0.0
    dist_kpc: float = float("nan")
    ra: str = ""
    dec: str = ""
    model: str = "DD"
    errs_are_illustrative: bool = False

    @property
    def f_gw(self) -> float:
        """GW frequency for the l=m=2 mass-quadrupole mode."""
        return 2.0 * self.f0

    @property
    def params(self) -> np.ndarray:
        return np.array([self.pb, self.x, self.e, self.om, self.t0], dtype=float)

    @property
    def param_errs(self) -> np.ndarray:
        return np.array([self.pb_err, self.x_err, self.e_err,
                         self.om_err, self.t0_err], dtype=float)

    def covariance(self) -> np.ndarray:
        """Diagonal covariance from 1-sigma errors (see module docstring caveat)."""
        v = np.maximum(self.param_errs, 0.0) ** 2
        return np.diag(v)


def _value_and_error(tokens):
    """Parse 'value [fitflag] [error]' from TEMPO2 par tokens after the key."""
    def f(s):
        return float(s.replace("D", "E").replace("d", "E"))
    try:
        value = f(tokens[0])
    except (ValueError, IndexError):
        return None, None
    rest = tokens[1:]
    err = None
    if len(rest) >= 2 and rest[0] in ("0", "1"):
        try:
            err = f(rest[1])
        except ValueError:
            err = None
    elif len(rest) == 1 and rest[0] not in ("0", "1"):
        try:
            err = f(rest[0])
        except ValueError:
            err = None
    return value, err


def parse_par(path: str) -> OrbitalEphemeris:
    """Minimal TEMPO2 .par parser supporting DD ({ECC,OM,T0}) and ELL1
    ({EPS1,EPS2,TASC}) binary models. Good enough for prototyping; for the real
    pipeline cwinpy/lalpulsar ingest .par directly."""
    raw: dict[str, tuple[float, float | None]] = {}
    text: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            key = parts[0].upper()
            val, err = _value_and_error(parts[1:])
            if val is None:
                text[key] = parts[1]
            else:
                raw[key] = (val, err)

    def get(k, default=0.0):
        return raw[k][0] if k in raw else default

    def gerr(k):
        return raw[k][1] if (k in raw and raw[k][1] is not None) else 0.0

    name = text.get("PSRJ", text.get("PSR", "unknown"))
    f0 = get("F0")
    pb = get("PB") * 86400.0
    pb_err = gerr("PB") * 86400.0
    x = get("A1")
    x_err = gerr("A1")

    if "ECC" in raw or "E" in raw:  # DD-like
        ek = "ECC" if "ECC" in raw else "E"
        e = get(ek)
        e_err = gerr(ek)
        om = np.deg2rad(get("OM"))
        om_err = np.deg2rad(gerr("OM"))
        t0 = mjd_to_gps(get("T0"))
        t0_err = gerr("T0") * 86400.0
        model = "DD"
    else:  # ELL1
        eps1, eps2 = get("EPS1"), get("EPS2")
        e = float(np.hypot(eps1, eps2))
        om = float(np.arctan2(eps1, eps2)) if e > 0 else 0.0
        # propagate eps errors -> (e, om) approximately (diagonal)
        de1, de2 = gerr("EPS1"), gerr("EPS2")
        e_err = float(np.hypot(eps1 * de1, eps2 * de2) / e) if e > 0 else max(de1, de2)
        om_err = float(np.hypot(de1, de2) / e) if e > 0 else 0.0
        tasc = mjd_to_gps(get("TASC"))
        t0 = tasc + om / (2.0 * np.pi) * pb
        t0_err = gerr("TASC") * 86400.0
        model = "ELL1"

    dist = float("nan")
    if "PX" in raw and raw["PX"][0] > 0:
        dist = 1.0 / raw["PX"][0]  # parallax [mas] -> distance [kpc]

    return OrbitalEphemeris(
        name=name, f0=f0, pb=pb, x=x, e=e, om=om, t0=t0,
        pb_err=pb_err, x_err=x_err, e_err=e_err, om_err=om_err, t0_err=t0_err,
        f0_err=gerr("F0"), dist_kpc=dist,
        ra=text.get("RAJ", ""), dec=text.get("DECJ", ""), model=model,
    )
