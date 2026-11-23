"""Demodulators: map an orbital-parameter vector -> a detection statistic.

`SyntheticDemodulator` is a faithful *miniature* of the real five-vector /
heterodyne resampling: it evaluates the binary Rømer phase over the time series
(Kepler solve + trig) and forms a complex matched filter. Crucially, its
per-call cost has the SAME structure as the real per-sample demodulation cost
(O(N_samples) phase evaluation + reduction), so the benchmark numbers it
produces are meaningful for the compute plan. It needs no GW data, no cwinpy,
and no LAL — so the whole marginalisation logic is testable in a sandbox.

The real path is `gravfold.cwinpy_adapter.CwinpyDemodulator`.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .orbit import bsb_delay


@dataclass
class DemodResult:
    twoF: float          # 2F detection statistic (2 d.o.f.; ~chi^2_2 under H0)
    complex_sum: complex # matched-filter complex amplitude (for combining)
    n: int


class SyntheticDemodulator:
    """Matched filter on heterodyned complex data, residual binary phase only.

    data[t] models a heterodyned series where Earth motion and spin evolution
    have already been removed, leaving the binary-orbit phase modulation
    exp(+i 2*pi*f_gw*b(t)) with b = -R the pulsar-timing binary delay (see
    gravfold.orbit.bsb_delay -- the sign a real signal carries; using +R here
    would double the modulation on real data, invisibly to any self-consistent
    test). Demodulating with a trial orbit cancels that phase when (and only
    when) the trial matches truth; a wrong orbit decorrelates the sum, which is
    exactly the mechanism behind the T_obs^max exclusion criterion.
    """

    def __init__(self, t, data, f_gw, sigma=1.0):
        self.t = np.asarray(t, dtype=float)
        self.data = np.asarray(data, dtype=complex)
        self.f_gw = float(f_gw)
        self.sigma = float(sigma)
        self.n = self.t.size

    def __call__(self, theta) -> DemodResult:
        pb, x, e, om, t0 = theta
        dly = bsb_delay(self.t, pb, x, e, om, t0)
        template = np.exp(2j * np.pi * self.f_gw * dly)
        # vdot(template, data) = sum(conj(template) * data); at theta=truth the
        # conj(+i phase) cancels the data's +i phase -> coherent sum.
        y = np.vdot(template, self.data)
        twoF = 2.0 * (np.abs(y) ** 2) / (self.n * self.sigma ** 2)
        return DemodResult(twoF=float(twoF), complex_sum=complex(y), n=self.n)


def make_synthetic_data(eph, true_theta, tobs_days, fs=1.0 / 60.0,
                        h0=1.0, sigma=1.0, rng=None):
    """Generate a heterodyned series with an injected CW at `true_theta`.

    fs defaults to 1/60 Hz, matching cwinpy's binned heterodyne output rate.
    Returns (t, data) with t in seconds from the start of the span.

    WARNING: h0 defaults to 1.0 (per-sample SNR 1) -- a "noise-only" series
    requires h0=0 EXPLICITLY. A forgotten default leaves a decorrelated
    injection that inflates the node-2F mean above the chi^2_2 value of 2.
    """
    rng = np.random.default_rng() if rng is None else rng
    n = int(round(tobs_days * 86400.0 * fs))
    t = np.arange(n) / fs
    dly = bsb_delay(t, *np.asarray(true_theta, dtype=float))   # = -R
    signal = h0 * np.exp(2j * np.pi * eph.f_gw * dly)
    noise = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * sigma / np.sqrt(2.0)
    return t, signal + noise


def center_epoch(eph, tobs_days):
    """Return a copy of the central params with the periastron epoch t0 placed at
    mid-span, so the synthetic series actually samples the orbit."""
    theta = eph.params.copy()
    theta[4] = 0.5 * tobs_days * 86400.0
    return theta
