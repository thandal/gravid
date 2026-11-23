"""Binary orbit phase model: Kepler solver and the binary Rømer delay.

The Rømer delay R(t)/c is the light-travel-time modulation of the pulsar's
signal caused by its orbital motion about the binary barycentre. It is the
*only* binary term that matters at the v/c << 1, fully-Keplerian level used by
the five-vector resampling demodulator (Amicucci et al. 2025, arXiv:2503.07863,
Eq. 12 — note their Einstein/Shapiro terms are the *solar* ones, not binary).

Canonical orbital parameter vector throughout gravid:
    theta = (pb, x, e, om, t0)
        pb : orbital period            [seconds]
        x  : a_p sin(i) / c            [light-seconds]   (== TEMPO2 "A1")
        e  : eccentricity              [dimensionless]
        om : argument of periapsis     [radians]
        t0 : epoch of periastron       [seconds]         (same clock as t)
"""
from __future__ import annotations
import numpy as np

PARAM_NAMES = ("pb", "x", "e", "om", "t0")


def kepler_solve(M, e, tol: float = 1e-13, maxiter: int = 50,
                 fixed_iters: int | None = None):
    """Solve Kepler's equation E - e*sin(E) = M (vectorised Newton iteration).

    ``fixed_iters`` runs EXACTLY that many Newton steps and skips the
    convergence test.  **Any caller that processes a series in blocks must use
    it.**  The adaptive path below stops on ``np.max(np.abs(dE)) < tol`` — a
    BLOCK-WIDE reduction — so every sample receives however many iterations
    its worst neighbour needed, and a sample's answer depends on which other
    samples happen to share its block.  Re-block the same series and the last
    bits move: for an orbit with e=0.584, 16.7M-sample blocks took 3
    iterations everywhere while 32k blocks converged in 2 for 115 of 6,104
    blocks, changing the injected strain.

    Iteration budget at tol=1e-13 from that same sweep: e=0.584 needs 4,
    e=0.7 needs 5, e=0.9 needs 7.  K=5 is at the residual floor (8.9e-16) for
    everything up to e~0.7.
    """
    M = np.asarray(M, dtype=float)
    Mr = np.mod(M, 2.0 * np.pi)
    # Good starting guess; for the near-circular MSP orbits here e < 1e-4 so this
    # converges in 2-3 iterations.
    E = Mr + e * np.sin(Mr) / (1.0 - np.sin(Mr + e) + np.sin(Mr))
    E = np.where(np.isfinite(E), E, Mr)
    if fixed_iters is not None:
        # Same arithmetic as the adaptive branch, minus the reduction that
        # couples samples to their block-mates.
        for _ in range(int(fixed_iters)):
            f = E - e * np.sin(E) - Mr
            fp = 1.0 - e * np.cos(E)
            E = E - f / fp
        return E
    for _ in range(maxiter):
        f = E - e * np.sin(E) - Mr
        fp = 1.0 - e * np.cos(E)
        dE = f / fp
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def roemer_delay(t, pb, x, e, om, t0, fixed_iters: int | None = None):
    """Binary Rømer delay R(t)/c, in seconds (Damour-Deruelle / Blandford-Teukolsky).

    R(t)/c = x * [ sin(om)*(cos E - e) + cos(om)*sin E*sqrt(1-e^2) ]
    with the eccentric anomaly E from M = 2*pi*(t-t0)/pb.
    """
    t = np.asarray(t, dtype=float)
    M = 2.0 * np.pi * (t - t0) / pb
    E = kepler_solve(M, e, fixed_iters=fixed_iters)
    return x * (np.sin(om) * (np.cos(E) - e)
                + np.cos(om) * np.sin(E) * np.sqrt(1.0 - e * e))


def bsb_delay(t, pb, x, e, om, t0, fixed_iters: int | None = None):
    """The BINARY DELAY in the pulsar-timing (lalpulsar/cwinpy) sense: the
    quantity that ADDS to the barycentric arrival time to give the source's
    emission proper time,

        tau_emission = t_SSB + bsb(t),          bsb = -R

    where R is the classical Roemer delay above (positive when the pulsar is
    on the far side of its orbit, so its light arrives LATER — hence the
    emission time is EARLIER by R).

    THIS is the delay that enters a signal's phase: a CW at f_gw observed at
    barycentric time t has phase 2*pi*f_gw*(t + bsb(t)). Everything that
    models or demodulates a signal must use this, not R.

    This is the sign convention of lalpulsar's own binary delay.  The sign
    matters: for months of coherent integration a flipped orbital delay does
    not merely dephase the template, it DOUBLES the modulation it was
    supposed to remove.
    """
    return -roemer_delay(t, pb, x, e, om, t0, fixed_iters=fixed_iters)


def orbital_phase(t, f_gw, theta):
    """GW phase modulation (radians) imprinted by the binary orbit on a signal
    at f_gw: 2*pi*f_gw*bsb(t) = -2*pi*f_gw*R(t)/c. See bsb_delay for why the
    sign is what it is."""
    pb, x, e, om, t0 = theta
    return 2.0 * np.pi * f_gw * bsb_delay(t, pb, x, e, om, t0)
