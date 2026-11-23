"""Small, dependency-light utilities for injection efficiency and h0 limits."""
from __future__ import annotations

import numpy as np


def heterodyne_statistic_response(offset_hz, *, cutoff_hz=0.25,
                                  butterworth_order=9,
                                  averaging_seconds=1.0):
    """Power-statistic response of the J0437 heterodyne reduction.

    The raw complex heterodyne is passed through a Butterworth filter in both
    time directions and then block-averaged.  For a monochromatic residual at
    ``offset_hz``, the forward/backward Butterworth contributes

        [1 + (|df| / fc) ** (2 n)] ** -2

    to a squared-amplitude statistic, while an averaging window of duration
    ``tau`` contributes ``sinc(df * tau) ** 2``.  The result is normalized to
    unity at zero residual frequency and can therefore put off-source control
    statistics on the target-frequency scale.
    """
    offset = np.asarray(offset_hz, dtype=np.float64)
    cutoff = float(cutoff_hz)
    order = int(butterworth_order)
    duration = float(averaging_seconds)
    if cutoff <= 0 or order <= 0 or duration <= 0:
        raise ValueError("cutoff, Butterworth order, and averaging time must be positive")
    butterworth = 1.0 / (1.0 + (np.abs(offset) / cutoff) ** (2 * order)) ** 2
    averaging = np.sinc(offset * duration) ** 2
    response = butterworth * averaging
    return float(response) if response.ndim == 0 else response


def linear_grid_crossing(x, y, target=0.95):
    """First linearly interpolated crossing on a precomputed efficiency grid."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 2:
        raise ValueError("x and y must be matching one-dimensional grids")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing")
    for i in range(len(x) - 1):
        if y[i] < target <= y[i + 1]:
            fraction = (target - y[i]) / (y[i + 1] - y[i])
            return float(x[i] + fraction * (x[i + 1] - x[i]))
    return None


def paired_grid_bootstrap_crossing(x, detected, *, target=0.95,
                                   confidence=0.90, n_bootstrap=20000,
                                   seed=20260825,
                                   min_valid_fraction=0.90):
    """Bootstrap a crossing when the same injections span every amplitude.

    ``detected`` has shape ``(n_amplitudes, n_injections)``.  Resampling its
    columns preserves the identity of a simulated source across amplitudes,
    as required by the J0437 campaign.  The interval quantifies Monte Carlo
    sampling of the injected source population conditional on the fixed data
    realization and fixed detection threshold; it is not a calibration or
    repeated-background interval.
    """
    x = np.asarray(x, dtype=np.float64)
    detected = np.asarray(detected, dtype=bool)
    if detected.ndim != 2 or detected.shape[0] != len(x):
        raise ValueError("detected must have shape (len(x), n_injections)")
    if detected.shape[1] < 2:
        raise ValueError("at least two paired injections are required")
    confidence = float(confidence)
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between zero and one")
    n_bootstrap = int(n_bootstrap)
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be positive")

    efficiency = detected.mean(axis=1)
    point = linear_grid_crossing(x, efficiency, target)
    rng = np.random.default_rng(seed)
    draws = []
    n_injections = detected.shape[1]
    for _ in range(n_bootstrap):
        take = rng.integers(0, n_injections, n_injections)
        crossing = linear_grid_crossing(
            x, detected[:, take].mean(axis=1), target)
        if crossing is not None:
            draws.append(crossing)
    valid_fraction = len(draws) / n_bootstrap
    enough = (len(draws) >= 20
              and valid_fraction >= float(min_valid_fraction))
    alpha = (1.0 - confidence) / 2.0
    interval = (None if point is None or not enough else
                [float(value) for value in
                 np.quantile(draws, [alpha, 1.0 - alpha])])
    return {
        "point": point,
        "interval": interval,
        "confidence": confidence,
        "target_efficiency": float(target),
        "n_bootstrap": n_bootstrap,
        "seed": int(seed),
        "valid_crossings": len(draws),
        "valid_fraction": valid_fraction,
        "efficiency": efficiency,
    }


def _pava(y, weight):
    """Weighted nondecreasing pool-adjacent-violators fit."""
    y = np.asarray(y, dtype=np.float64)
    weight = np.asarray(weight, dtype=np.float64)
    if y.ndim != 1 or weight.shape != y.shape or np.any(weight <= 0):
        raise ValueError("PAVA needs matching 1-D values and positive weights")
    levels, weights, starts, stops = [], [], [], []
    for i, (value, w) in enumerate(zip(y, weight)):
        levels.append(float(value)); weights.append(float(w))
        starts.append(i); stops.append(i + 1)
        while len(levels) >= 2 and levels[-2] > levels[-1]:
            total = weights[-2] + weights[-1]
            mean = (levels[-2] * weights[-2] + levels[-1] * weights[-1]) / total
            levels[-2:] = [mean]; weights[-2:] = [total]
            stops[-2:] = [stops[-1]]; starts.pop()
    out = np.empty_like(y)
    for level, start, stop in zip(levels, starts, stops):
        out[start:stop] = level
    return out


def binned_isotonic_efficiency(h0, detected, *, n_bins=10, edges=None,
                               weights=None):
    """Isotonic efficiency curve; optional per-trial importance weights.

    Weights support reweighting one injection population to a different
    target posterior (self-normalized within each bin): efficiency becomes
    sum(w * detected) / sum(w).  With weights=None every trial counts once.
    """
    h0 = np.asarray(h0, dtype=np.float64)
    detected = np.asarray(detected, dtype=bool)
    if h0.ndim != 1 or detected.shape != h0.shape or h0.size < 2:
        raise ValueError("h0/detected must be matching vectors with >=2 trials")
    if np.any(~np.isfinite(h0)) or np.any(h0 <= 0):
        raise ValueError("h0 values must be finite and positive")
    if weights is None:
        weights = np.ones_like(h0)
    else:
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != h0.shape or np.any(~np.isfinite(weights)) \
                or np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError("weights must be finite, nonnegative, not all 0")
    if edges is None:
        edges = np.geomspace(h0.min() * (1 - 1e-12),
                             h0.max() * (1 + 1e-12), int(n_bins) + 1)
    else:
        edges = np.asarray(edges, dtype=np.float64)
    index = np.clip(np.searchsorted(edges, h0, side="right") - 1,
                    0, len(edges) - 2)
    centers, counts, recovered, raw = [], [], [], []
    for i in range(len(edges) - 1):
        take = (index == i) & (weights > 0)
        if not np.any(take):
            continue
        w = weights[take]
        centers.append(float(np.exp(
            np.average(np.log(h0[take]), weights=w))))
        counts.append(float(w.sum()))
        recovered.append(float(np.sum(w * detected[take])))
        raw.append(float(np.sum(w * detected[take]) / w.sum()))
    counts = np.asarray(counts, dtype=np.float64)
    raw = np.asarray(raw, dtype=np.float64)
    fitted = _pava(raw, counts)
    return {
        "h0": np.asarray(centers), "count": counts,
        "recovered": np.asarray(recovered, dtype=np.float64),
        "raw_efficiency": raw, "isotonic_efficiency": fitted,
        "edges": edges,
    }


def efficiency_crossing(curve, target=0.95):
    """Log-interpolated target crossing, requiring a two-sided bracket.

    Returns None when the curve never reaches the target OR when the first
    bin is already at/above it: without a sub-target bin below the crossing
    there is no lower bracket, and reporting the first bin center would
    understate the crossing by an unquantified amount.
    """
    x = np.asarray(curve["h0"], dtype=np.float64)
    y = np.asarray(curve["isotonic_efficiency"], dtype=np.float64)
    above = np.flatnonzero(y >= float(target))
    if not len(above) or int(above[0]) == 0:
        return None
    i = int(above[0])
    if y[i] == y[i - 1]:
        return float(x[i])
    fraction = (float(target) - y[i - 1]) / (y[i] - y[i - 1])
    return float(np.exp(np.log(x[i - 1])
                        + fraction * (np.log(x[i]) - np.log(x[i - 1]))))


def bootstrap_crossing(h0, detected, *, target=0.95, n_bins=10,
                       n_bootstrap=1000, seed=0, min_valid_fraction=0.9,
                       weights=None):
    """Point crossing plus a 90% bootstrap interval with failure accounting.

    Replicates whose resampled curve has no bracketed crossing are failures,
    not discards: the interval is reported only when at least
    ``min_valid_fraction`` of replicates cross.  Conditioning an interval on
    the survivors alone biases boundary estimates.  Optional per-trial
    weights (importance reweighting to a different target posterior) ride
    along through the bootstrap resample.
    """
    h0 = np.asarray(h0, dtype=np.float64)
    detected = np.asarray(detected, dtype=bool)
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float64)
    base = binned_isotonic_efficiency(h0, detected, n_bins=n_bins,
                                      weights=weights)
    point = efficiency_crossing(base, target)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(int(n_bootstrap)):
        take = rng.integers(0, len(h0), len(h0))
        curve = binned_isotonic_efficiency(
            h0[take], detected[take], edges=base["edges"],
            weights=None if weights is None else weights[take])
        value = efficiency_crossing(curve, target)
        if value is not None:
            draws.append(value)
    enough = len(draws) >= max(20, min_valid_fraction * int(n_bootstrap))
    interval = (None if point is None or not enough else
                [float(x) for x in np.quantile(draws, [0.05, 0.95])])
    return point, interval, len(draws), base
