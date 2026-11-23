"""Per-sample noise weights for the reduced (heterodyned) J0437 strain.

Weighting every sample of a five-day file by the inverse of that file's plain
variance is far from optimal.  Three properties of the reduced data are the
reason, and this module addresses each:

1.  **Transient before each stride end.**  CWInPy reads and filters the data
    one stride (3600 s) at a time: forward, then backward.  The backward pass
    starts from the forward pass's final filter state, which is the wrong
    state for the reversed series, so the last ~30 s of every stride carry a
    start-up transient thousands of times the noise power.  Strides end
    ``stride`` seconds after the start of each segment and at the end of the
    segment; CWInPy's edge crop does not remove them.  The samples in the
    ``stride_end_seconds`` before each stride end are dropped.

2.  **Changing noise level and short glitches.**  The noise variance is
    estimated from a running median of the sample power, which follows the
    level over minutes and is insensitive to glitches.  Samples whose power
    exceeds ``power_cut`` times the running variance are dropped, with
    ``cut_margin_seconds`` either side.

3.  **Filtered noise is not white.**  The 1-Hz samples have passed a 0.25-Hz
    low-pass filter.  The noise power spectral density at the signal frequency
    is therefore the sample variance divided by the filter's equivalent noise
    bandwidth (0.446 Hz), not by the 1-Hz sample rate.  Using it makes
    ``2|Y|^2/NORM`` follow a chi-squared distribution with two degrees of
    freedom in Gaussian noise.

The weights depend only on the data's power, not its phase, and a signal at
the amplitudes of interest changes a sample's power by a part in a million.
"""
from __future__ import annotations

import numpy as np

LN2 = float(np.log(2.0))

PARAMETERS = dict(stride=3600.0, crop=60.0, stride_end_seconds=60.0,
                  median_seconds=301, power_cut=20.0, cut_margin_seconds=2,
                  min_segment_samples=300)


def equivalent_noise_bandwidth(*, cutoff_hz=0.25, butterworth_order=9,
                               averaging_seconds=1.0, sample_rate_hz=1.0):
    """Equivalent noise bandwidth of the reduction, in Hz.

    The reduction's power response integrated over the sampled band.  White
    input noise of power spectral density S leaves samples of variance
    S times this bandwidth.
    """
    from scipy.integrate import quad
    from .upper_limit import heterodyne_statistic_response

    half = 0.5 * float(sample_rate_hz)
    area, _ = quad(lambda f: heterodyne_statistic_response(
        f, cutoff_hz=cutoff_hz, butterworth_order=butterworth_order,
        averaging_seconds=averaging_seconds), -half, half, limit=400)
    return float(area)


def segments(times, dt=1.0):
    """Start and stop indices of runs of samples spaced by ``dt``."""
    times = np.asarray(times, dtype=np.float64)
    if times.size == 0:
        return np.array([], int), np.array([], int)
    breaks = np.flatnonzero(np.diff(times) > 1.5 * dt) + 1
    return np.r_[0, breaks], np.r_[breaks, times.size]


def stride_end_mask(times, *, stride=3600.0, crop=60.0, stride_end_seconds=60.0,
                    dt=1.0):
    """True for samples to keep; False just before a stride end.

    Each segment has ``crop`` seconds already removed at both ends.  Counting
    samples from the start of the segment, strides end at ``k*stride - crop``
    for every positive integer ``k`` that falls inside it, and at its end.
    """
    keep = np.ones(len(times), dtype=bool)
    window = int(round(stride_end_seconds / dt))
    first = int(round((stride - crop) / dt))
    step = int(round(stride / dt))
    if window <= 0:
        return keep
    for lo, hi in zip(*segments(times, dt)):
        n = hi - lo
        for end in [n] + list(range(first, n, step)):
            keep[lo + max(0, end - window):lo + end] = False
    return keep


def running_median_weights(times, data, *, stride=3600.0, crop=60.0,
                           stride_end_seconds=60.0, median_seconds=301,
                           power_cut=20.0, cut_margin_seconds=2,
                           min_segment_samples=300, dt=1.0, bandwidth_hz=None):
    """Sample spacing over noise power spectral density; zero where dropped.

    Returns ``(weights, summary)``.  ``weights[t] = bandwidth*dt/variance[t]``
    where ``variance`` is the running-median estimate of the complex sample
    variance and ``bandwidth`` the equivalent noise bandwidth.  ``summary``
    counts what was dropped and why.
    """
    from scipy.ndimage import median_filter

    data = np.asarray(data)
    power = data.real ** 2 + data.imag ** 2
    if bandwidth_hz is None:
        bandwidth_hz = equivalent_noise_bandwidth()
    scale = float(bandwidth_hz) * dt
    keep = stride_end_mask(times, stride=stride, crop=crop,
                           stride_end_seconds=stride_end_seconds, dt=dt)
    n_stride_ends = int(np.sum(~keep))
    variance = np.zeros(power.size)
    n_short = 0
    for lo, hi in zip(*segments(times, dt)):
        index = lo + np.flatnonzero(keep[lo:hi])
        if index.size < min_segment_samples:
            n_short += index.size
            keep[lo:hi] = False
            continue
        variance[index] = median_filter(
            power[index], size=int(round(median_seconds / dt)),
            mode="reflect") / LN2
    ratio = np.zeros(power.size)
    np.divide(power, variance, out=ratio, where=keep & (variance > 0))
    keep &= variance > 0
    above = keep & (ratio > power_cut)
    margin = int(round(cut_margin_seconds / dt))
    if margin > 0 and above.any():
        above = np.convolve(above, np.ones(2 * margin + 1), "same") > 0
    n_cut = int(np.sum(keep & above))
    keep &= ~above
    weights = np.zeros(power.size)
    weights[keep] = scale / variance[keep]
    kept = int(keep.sum())
    weighted_noise = float(np.sum(weights ** 2 * power)) if kept else 0.0
    summary = dict(
        samples=int(power.size), kept=kept,
        dropped_stride_ends=n_stride_ends,
        dropped_short_segments=int(n_short),
        dropped_above_power_cut=n_cut,
        # Realised noise power in the weighted sum per unit weight; 1 when
        # the weights are right.
        noise_power_per_unit_weight=(weighted_noise / (scale * np.sum(weights))
                                     if kept else None),
        median_variance=float(np.median(variance[keep])) if kept else None,
        # Signal-to-noise power a constant unit-amplitude signal would get,
        # using the realised noise: (sum w)^2 / sum(w^2 power).
        signal_to_noise_power=(float(np.sum(weights) ** 2 / weighted_noise)
                               if kept else 0.0),
    )
    return weights, summary
