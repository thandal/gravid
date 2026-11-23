"""Noise weights for the reduced J0437 strain: masks, variance, normalisation."""
import numpy as np
import pytest

from gravfold.noise_weights import (equivalent_noise_bandwidth,
                                    running_median_weights, segments,
                                    stride_end_mask)


def _times(lengths, gap=121.0, start=1000.5):
    out, t = [], start
    for n in lengths:
        out.append(t + np.arange(n))
        t = out[-1][-1] + gap
    return np.concatenate(out)


def _noise(rng, n, variance=1.0):
    return np.sqrt(variance / 2) * (rng.normal(size=n) + 1j * rng.normal(size=n))


def test_equivalent_noise_bandwidth_matches_reduction_response():
    # 9th-order forward/backward Butterworth at 0.25 Hz times a 1-s average.
    assert equivalent_noise_bandwidth() == pytest.approx(0.44607, abs=2e-5)
    # A much wider filter leaves only the 1-s averaging response.
    assert equivalent_noise_bandwidth(cutoff_hz=50.0) == pytest.approx(0.7737, rel=2e-3)


def test_segments_split_on_gaps():
    lo, hi = segments(_times([5, 7]))
    assert lo.tolist() == [0, 5] and hi.tolist() == [5, 12]


def test_stride_end_mask_positions():
    # A full segment: 4096 s cropped by 60 s at both ends is 3976 samples.
    # Strides end 3540 samples in (3600 s from the uncropped start) and at
    # the end of the segment.
    keep = stride_end_mask(_times([3976, 1000, 40]))
    full, short, tiny = keep[:3976], keep[3976:4976], keep[4976:]
    dropped = np.flatnonzero(~full)
    assert dropped.tolist() == list(range(3480, 3540)) + list(range(3916, 3976))
    assert np.flatnonzero(~short).tolist() == list(range(940, 1000))
    assert not tiny.any()


def test_weights_are_inverse_noise_density_in_stationary_noise():
    rng = np.random.default_rng(1)
    times = _times([3976] * 6)
    data = _noise(rng, times.size, variance=4.0)
    weights, summary = running_median_weights(times, data, bandwidth_hz=1.0)
    kept = weights > 0
    # 120 of 3976 samples dropped per full segment, and almost nothing else.
    assert summary["dropped_stride_ends"] == 6 * 120
    assert summary["dropped_above_power_cut"] <= 5
    assert np.median(1 / weights[kept]) == pytest.approx(4.0, rel=0.03)
    assert summary["noise_power_per_unit_weight"] == pytest.approx(1.0, abs=0.05)


def test_transients_and_level_changes_do_not_leak_into_the_weighted_sum():
    rng = np.random.default_rng(2)
    times = _times([3976] * 8)
    variance = np.ones(times.size)
    variance[2 * 3976:4 * 3976] = 9.0                 # a noisier stretch
    data = _noise(rng, times.size) * np.sqrt(variance)
    clean = data.copy()
    for segment in range(8):                          # stride-end transients
        for end in (3540, 3976):
            stop = segment * 3976 + end
            data[stop - 25:stop] += 300.0 * np.exp(-np.arange(25)[::-1] / 4.0)
    glitches = rng.choice(times.size, 60, replace=False)
    data[glitches] += 40.0                            # isolated loud samples
    weights, summary = running_median_weights(times, data, bandwidth_hz=1.0)
    assert summary["dropped_above_power_cut"] >= 50
    # The weighted sum's noise matches what clean data would give.
    realised = np.sum(weights ** 2 * np.abs(data) ** 2) / np.sum(weights)
    assert realised == pytest.approx(1.0, abs=0.08)
    # Weights follow the level change: the noisy stretch counts ~9 times less.
    quiet = np.median(weights[weights > 0][:3000])
    noisy = np.median(weights[2 * 3976 + 500:2 * 3976 + 3000])
    assert quiet / noisy == pytest.approx(9.0, rel=0.15)
    # Plain per-file variance would have been dominated by the transients.
    assert np.var(data) / np.var(clean) > 20


def test_statistic_is_chi_squared_with_equivalent_noise_bandwidth():
    # Filter white noise the way the reduction does, then check that the
    # matched-filter statistic with these weights averages 2.
    from scipy.signal import butter, sosfiltfilt
    rng = np.random.default_rng(3)
    rate, count, n = 16, 40, 3976
    sos = butter(9, 0.25, fs=rate, btype="low", output="sos")
    times = _times([n] * count)
    # Templates at small frequency offsets (response within 0.3% of unity),
    # spaced far more than 1/duration apart, give independent noise draws.
    offsets = np.linspace(-0.03, 0.03, 61)
    values = []
    for _ in range(2):
        chunks = []
        for _segment in range(count):
            raw = _noise(rng, (n + 200) * rate)
            filtered = sosfiltfilt(sos, raw)
            averaged = filtered.reshape(-1, rate).mean(axis=1)
            chunks.append(averaged[100:100 + n])
        data = np.concatenate(chunks)
        weights, _ = running_median_weights(times, data)
        for offset in offsets:
            y = np.sum(weights * data * np.exp(-2j * np.pi * offset * times))
            values.append(2 * abs(y) ** 2 / np.sum(weights))
    # Mean of 122 chi-squared(2) draws: 2 +- 0.18.  Dividing the variance by
    # the 1-Hz sample rate instead of the bandwidth would give about 4.5.
    assert 1.5 < np.mean(values) < 2.5
