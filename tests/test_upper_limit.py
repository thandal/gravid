import numpy as np

from gravfold.upper_limit import (binned_isotonic_efficiency,
                                   bootstrap_crossing, efficiency_crossing,
                                   heterodyne_statistic_response,
                                   linear_grid_crossing,
                                   paired_grid_bootstrap_crossing)


def test_heterodyne_response_includes_filter_and_block_average():
    response = heterodyne_statistic_response(np.array([0.0, 0.198, 0.2]))
    assert response[0] == 1.0
    assert np.isclose(response[1], 0.8516971020843417)
    assert np.isclose(response[2], 0.8444419359046818)
    assert np.all(np.diff(response) < 0)


def test_paired_grid_bootstrap_resamples_source_identities():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    # Each column is one source evaluated at all four amplitudes.
    detected = np.array([
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1],
        [0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1],
    ], dtype=bool)
    result = paired_grid_bootstrap_crossing(
        x, detected, target=0.75, confidence=0.8,
        n_bootstrap=500, seed=9, min_valid_fraction=0.5)
    assert result["point"] == linear_grid_crossing(
        x, detected.mean(axis=1), 0.75)
    assert result["interval"] is not None
    assert result["interval"][0] <= result["point"] <= result["interval"][1]
    assert result["valid_crossings"] > 250


def test_isotonic_efficiency_is_monotone_and_weighted():
    h0 = np.geomspace(1, 10, 100)
    detected = h0 > 4
    # Deliberately introduce local reversals.
    detected[45:50] = True
    detected[50:55] = False
    curve = binned_isotonic_efficiency(h0, detected, n_bins=10)
    assert np.all(np.diff(curve["isotonic_efficiency"]) >= 0)
    assert curve["count"].sum() == len(h0)


def test_crossing_and_bootstrap_bracket_step():
    rng = np.random.default_rng(5)
    h0 = np.exp(rng.uniform(np.log(1), np.log(10), 500))
    probability = 1 / (1 + np.exp(-8 * (np.log(h0) - np.log(3))))
    detected = rng.uniform(size=len(h0)) < probability
    point, interval, valid, curve = bootstrap_crossing(
        h0, detected, n_bins=12, n_bootstrap=200, seed=8)
    assert point is not None and 3 < point < 8
    assert interval is not None and interval[0] <= point <= interval[1]
    assert valid > 150
    assert efficiency_crossing(curve) == point


def test_unbracketed_returns_none():
    h0 = np.geomspace(1, 10, 100)
    curve = binned_isotonic_efficiency(h0, np.zeros(100, dtype=bool), n_bins=10)
    assert efficiency_crossing(curve) is None


def test_zero_one_weights_equal_subsetting():
    rng = np.random.default_rng(11)
    h0 = np.exp(rng.uniform(np.log(1), np.log(10), 400))
    detected = h0 > 3.5
    keep = rng.uniform(size=400) < 0.5
    weighted = binned_isotonic_efficiency(
        h0, detected, n_bins=8, weights=keep.astype(float))
    subset = binned_isotonic_efficiency(
        h0[keep], detected[keep], n_bins=8,
        edges=weighted["edges"])
    assert np.allclose(weighted["h0"], subset["h0"])
    assert np.allclose(weighted["isotonic_efficiency"],
                       subset["isotonic_efficiency"])


def test_importance_weights_shift_the_crossing():
    # Two subpopulations with different detection curves: reweighting toward
    # the harder one must raise the estimated crossing.
    rng = np.random.default_rng(12)
    h0 = np.exp(rng.uniform(np.log(1), np.log(10), 2000))
    hard = rng.uniform(size=2000) < 0.5
    threshold = np.where(hard, 5.0, 2.0)
    detected = h0 > threshold
    easy_w = np.where(hard, 0.2, 1.8)
    hard_w = np.where(hard, 1.8, 0.2)
    easy_point, _, _, _ = bootstrap_crossing(
        h0, detected, n_bins=12, n_bootstrap=50, seed=1, weights=easy_w)
    hard_point, _, _, _ = bootstrap_crossing(
        h0, detected, n_bins=12, n_bootstrap=50, seed=1, weights=hard_w)
    assert easy_point is not None and hard_point is not None
    assert hard_point > easy_point


def test_first_bin_above_target_has_no_lower_bracket():
    # Every trial detected: the curve starts at 100%, so there is no
    # sub-target bin below the crossing and no defensible point estimate.
    h0 = np.geomspace(1, 10, 100)
    curve = binned_isotonic_efficiency(h0, np.ones(100, dtype=bool), n_bins=10)
    assert efficiency_crossing(curve) is None
    point, interval, valid, _ = bootstrap_crossing(
        h0, np.ones(100, dtype=bool), n_bins=10, n_bootstrap=50, seed=3)
    assert point is None and interval is None and valid == 0
