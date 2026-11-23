import numpy as np
from gravid.posterior import OrbitalPosterior

MEAN = np.array([5.0 * 86400.0, 3.3, 1e-5, 0.0, 3.0e6])
COV = np.diag([0.05, 2e-4, 1e-6, 0.05, 2.0]) ** 2


def test_mc_recovers_mean_cov():
    post = OrbitalPosterior(MEAN, COV)
    pts, w = post.sample_mc(200_000, rng=np.random.default_rng(0))
    assert np.allclose(pts.mean(0), MEAN, rtol=0, atol=5e-2 * np.sqrt(np.diag(COV)) + 1e-9)
    emp = np.cov(pts.T)
    assert np.allclose(np.diag(emp), np.diag(COV), rtol=0.1)


def test_unscented_weights_and_mean():
    post = OrbitalPosterior(MEAN, COV)
    pts, w = post.nodes_unscented()
    assert pts.shape == (11, 5)
    assert np.isclose(w.sum(), 1.0)
    assert np.allclose((w[:, None] * pts).sum(0), MEAN, rtol=0, atol=1e-6 * np.abs(MEAN) + 1e-9)


def test_gauss_hermite_weights_sum_to_one():
    post = OrbitalPosterior(MEAN, COV)
    pts, w = post.nodes_gauss_hermite(k=3)
    assert pts.shape == (3 ** 5, 5)
    assert np.isclose(w.sum(), 1.0, atol=1e-9)
    assert np.allclose((w[:, None] * pts).sum(0), MEAN, rtol=0, atol=1e-6 * np.abs(MEAN) + 1e-9)


def test_zero_variance_axis_is_invertible():
    cov = COV.copy()
    cov[2, 2] = 0.0  # perfectly known eccentricity
    post = OrbitalPosterior(MEAN, cov)
    pts, w = post.nodes_unscented()
    assert np.all(np.isfinite(pts))
