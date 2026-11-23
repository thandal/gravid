import numpy as np
from gravfold.orbital_metric import (whitened_reduced_metric,
                                     worst_case_mismatch, tobs_max)

# eccentric test point (e=0.1): om is well-defined, so a diagonal covariance is
# not degenerate -- near-circular orbits need the full covariance (see
# theta_cov_from_parfile), which is exercised on real data, not here.
THETA = np.array([5.0 * 86400, 3.3, 0.1, 0.6, 3.0e8])    # pb,x,e,om,t0
COV = np.diag([0.3, 3e-5, 1e-5, 1e-4, 1.0]) ** 2         # illustrative 1-sigma


def test_reduced_metric_is_symmetric_psd():
    g = whitened_reduced_metric(THETA, COV, 350.0, 180.0 * 86400)
    assert g.shape == (5, 5)
    assert np.allclose(g, g.T, atol=1e-12 * (np.abs(g).max() + 1))
    w = np.linalg.eigvalsh(g)          # Gram matrix; PSD up to (severe) conditioning
    assert w.max() > 0
    assert w.min() > -1e-3 * w.max()


def test_worst_mismatch_grows_with_tobs():
    ms = [worst_case_mismatch(THETA, COV, 350.0, T * 86400)
          for T in (10, 40, 160)]
    assert ms[0] < ms[1] < ms[2]


def test_tobs_max_shrinks_with_bigger_covariance():
    lo, _ = tobs_max(THETA, COV, 350.0)
    hi, _ = tobs_max(THETA, COV * 100.0, 350.0)   # 10x worse 1-sigma
    assert hi < lo


def test_tobs_max_shrinks_with_frequency():
    # higher f_gw -> higher modulation index -> more orbital sensitivity
    lo, _ = tobs_max(THETA, COV, 50.0)
    hi, _ = tobs_max(THETA, COV, 800.0)
    assert hi < lo
