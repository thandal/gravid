import json

import numpy as np
import pytest

from gravfold.orbital_metric import (EPOCH_GPS, J0437_FIT_COVARIANCE,
                                     load_orbital_covariance, covariance_sha256,
                                     require_covariance_match)


def test_j0437_fit_covariance_is_canonical_and_positive_definite():
    cov, metadata = load_orbital_covariance(J0437_FIT_COVARIANCE)
    assert cov.shape == (5, 5)
    assert metadata["target"] == "J0437-4715"
    assert metadata["epoch_convention"] == EPOCH_GPS
    assert np.allclose(cov, cov.T)
    np.linalg.cholesky(cov)


def test_j0437_fit_covariance_preserves_native_qr_covariance():
    cov, metadata = load_orbital_covariance(J0437_FIT_COVARIANCE)
    sigma = np.sqrt(np.diag(cov))
    assert metadata["fit"]["arguments"] == ["-qrfit", "-norescale"]
    raw = np.asarray(metadata["covariance_tcb_internal"])
    scale = metadata["tcb_to_tdb_linear_scale"]
    units = np.array([86400*scale, scale, 1, 1, 86400*scale])
    assert np.allclose(cov / np.outer(units, units), raw, rtol=1e-12, atol=0)
    # Regression against the native fitter, independently checked by QR.
    assert sigma[0] == pytest.approx(0.01586, rel=0.001)
    assert metadata["correlation"][0][3] == pytest.approx(-0.331384, abs=1e-5)
    assert metadata["correlation"][1][2] == pytest.approx(0.0625852, abs=1e-5)
    corr = cov / np.outer(sigma, sigma)
    # This near-degeneracy is the central reason a diagonal covariance is not
    # acceptable for J0437's almost circular orbit.
    assert corr[3, 4] > 0.9999999999
    assert np.allclose(corr, metadata["correlation"], rtol=2e-12, atol=2e-12)


def test_search_extension_requires_exact_covariance_identity():
    require_covariance_match({"covariance_sha256": covariance_sha256()})
    for historical in ({}, {"covariance_sha256": "55b5d988e2e6571517961914c4e17ad2bfb801e282a07a279523b7556cad7bf6"}):
        with pytest.raises(ValueError, match="rerun the search"):
            require_covariance_match(historical)


def test_covariance_loader_rejects_asymmetric_small_cross_term(tmp_path):
    _, metadata = load_orbital_covariance(J0437_FIT_COVARIANCE)
    metadata["covariance"][2][3] += 1e-20
    path = tmp_path / "asymmetric.json"
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="not symmetric"):
        load_orbital_covariance(path)
