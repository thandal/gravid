import json
from pathlib import Path

import numpy as np
import pytest

from gravfold.j0437_model import COVARIANCE, SEARCH_PAR, load_covariance


def test_expanded_covariance_retains_period_advance_correlation():
    cov, artifact = load_covariance()
    corr = cov[0,6]/np.sqrt(cov[0,0]*cov[6,6])
    assert corr > .999999
    # Derived sidereal angular frequency includes both terms with opposite signs.
    # Fixing OMDOT destroys this cancellation even with the correct PB sigma.
    from gravfold.ephemeris import parse_par
    pb = parse_par(SEARCH_PAR).params[0]
    gradient = np.zeros(10)
    gradient[0] = -2*np.pi/pb**2
    gradient[6] = np.pi/(180*365.25*86400)
    joint = gradient @ cov @ gradient
    fixed = gradient[0]**2*cov[0,0]
    assert joint < fixed*1e-5
    assert artifact['independent_verification']['max_correlation_difference'] < 1e-7


def test_lal_binary_coordinates_roundtrip_and_geometry():
    pytest.importorskip('cwinpy')
    from cwinpy.parfile import PulsarParameters
    from gravfold.j0437_model import binary_mean, apply_binary, KEYS
    p = PulsarParameters(str(SEARCH_PAR))
    theta = binary_mean(p)
    q = apply_binary(PulsarParameters(str(SEARCH_PAR)), theta)
    for key in KEYS:
        assert q[key] == pytest.approx(p[key], rel=2e-15)
    assert q['SINI'] == pytest.approx(np.sin(p['KIN']))
    with pytest.raises(ValueError):
        apply_binary(q, theta[:5])


def test_campaign_loader_matches_cwinpy():
    pytest.importorskip('cwinpy')
    from cwinpy import HeterodynedData
    from gravfold.j0437_campaign import load_chunk
    paths = sorted(Path('results/campaign_het_v2').glob('het_*_g*.hdf5'))
    if not paths:
        pytest.skip('restored campaign unavailable')
    for detector in ('H1','L1'):
        path = next(p for p in paths if f'_{detector}_' in p.name)
        det, times, data, _ = load_chunk(path)
        expected = HeterodynedData.read(str(path))
        assert det == expected.detector
        assert np.array_equal(times, np.asarray(expected.times))
        assert np.array_equal(data, np.asarray(expected.data))


def test_search_astrometric_errors_preserve_native_units():
    source = Path('scripts/timing/qr_solve/fit.par')
    rows = lambda p: {l.split()[0]: l.split()[1:] for l in p.read_text().splitlines() if l.split()}
    original, converted = rows(source), rows(SEARCH_PAR)
    for key in ('RAJ', 'DECJ'):
        assert converted[key] == original[key]
