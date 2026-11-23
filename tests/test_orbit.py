import numpy as np
from gravid.orbit import kepler_solve, roemer_delay


def test_kepler_circular():
    M = np.linspace(0, 2 * np.pi, 50, endpoint=False)
    E = kepler_solve(M, 0.0)
    assert np.allclose(np.mod(E, 2 * np.pi), np.mod(M, 2 * np.pi), atol=1e-10)


def test_kepler_roundtrip():
    M = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    for e in (1e-5, 0.1, 0.6, 0.617):  # incl. the Hulse-Taylor eccentricity
        E = kepler_solve(M, e)
        assert np.allclose(E - e * np.sin(E), np.mod(M, 2 * np.pi), atol=1e-10)


def test_roemer_circular_amplitude():
    # e=0: R/c = x*sin(E+om), so peak-to-peak amplitude is 2*x.
    t = np.linspace(0, 5 * 86400.0, 5000)
    R = roemer_delay(t, pb=5.0 * 86400.0, x=3.3, e=0.0, om=0.4, t0=0.0)
    assert np.isclose(R.max() - R.min(), 2 * 3.3, rtol=1e-3)


def test_roemer_eccentric_finite():
    t = np.linspace(0, 80 * 86400.0, 10000)
    R = roemer_delay(t, pb=79.5 * 86400.0, x=9.3, e=0.4, om=1.1, t0=1e6)
    assert np.all(np.isfinite(R))
    assert R.max() - R.min() < 4 * 9.3  # bounded by orbit size
