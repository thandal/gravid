"""Multi-parameter orbital coherence metric and T_obs^max (Leaci-Prix / Mirasola
style), plus a covariance extractor from a .par file.

How long can a fully-coherent demodulation at the *central* orbit run before the
orbital-parameter uncertainty spoils it? The binary orbit imprints a GW phase
modulation phi(t) = 2*pi*f_gw*R(t;theta)/c, theta = (pb, x, e, om, t0). A CW
search maximises over the intrinsic GW frequency and phase, so a constant and a
linear-in-t phase drift are free; the *reduced* mismatch from an orbital offset
d(theta) is the variance of the residual phase.

Two things the implementation gets right (both learned the hard way from real
data):

  * WHITENED perturbations. We differentiate along the columns of L = chol(Sigma)
    with an O(1e-3) step, i.e. steps are a small fraction of 1 sigma -- NOT a
    fraction of the parameter's absolute value. For a high-modulation-index MSP
    (t0 ~ 1e9 s) an absolute-relative step would wrap the phase hundreds of times
    and give nonsense. The whitened metric g_z directly yields the worst-case
    1-sigma mismatch as its largest eigenvalue.

  * FULL covariance, good basis. For near-circular orbits om and the orbit epoch
    are nearly degenerate, so DIAGONAL errors in (e, om, t0) fabricate huge false
    decorrelation. `theta_cov_from_parfile` provides an ELL1-Jacobian proxy when
    only a .par is available. This module's five-parameter metric is a reduced
    model diagnostic. The J0437 production search uses `j0437_model` to
    propagate all ten fitted binary parameters through the full T2 waveform,
    retaining the native QR PB-OMDOT correlation.
"""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import numpy as np

from .orbit import roemer_delay

_VAR_FLOOR = 1e-30

# ---------------------------------------------------------------------------
# Epoch conventions. theta[4] (the orbit epoch t0) is in DIFFERENT time units
# depending on which extractor produced it:
#   * theta_cov_from_pint     -> t0 = MJD * 86400   ("MJD-seconds", EPOCH_MJD_SECONDS)
#   * theta_cov_from_parfile  -> t0 = GPS seconds   (cwinpy TASC/T0, EPOCH_GPS)
# Any observation window handed to window_metric/window_mismatch MUST be in the
# same convention as theta[4]; use window_from_mjd() to convert unambiguously.
# ---------------------------------------------------------------------------
EPOCH_MJD_SECONDS = "mjd_s"
EPOCH_GPS = "gps"

# O4 a+b observation window (MJD), the baseline the LVK narrowband search used.
O4AB_MJD = (60088.0, 60698.0)

J0437_FIT_COVARIANCE = (
    Path(__file__).with_name("resources") /
    "j0437_reardon2024_covariance.json"
)


def covariance_sha256(path=J0437_FIT_COVARIANCE):
    """Identify the exact covariance used by a search, not just its filename."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_covariance_match(result, path=J0437_FIT_COVARIANCE):
    """Prevent extending historical search products with a changed covariance."""
    if result.get("covariance_sha256") != covariance_sha256(path):
        raise ValueError("search result covariance hash is missing or differs; "
                         "rerun the search and injections with the current covariance")


def load_orbital_covariance(path=J0437_FIT_COVARIANCE):
    """Load and validate a canonical ``(PB, A1, ECC, OM, T0)`` covariance.

    The returned units are ``(s, light-s, 1, rad, s)``.  Artifacts consumed by
    the real cwinpy path must use the GPS epoch convention for T0 differentials.
    Validation is deliberately strict: a mislabeled or rounded near-singular
    covariance is more dangerous here than a loud failure.
    """
    path = Path(path)
    with path.open() as fh:
        artifact = json.load(fh)
    expected_order = ["PB", "A1", "ECC", "OM", "T0"]
    expected_units = ["s", "light-s", "1", "rad", "s"]
    if artifact.get("parameter_order") != expected_order:
        raise ValueError(f"{path}: expected parameter_order={expected_order}")
    if artifact.get("canonical_units") != expected_units:
        raise ValueError(f"{path}: expected canonical_units={expected_units}")
    if artifact.get("epoch_convention") != EPOCH_GPS:
        raise ValueError(f"{path}: covariance T0 must use epoch_convention='gps'")
    cov = np.asarray(artifact.get("covariance"), dtype=float)
    if cov.shape != (5, 5) or not np.all(np.isfinite(cov)):
        raise ValueError(f"{path}: covariance must be a finite 5x5 matrix")
    # Do not scale the absolute tolerance by the largest (T0) variance: doing
    # so would make the tiny ECC cross-terms effectively exempt from symmetry.
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=1e-30):
        raise ValueError(f"{path}: covariance is not symmetric")
    if np.any(np.diag(cov) <= 0):
        raise ValueError(f"{path}: covariance has non-positive variance")
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{path}: covariance is not positive definite") from exc
    return cov, artifact


def window_from_mjd(mjd_start, mjd_end, convention=EPOCH_MJD_SECONDS):
    """Convert an MJD observation window to the time units of theta[4]."""
    if convention == EPOCH_MJD_SECONDS:
        return float(mjd_start) * 86400.0, float(mjd_end) * 86400.0
    if convention == EPOCH_GPS:
        from .ephemeris import mjd_to_gps
        return mjd_to_gps(mjd_start), mjd_to_gps(mjd_end)
    raise ValueError(f"unknown epoch convention {convention!r}")


def _regularise(cov):
    cov = np.asarray(cov, float).copy()
    for i in range(cov.shape[0]):
        if cov[i, i] <= _VAR_FLOOR:
            cov[i, i] = _VAR_FLOOR
    return cov


_GOLDEN = 0.6180339887498949


def _window_times(t_start, t_end, n):
    """Deterministic low-discrepancy (golden-ratio) time samples over the window.

    A uniform grid over a long window can alias against short orbital periods
    (n=2048 over 610 d is one sample per ~7 h, while PB can be ~2 h); the golden
    sequence is equidistributed and never resonates with a periodic signal.
    """
    u = np.sort(((np.arange(int(n)) + 0.5) * _GOLDEN) % 1.0)
    return float(t_start) + (float(t_end) - float(t_start)) * u


def window_metric(theta, cov, f_gw, t_start, t_end, n=2048, eps=1e-3):
    """Whitened reduced phase metric g_z over the ABSOLUTE window [t_start, t_end].

    [t_start, t_end] must be in the same time convention as theta[4] (see the
    EPOCH_* constants / window_from_mjd). This is the physically-correct object
    for search planning: parameter errors (esp. sigma_PB) accumulate orbital
    phase with time since the fit epoch, so a window years after the epoch has a
    much larger mismatch than an epoch-centred one of the same duration.
    Worst-case mismatch from a 1-sigma offset = max eigenvalue(g_z).
    """
    theta = np.asarray(theta, float)
    L = np.linalg.cholesky(_regularise(cov))
    t = _window_times(t_start, t_end, n)
    dim = L.shape[1]
    P = np.empty((dim, t.size))
    # NOTE: roemer_delay (not bsb_delay = -roemer_delay) is deliberate and
    # harmless here. The metric is a quadratic form in the phase DERIVATIVES
    # (g ~ <d_i phi d_j phi> - <d_i phi><d_j phi>), so flipping the sign of
    # phi leaves it identically unchanged: node placement and every template
    # count are sign-invariant. The sign only matters where a signal is
    # modelled or demodulated (gravfold.orbit.bsb_delay) — don't "unify" this.
    for j in range(dim):
        d = L[:, j] * eps                       # step = eps * (1-sigma direction)
        phip = roemer_delay(t, *(theta + d))
        phim = roemer_delay(t, *(theta - d))
        P[j] = 2 * np.pi * f_gw * (phip - phim) / (2 * eps)
    # project out {1, t} (absorbed by GW phase and frequency)
    tc = t - t.mean()
    B = np.vstack([np.ones_like(t), tc / (np.std(tc) + 1e-30)]).T
    Q, _ = np.linalg.qr(B)
    P = P - (P @ Q) @ Q.T
    return (P @ P.T) / t.size


def window_mismatch(theta, cov, f_gw, t_start, t_end, n=2048):
    """Worst-case 1-sigma reduced mismatch over the absolute window."""
    g = window_metric(theta, cov, f_gw, t_start, t_end, n=n)
    return float(np.max(np.linalg.eigvalsh(g)))


def searchable_over_window(theta, cov, f_gw, t_start, t_end, mu_max=0.2, n=2048):
    """(worst-case mismatch over the window, mismatch < mu_max)."""
    mm = window_mismatch(theta, cov, f_gw, t_start, t_end, n=n)
    return mm, bool(mm < mu_max)


def whitened_reduced_metric(theta, cov, f_gw, tobs_s, n=2048, eps=1e-3):
    """DEPRECATED for ranking: epoch-centred metric (window centred on theta[4]).

    Kept for backward compatibility; it understates the mismatch for any
    observation window that does not straddle the fit epoch. Use window_metric.
    """
    theta = np.asarray(theta, float)
    return window_metric(theta, cov, f_gw,
                         theta[4] - 0.5 * tobs_s, theta[4] + 0.5 * tobs_s,
                         n=n, eps=eps)


def worst_case_mismatch(theta, cov, f_gw, tobs_s, n=2048):
    g = whitened_reduced_metric(theta, cov, f_gw, tobs_s, n=n)
    return float(np.max(np.linalg.eigvalsh(g)))


def tobs_max(theta, cov, f_gw, mu_max=0.2, tobs_lo_d=0.25, tobs_hi_d=3650.0,
             n_grid=44):
    """DEPRECATED for ranking: longest EPOCH-CENTRED coherent span (days) before a
    1-sigma orbital offset gives reduced mismatch > mu_max. Understates the
    mismatch over a real observation window displaced from the fit epoch --
    use window_mismatch/searchable_over_window with the actual window instead."""
    grid = np.geomspace(tobs_lo_d, tobs_hi_d, n_grid) * 86400.0
    mm = np.array([worst_case_mismatch(theta, cov, f_gw, T) for T in grid])
    if mm[0] >= mu_max:
        return float(tobs_lo_d), float(mm[-1])
    if mm[-1] < mu_max:
        return float(tobs_hi_d), float(mm[-1])
    k = int(np.searchsorted(mm, mu_max))
    x0, x1 = np.log(grid[k - 1]), np.log(grid[k])
    xc = x0 + (mu_max - mm[k - 1]) * (x1 - x0) / (mm[k] - mm[k - 1])
    return float(np.exp(xc) / 86400.0), float(mm[-1])


def theta_cov_from_pint(fitmodel, cvm):
    """(theta, full_cov, f_gw) in canonical (pb,x,e,om,t0) units from a fitted
    PINT model + its parameter_covariance_matrix. Uses the FULL fit covariance
    (correlations + red noise), which is mandatory for near-circular binaries
    where a diagonal treatment fabricates decorrelation. Handles ELL1 (EPS1,
    EPS2, TASC) and BT/DD/T2 (ECC, OM, T0); fixed params contribute zero variance.

    EPOCH CONVENTION: t0 = MJD * 86400 ("MJD-seconds", EPOCH_MJD_SECONDS).
    Windows for window_metric must use window_from_mjd(..., EPOCH_MJD_SECONDS).
    """
    import numpy as _np
    names = list(cvm.get_label_names()[0])
    C = _np.asarray(cvm.matrix)
    m = fitmodel
    f_gw = 2.0 * float(m.F0.value)
    pb = float(m.PB.value) * 86400.0
    x = float(m.A1.value)
    def unc(p):
        par = getattr(m, p, None)
        v = getattr(par, "uncertainty_value", None) if par is not None else None
        return float(v) if v else 0.0

    eps1_par = getattr(m, "EPS1", None)
    is_ell1 = eps1_par is not None and eps1_par.value is not None
    if is_ell1:
        native = ["PB", "A1", "EPS1", "EPS2", "TASC"]
        eps1, eps2, tasc = float(m.EPS1.value), float(m.EPS2.value), float(m.TASC.value)
        e = float(_np.hypot(eps1, eps2)); e2 = e * e if e > 0 else 1e-30
        om = float(_np.arctan2(eps1, eps2)) if e > 0 else 0.0
        t0 = tasc * 86400.0 + om / (2 * _np.pi) * pb
        sig = _np.array([unc("PB") * 86400, unc("A1"), unc("EPS1"),
                         unc("EPS2"), unc("TASC") * 86400])
        J = _np.zeros((5, 5)); J[0, 0] = 1; J[1, 1] = 1
        J[2, 2] = eps1 / e if e > 0 else 0; J[2, 3] = eps2 / e if e > 0 else 0
        J[3, 2] = eps2 / e2; J[3, 3] = -eps1 / e2
        J[4, 0] = om / (2 * _np.pi)
        J[4, 2] = pb / (2 * _np.pi) * (eps2 / e2)
        J[4, 3] = pb / (2 * _np.pi) * (-eps1 / e2)
        J[4, 4] = 1.0
    else:                                            # BT / DD / T2 (ECC,OM,T0)
        native = ["PB", "A1", "ECC", "OM", "T0"]
        e = float(m.ECC.value); om = float(_np.deg2rad(m.OM.value))
        t0 = float(m.T0.value) * 86400.0
        sig = _np.array([unc("PB") * 86400, unc("A1"), unc("ECC"),
                         _np.deg2rad(unc("OM")), unc("T0") * 86400])
        J = _np.eye(5)

    corr = _np.eye(5)
    for a in range(5):
        for b in range(5):
            if native[a] in names and native[b] in names:
                ia, ib = names.index(native[a]), names.index(native[b])
                den = _np.sqrt(C[ia, ia] * C[ib, ib])
                if den > 0:
                    corr[a, b] = C[ia, ib] / den
    cov_native = corr * _np.outer(sig, sig)
    cov = J @ cov_native @ J.T
    theta = _np.array([pb, x, e, om, t0])
    return theta, cov, f_gw


def theta_cov_from_parfile(parfile):
    """Return (theta, cov, f_gw) from a TEMPO2 .par, propagating the .par 1-sigma
    errors to the canonical (pb, x, e, om, t0) basis with correlations.

    For ELL1/T2 pulsars the native basis is (PB, A1, EPS1, EPS2, TASC); we build a
    diagonal covariance there and transform via the analytic Jacobian, which
    correctly couples (e, om, t0) for near-circular orbits. For DD the canonical
    basis is native. Requires cwinpy for exact epoch/units.

    EPOCH CONVENTION: t0 is in GPS seconds (EPOCH_GPS, cwinpy's TASC/T0 units).
    Windows for window_metric must use window_from_mjd(..., EPOCH_GPS).
    """
    from cwinpy.parfile import PulsarParameters
    par = PulsarParameters(parfile)
    if par["BINARY"] is None or par["A1"] is None:
        raise ValueError(f"{parfile}: isolated pulsar (no orbit)")
    f_gw = 2.0 * float(par["F0"])
    pb = float(par["PB"]); x = float(par["A1"])

    errs = {}
    for line in open(parfile):
        p = line.split()
        if len(p) >= 4 and p[2] in ("0", "1"):
            try:
                errs[p[0].upper()] = float(p[3].replace("D", "E"))
            except ValueError:
                pass

    if par["EPS1"] is not None:                         # ELL1 / T2
        eps1, eps2 = float(par["EPS1"]), float(par["EPS2"])
        tasc = float(par["TASC"])
        e = float(np.hypot(eps1, eps2))
        om = float(np.arctan2(eps1, eps2)) if e > 0 else 0.0
        t0 = tasc + om / (2 * np.pi) * pb
        theta = np.array([pb, x, e, om, t0])
        # measured basis m = (pb, x, eps1, eps2, tasc); diagonal 1-sigma
        sm = np.array([errs.get("PB", 0) * 86400, errs.get("A1", 0),
                       errs.get("EPS1", 0), errs.get("EPS2", 0),
                       errs.get("TASC", 0) * 86400])
        Cm = np.diag(sm ** 2)
        e2 = e * e if e > 0 else 1e-30
        # Jacobian J = d(pb,x,e,om,t0)/d(pb,x,eps1,eps2,tasc)
        J = np.zeros((5, 5))
        J[0, 0] = 1.0; J[1, 1] = 1.0
        J[2, 2] = eps1 / e if e > 0 else 0.0            # de/deps1
        J[2, 3] = eps2 / e if e > 0 else 0.0            # de/deps2
        dom_de1 = eps2 / e2; dom_de2 = -eps1 / e2       # dom/deps1, dom/deps2
        J[3, 2] = dom_de1; J[3, 3] = dom_de2
        J[4, 0] = om / (2 * np.pi)                      # dt0/dpb
        J[4, 2] = pb / (2 * np.pi) * dom_de1            # dt0/deps1 (via om)
        J[4, 3] = pb / (2 * np.pi) * dom_de2            # dt0/deps2 (via om)
        J[4, 4] = 1.0                                   # dt0/dtasc
        cov = J @ Cm @ J.T
    else:                                               # DD / BT
        e = float(par["ECC"] or 0.0); om = float(par["OM"] or 0.0)
        t0 = float(par["T0"]); theta = np.array([pb, x, e, om, t0])
        s = np.array([errs.get("PB", 0) * 86400, errs.get("A1", 0),
                      errs.get("ECC", 0), np.deg2rad(errs.get("OM", 0)),
                      errs.get("T0", 0) * 86400])
        cov = np.diag(s ** 2)
    return theta, cov, f_gw
