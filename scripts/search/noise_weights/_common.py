"""Shared helpers for the search noise-weight diagnostics (one detector at a time)."""
import glob
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).absolute().parents[3]
sys.path.insert(0, str(ROOT))
from gravfold.j0437_campaign import load_chunk  # noqa: E402

HERE = Path(__file__).absolute().parent
BULK = ROOT / "results/timing/noise_weights"
LN2 = float(np.log(2))


def files(det):
    return sorted(glob.glob(str(ROOT / f"results/campaign_het_v2/het_*_{det}_*.hdf5")))


def load(det):
    """Times and complex samples of one detector, all files (about 0.4 GB)."""
    parts = [load_chunk(f)[1:3] for f in files(det)]
    return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])


def segments(times):
    breaks = np.flatnonzero(np.diff(times) > 1.5) + 1
    return np.r_[0, breaks], np.r_[breaks, times.size]


def cplx(a):
    a = np.asarray(a, float)
    return a[..., 0] + 1j * a[..., 1]


def statistic_coefficients(f_run):
    """Coefficients (a, b, c) of each statistic as a function of amplitude.

    With injection j added at amplitude h, the statistic is exactly
    a_j h^2 + b_j h + c, where c is the statistic of the data alone.  Returns
    the saved sums too: data vector, matrix and per-injection signal vectors.
    """
    products = f_run["products"]
    yd, m, ys = cplx(products["y_data"]), cplx(products["m"]), cplx(products["y_signal"])
    mi = np.linalg.inv(m)
    f_statistic = (2 * np.real(np.einsum("ji,ik,jk->j", np.conj(ys), mi, ys)),
                   4 * np.real(np.einsum("i,ik,jk->j", np.conj(yd), mi, ys)),
                   2 * float(np.real(np.conj(yd) @ mi @ yd)))
    fixed_orientation = (2 * np.abs(ys[:, 0]) ** 2 / m[0, 0].real,
                         4 * np.real(np.conj(yd[0]) * ys[:, 0]) / m[0, 0].real,
                         2 * abs(yd[0]) ** 2 / m[0, 0].real)
    return dict(fixed_orientation=fixed_orientation, f_statistic=f_statistic,
                y_signal=ys)
