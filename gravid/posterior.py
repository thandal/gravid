"""Orbital-parameter posterior and the node sets used to marginalise over it.

This is the genuinely new piece of the contribution: instead of evaluating the
detection statistic at a single central ephemeris (and excluding the target when
the worst-case orbital error spoils a full-T_obs coherent demodulation), we place
nodes across the orbital posterior and combine.

Three node strategies are provided, cheapest first:
  * unscented   : 2n+1 sigma points (n=5 -> 11 evals). Captures mean+covariance
                  exactly for a Gaussian; the default fast estimator.
  * gauss_hermite: tensor-product Gauss-Hermite, k^n nodes. Higher accuracy,
                  curse-of-dimensionality (k=3 -> 243 nodes).
  * mc          : plain Monte-Carlo draws. The validated-final-number path.
"""
from __future__ import annotations
import numpy as np

from .orbit import PARAM_NAMES

_VAR_FLOOR = 1e-30  # so a perfectly-known parameter (sigma=0) is still invertible


class OrbitalPosterior:
    def __init__(self, mean, cov, names=PARAM_NAMES):
        self.mean = np.asarray(mean, dtype=float)
        self.cov = self._regularise(np.asarray(cov, dtype=float))
        self.names = tuple(names)
        self.dim = self.mean.size
        self.chol = np.linalg.cholesky(self.cov)

    @classmethod
    def from_ephemeris(cls, eph):
        return cls(eph.params, eph.covariance(), PARAM_NAMES)

    @staticmethod
    def _regularise(cov):
        cov = cov.copy()
        for i in range(cov.shape[0]):
            if cov[i, i] <= _VAR_FLOOR:
                cov[i, i] = _VAR_FLOOR
        return cov

    # ---- node sets -------------------------------------------------------
    def sample_mc(self, n, rng=None):
        rng = np.random.default_rng() if rng is None else rng
        z = rng.standard_normal((int(n), self.dim))
        pts = self.mean + z @ self.chol.T
        w = np.full(int(n), 1.0 / int(n))
        return pts, w

    def nodes_unscented(self, alpha=1.0, beta=2.0, kappa=0.0):
        n = self.dim
        lam = alpha ** 2 * (n + kappa) - n
        scale = np.linalg.cholesky((n + lam) * self.cov)
        pts = np.empty((2 * n + 1, n))
        wm = np.empty(2 * n + 1)
        pts[0] = self.mean
        wm[0] = lam / (n + lam)
        for i in range(n):
            pts[1 + i] = self.mean + scale[:, i]
            pts[1 + n + i] = self.mean - scale[:, i]
            wm[1 + i] = wm[1 + n + i] = 1.0 / (2.0 * (n + lam))
        return pts, wm  # weights for mean (expectation) integration

    def nodes_gauss_hermite(self, k=3):
        # Probabilists' Hermite -> nodes/weights for the N(0,1) density.
        x, w = np.polynomial.hermite_e.hermegauss(k)
        w = w / np.sqrt(2.0 * np.pi)
        grids = np.meshgrid(*([x] * self.dim), indexing="ij")
        wgrids = np.meshgrid(*([w] * self.dim), indexing="ij")
        z = np.stack([g.ravel() for g in grids], axis=1)
        weights = np.prod(np.stack([g.ravel() for g in wgrids], axis=1), axis=1)
        pts = self.mean + z @ self.chol.T
        return pts, weights

    def nodes(self, method="unscented", **kw):
        if method == "mc":
            return self.sample_mc(kw.get("n", 1000), kw.get("rng"))
        if method == "gauss_hermite":
            return self.nodes_gauss_hermite(kw.get("k", 3))
        if method == "unscented":
            return self.nodes_unscented()
        raise ValueError(f"unknown node method {method!r}")


def metric_node_count(g_z, mu_max=0.2, z_range=3.0):
    """How many templates `nodes_metric` would place under whitened metric
    `g_z` — WITHOUT building them (the grid can be millions of nodes).

    This is the SAME counting rule nodes_metric uses (per eigendirection,
    n_i = ceil(2*z_range / sqrt(mu/lam_i)), collapsing to 1 where the spacing
    already covers the box), factored out so that anything reasoning about
    cost — the adaptive-segmentation planner, the template-count tables —
    counts what the search will actually place. A separate
    re-derivation of this arithmetic is how a planner's cap silently stops
    meaning what it says.
    """
    g = np.asarray(g_z, dtype=float)
    g = 0.5 * (g + g.T)
    lam = np.clip(np.linalg.eigvalsh(g), 0.0, None)
    n = 1.0
    for lam_i in lam:
        if lam_i <= 0.0:
            continue
        dz = np.sqrt(float(mu_max) / lam_i)
        if dz < 2.0 * z_range:
            n *= float(np.ceil(2.0 * z_range / dz))
    return n


def nodes_metric(posterior: OrbitalPosterior, g_z, mu_max=0.2, z_range=3.0,
                 max_nodes=20000):
    """Deterministic metric-spaced posterior nodes (no randomness).

    Works in the whitened space z (theta = mean + L z, prior N(0, I)). Given the
    whitened phase metric g_z over the SEARCH WINDOW (window_metric with this
    posterior's covariance), eigendecompose g_z = V diag(lam) V^T; in the rotated
    coordinates y = V^T z the mismatch is sum_i lam_i y_i^2 and the prior is
    still N(0, I). Along eigendirection i the node spacing is
    dz_i = sqrt(mu_max / lam_i) (directions with negligible lam collapse to a
    single node), tensor-gridded over [-z_range, +z_range]; the weight of each
    node is the Gaussian prior mass of its cell (product of 1-D cell integrals,
    i.e. a truncated-at-z_range prior). If the total node count would exceed
    max_nodes, the mismatch budget is raised uniformly (reported as
    info['effective_mu']).

    Coverage: the worst-case mismatch from any point in the grid box to its
    nearest node is ~ (n_active/4) * effective_mu.

    Returns (pts [n,dim], weights [n], info dict with n_nodes, effective_mu,
    n_per_dim, eigvals, eigvecs).
    """
    from scipy.special import ndtr
    g = np.asarray(g_z, dtype=float)
    g = 0.5 * (g + g.T)
    lam, V = np.linalg.eigh(g)
    # Canonicalize eigenvector signs (largest-|component| positive). LAPACK
    # makes NO promise about the sign of an eigenvector, and a flipped column
    # would mirror the node grid along that axis — same node set, permuted
    # order — which would silently break the nodes[shard::nshards] stride
    # assignment across machines. Defence in depth only: in the measured
    # cross-CPU behaviour (an OPENBLAS_CORETYPE sweep
    # over ZEN/HASWELL/SANDYBRIDGE/NEHALEM/PRESCOTT on the production image)
    # the signs never actually flipped; what varies across BLAS kernel
    # families is the LAST BITS of the metric/eigenvectors (~1e-15 relative),
    # which moves node coordinates infinitesimally but never reorders them.
    # That byte-level jitter is handled at merge time by a numeric grid
    # comparison (gpu_search.merge_shards), not here.
    for i in range(V.shape[1]):
        k = int(np.argmax(np.abs(V[:, i])))
        if V[k, i] < 0.0:
            V[:, i] = -V[:, i]
    lam = np.clip(lam, 0.0, None)
    dim = posterior.dim

    def counts_for(mu):
        n_i = np.ones(dim, dtype=np.int64)
        for i in range(dim):
            if lam[i] <= 0.0:
                continue
            dz = np.sqrt(mu / lam[i])
            if dz < 2.0 * z_range:
                n_i[i] = int(np.ceil(2.0 * z_range / dz))
        return n_i

    mu = float(mu_max)
    n_i = counts_for(mu)
    # inflate the mismatch budget until the tensor grid fits under max_nodes
    for _ in range(200):
        total = float(np.prod(n_i.astype(float)))
        if total <= max_nodes:
            break
        n_active = max(int(np.sum(n_i > 1)), 1)
        mu *= (total / max_nodes) ** (2.0 / n_active) * 1.05
        n_i = counts_for(mu)

    axes_c, axes_w = [], []
    for i in range(dim):
        edges = np.linspace(-z_range, z_range, int(n_i[i]) + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        w = ndtr(edges[1:]) - ndtr(edges[:-1])   # Gaussian mass per cell
        axes_c.append(centers)
        axes_w.append(w)
    mesh_c = np.meshgrid(*axes_c, indexing="ij")
    mesh_w = np.meshgrid(*axes_w, indexing="ij")
    y = np.stack([m.ravel() for m in mesh_c], axis=1)
    w = np.prod(np.stack([m.ravel() for m in mesh_w], axis=1), axis=1)
    z = y @ V.T                                   # z = V y (rows)
    pts = posterior.mean + z @ posterior.chol.T
    info = dict(n_nodes=int(y.shape[0]), effective_mu=float(mu),
                n_per_dim=n_i.tolist(), eigvals=lam, eigvecs=V)
    return pts, w, info
