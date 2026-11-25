"""Exact linear GLS/REML with the native full design and epoch noise.
Whitening uses the exported TOA errors. Low-rank identities avoid dense nTOA matrices.
All epoch-process amplitudes are in microseconds. No posterior interpretation.
"""

import numpy as np
from paths import BASE
from scipy.linalg import cho_factor, cho_solve, qr, solve_triangular
from scipy.optimize import minimize


def load(name):
    wd = BASE / (name + "_ELL1")
    n, p, nc = map(int, (wd / "design.meta").read_text().splitlines()[0].split())
    assert nc == 0
    data = np.fromfile(wd / "toas.bin").reshape(n, 4)
    M = np.fromfile(wd / "design.bin").reshape(n, p)
    meta = np.loadtxt(wd / "design.meta", skiprows=1, dtype=int)
    return wd, data, M, meta


def kernel(t, u, theta, kind, span, diagonal=False):
    dt = np.abs(t[:, None] - u[None, :])
    K = np.zeros_like(dt)
    if kind != "white" and diagonal:
        K += np.eye(len(t)) * np.exp(2 * theta[1])
    if kind == "red":
        amp = np.exp(theta[2])
        gamma = theta[3]
        freq = np.arange(1, 31) / span
        weights = np.arange(1, 31, dtype=float) ** (-gamma)
        weights /= weights.sum()
        K += amp**2 * np.sum(
            np.cos(2 * np.pi * dt[:, :, None] * freq) * weights, axis=2
        )
    if kind == "exp":
        K += np.exp(2 * theta[2]) * np.exp(-dt / np.exp(theta[3]))
    return K


def epoch_ids(times):
    # Separate observations on the same day; channel TOAs are nearly simultaneous.
    order = np.argsort(times)
    groups = np.cumsum(np.r_[True, np.diff(times[order]) * 86400 > 300]) - 1
    inverse = np.empty(len(times), dtype=int)
    inverse[order] = groups
    return inverse


class Fit:
    def __init__(self, data, M, span):
        self.data = data
        self.span = span
        self.spectrum_key = None
        self.inverse = epoch_ids(data[:, 0])
        self.days = np.unique(self.inverse)
        self.t = np.array(
            [data[self.inverse == i, 0].mean() for i in range(len(self.days))]
        )
        self.err = data[:, 2]
        self.y = data[:, 1] / self.err
        norm = np.linalg.norm(M / self.err[:, None], axis=0)
        valid = np.flatnonzero(norm > 0)
        A = M[:, valid] / self.err[:, None] / norm[valid]
        _, r, piv = qr(A, mode="economic", pivoting=True)
        rank = np.sum(abs(np.diag(r)) > abs(r[0, 0]) * 1e-11)
        if rank != len(valid):
            raise ValueError(("rank deficient", rank, len(valid)))
        self.cols = valid
        self.scale = norm[valid]
        self.A = A
        Q, R = np.linalg.qr(A, mode="reduced")
        Ri = solve_triangular(R, np.eye(len(valid)))
        self.Q = Q
        self.Ri = Ri
        self.C0 = Ri @ Ri.T
        self.beta0 = Ri @ Q.T @ self.y
        residual = self.y - Q @ (Q.T @ self.y)
        self.q = float(residual @ residual)
        self.df = len(data) - rank
        W = np.zeros((len(data), len(self.days)))
        W[np.arange(len(data)), self.inverse] = 1e-6 / self.err
        self.W = W
        self.WTW = W.T @ W
        self.D = Ri @ (Q.T @ W)
        PW = W - Q @ (Q.T @ W)
        self.RW = PW.T @ PW
        self.b = W.T @ residual

    def set_y(self, y):
        self.y = y
        self.beta0 = self.Ri @ self.Q.T @ y
        r = y - self.Q @ (self.Q.T @ y)
        self.q = float(r @ r)
        self.b = self.W.T @ r

    def evaluate(self, theta, kind, details=False):
        s2 = np.exp(2 * theta[0])
        if kind == "red":
            if self.spectrum_key != (self.span, theta[3]):
                frequencies = np.arange(1, 31) / self.span
                if getattr(self, "cached_span", None) != self.span:
                    dt = np.abs(self.t[:, None] - self.t[None, :])
                    self.cosines = np.cos(2 * np.pi * dt[:, :, None] * frequencies)
                    self.cached_span = self.span
                weights = np.arange(1, 31, dtype=float) ** (-theta[3])
                weights /= weights.sum()
                self.spectrum = np.sum(self.cosines * weights, axis=2)
                self.spectrum_key = (self.span, theta[3])
            K = (
                np.eye(len(self.t)) * np.exp(2 * theta[1])
                + np.exp(theta[2]) ** 2 * self.spectrum
            )
        else:
            K = kernel(self.t, self.t, theta, kind, self.span, True)
        if kind == "white":
            Z = K
            logdet = 0.0
        else:
            L = np.linalg.cholesky(K)
            cf = cho_factor(np.eye(len(K)) + L.T @ self.RW @ L / s2, lower=True)
            Z = L @ cho_solve(cf, L.T)
            logdet = 2 * np.log(np.diag(cf[0])).sum()
        obj = self.df * np.log(s2) + logdet + (self.q - self.b @ Z @ self.b / s2) / s2
        if not details:
            return float(obj)
        beta = self.beta0 - self.D @ Z @ self.b / s2
        cov = s2 * self.C0 + self.D @ Z @ self.D.T
        return float(obj), beta, cov, K

    def optimize(self, kind, start=None, quick=False):
        if kind == "white":
            theta = np.array([0.5 * np.log(self.q / self.df)])
        else:
            bounds = [(np.log(0.3), np.log(4)), (np.log(1e-4), np.log(100))]
            starts = [[0, np.log(0.3)], [0, np.log(2)]]
            if kind == "red":
                bounds += [(np.log(1e-4), np.log(300)), (0, 7)]
                starts = [
                    [0, np.log(0.3), np.log(1), 2],
                    [0, np.log(1), np.log(10), 4],
                    [0, np.log(0.1), np.log(0.1), 1],
                ]
            if kind == "exp":
                bounds += [(np.log(1e-4), np.log(300)), (np.log(5), np.log(5000))]
                starts = [
                    [0, np.log(0.3), np.log(1), np.log(100)],
                    [0, np.log(1), np.log(10), np.log(500)],
                ]
            if start is not None:
                starts = [start] if quick else [start] + starts
            fits = [
                minimize(
                    self.evaluate,
                    x,
                    args=(kind,),
                    method="L-BFGS-B",
                    bounds=bounds,
                    options={"ftol": 1e-10, "maxiter": 300},
                )
                for x in starts
            ]
            best = min(fits, key=lambda x: x.fun)
            if not best.success and np.linalg.norm(best.jac) > 1e-2:
                raise RuntimeError(str(best))
            theta = best.x
        obj, beta, cov, K = self.evaluate(theta, kind, True)
        return {
            "theta": theta,
            "objective": obj,
            "beta": beta,
            "cov": cov,
            "K": K,
            "EFAC": float(np.exp(theta[0])),
        }

    def predictive(self, fit, kind, test, Mtest, E):
        theta = fit["theta"]
        s2 = np.exp(2 * theta[0])
        K = fit["K"]
        if kind == "white":
            Z = K
        else:
            L = np.linalg.cholesky(K)
            cf = cho_factor(np.eye(len(K)) + L.T @ self.WTW @ L / s2, lower=True)
            Z = L @ cho_solve(cf, L.T)

        def cinv(B):
            return (B - self.W @ (Z @ (self.W.T @ B)) / s2) / s2

        te = (
            E @ test[:, 0]
        )  # E sums to one but signed chromatic extrapolation; times nearly simultaneous
        Mt = (E @ Mtest[:, self.cols]) / self.scale * 1e6
        Ks = kernel(te, self.t, theta, kind, self.span)
        Kss = kernel(te, te, theta, kind, self.span, True)
        latent = Ks @ self.W.T
        pred = Mt @ fit["beta"] + latent @ cinv(self.y - self.A @ fit["beta"])
        G = Mt - latent @ cinv(self.A)
        white = s2 * (E * test[:, 2]) @ (E * test[:, 2]).T * 1e12
        V = white + Kss - latent @ cinv(latent.T) + G @ fit["cov"] @ G.T
        return E @ test[:, 1] * 1e6 - pred, (V + V.T) / 2, te


def epoch_projection(data):
    ids = epoch_ids(data[:, 0])
    days = np.unique(ids)
    E = []
    for day in days:
        idx = np.flatnonzero(ids == day)
        a = data[idx]
        X = np.column_stack([np.ones(len(a)), (1000 / a[:, 3]) ** 2])
        w = a[:, 2] ** -2
        row = np.zeros(len(data))
        row[idx] = np.linalg.solve(X.T @ (w[:, None] * X), X.T * w)[0]
        E.append(row)
    return np.array(E)


def serialize(fit):
    return {
        k: v.tolist() if isinstance(v, np.ndarray) else v
        for k, v in fit.items()
        if k not in ["K", "beta", "cov"]
    }
