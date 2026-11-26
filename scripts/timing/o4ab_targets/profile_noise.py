"""Noise-profile sensitivity envelope, not a hyperparameter posterior.
Grid includes noise amplitudes/gamma well beyond fitted values and two low-frequency cutoffs.
"""

import json

import numpy as np
from gls import Fit, load
from paths import NAMES, OUT
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

threshold = chi2.ppf(0.99, 4)


def amplitude_profiles(fit, base_objective, jitter, gamma, amplitudes):
    """Keep accepted grid points and approach each bracketed cutoff from inside."""
    def profile(amplitude):
        theta = np.array([0.0, np.log(jitter), np.log(amplitude), gamma])
        scale = minimize_scalar(
            lambda x: fit.evaluate(np.r_[x, theta[1:]], "red"),
            bounds=(np.log(0.6), np.log(1.6)),
            method="bounded",
            options={"xatol": 1e-6},
        )
        if not scale.success:
            raise RuntimeError(scale.message)
        theta[0] = scale.x
        return theta, scale.fun - base_objective

    grid = [(amp, profile(amp)) for amp in amplitudes]
    retained = {amp: point for amp, point in grid if point[1] <= threshold}
    for (a, left), (b, right) in zip(grid, grid[1:]):
        if (left[1] <= threshold) == (right[1] <= threshold):
            continue
        inside, outside = (a, b) if left[1] <= threshold else (b, a)
        point = retained[inside]
        while abs(np.log(outside / inside)) > 1e-4:
            mid = np.sqrt(inside * outside)
            trial = profile(mid)
            if trial[1] <= threshold:
                inside, point = mid, trial
            else:
                outside = mid
        retained[inside] = point
    return retained.items()


results = {}
for name in NAMES:
    wd, data, M, meta = load(name)
    rows = [
        s.split("\t")
        for s in (wd / "covariance.tsv").read_text().splitlines()
        if s.startswith("parameter\t")
    ]
    indices = [
        next(i for i, (a, c) in enumerate(meta) if a == int(r[2]) and c == 0)
        for r in rows
    ]
    formal = np.load(wd / "robust_covariances.npz")["formal_QR"]
    sd = np.sqrt(formal.diagonal())
    L = np.linalg.cholesky(formal / np.outer(sd, sd)) * sd[:, None]
    Li = np.linalg.inv(L)
    d = np.fromfile(wd / "axes.bin").reshape(19, 8192)
    D = (d[1:10] - d[10:]) / 0.2
    D -= D.mean(axis=1)[:, None]
    f0 = float(
        next(
            s.split()[1]
            for s in (wd / "native.par").read_text().splitlines()
            if s.startswith("F0 ")
        )
    )
    G = (2 * np.pi * 2 * f0 * 1.001) ** 2 * D @ D.T / 8192
    cases = []
    for extension in [1, 2]:
        fit = Fit(data, M, np.ptp(data[:, 0]) * extension)
        base = fit.optimize("red")
        pos = [list(fit.cols).index(i) for i in indices]
        fac = np.array([float(r[4]) / float(r[5]) for r in rows]) / fit.scale[pos]
        jittervals = sorted({
            1e-4, 0.25, 0.3, *np.arange(0.5, 2.125, 0.125),
            2.5, 3.0, 5.0, 10.0, float(np.exp(base["theta"][1])),
        })
        ampvals = sorted({
            1e-4, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, float(np.exp(base["theta"][2])),
        })
        gammas = sorted({
            *np.arange(0.0, 2.25, 0.25), 3.0, 4.0, 5.0, 6.0, 7.0,
            float(base["theta"][3]),
        })
        for gamma in gammas:
            for jitter in jittervals:
                profiles = amplitude_profiles(fit, base["objective"], jitter, gamma, ampvals)
                for amp, (th, delta_reml) in profiles:
                    _, b, C, _ = fit.evaluate(th, "red", True)
                    Cb = C[np.ix_(pos, pos)] * np.outer(fac, fac)
                    mu = b[pos] * fac
                    W = Li @ Cb @ Li.T
                    v, V = np.linalg.eigh((W + W.T) / 2)
                    S = V @ np.diag(np.sqrt(np.maximum(v, 0)))
                    ev = np.linalg.eigvalsh(S.T @ G @ S)
                    shift = mu @ Li.T
                    bound = (
                        np.sqrt(max(shift @ G @ shift, 0))
                        + np.sqrt(chi2.ppf(0.95, 9) * max(ev[-1], 0))
                    ) ** 2
                    cases.append(
                        {
                            "span_factor": extension,
                            "EFAC": float(np.exp(th[0])),
                            "jitter_us": jitter,
                            "red_rms_us": amp,
                            "gamma": gamma,
                            "delta_REML": delta_reml,
                            "nominal95_linear_loss_bound": float(bound),
                        }
                    )
    results[name] = {
        "profile_threshold": threshold,
        "threshold_description": "99% chi-square(4) reference, a sensitivity cutoff, not calibrated boundary-noise confidence coverage",
        "accepted_cases": len(cases),
        "maximum_nominal95_linear_loss_bound": max(
            x["nominal95_linear_loss_bound"] for x in cases
        ),
        "worst_case": max(cases, key=lambda x: x["nominal95_linear_loss_bound"]),
        "cases": cases,
    }
    (OUT / "noise_profile_results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(
        name, results[name]["accepted_cases"], results[name]["worst_case"], flush=True
    )
