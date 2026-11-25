"""Cover the phase ellipsoid with adaptive slabs and one/two templates each.

For an interval in the strongest coordinate, split the second coordinate by
sign when needed. The maximum squared phase error lies among four endpoint
cases: residual norm in the second or third principal direction. Convexity
proves these cases bound the full nine-dimensional region, not just samples.
"""

import json

import numpy as np
from paths import LOFAR
from scipy.optimize import minimize

WD = LOFAR / "lofar_fit"
OUT = WD.parent / "bank_efficiency"
base = np.load(OUT / "uniform_bank.npz")
lam = base["eigenvalues"]
R = np.sqrt(float(base["critical_q"]))
reserve = 1e-4
limit = (np.sqrt(0.03) - reserve) ** 2
grid = np.linspace(-R, R, 257)
cache = {}


def solve(i, j, n):
    key = (i, j, n)
    if key in cache:
        return cache[key]
    lo, hi = grid[i], grid[j]
    x = np.sqrt(lam[0]) * np.array([lo, hi])
    remaining = np.maximum(R * R - np.array([lo, hi]) ** 2, 0)
    if n == 1:
        if j == i:
            c = x[0]
        else:
            c = np.clip((1 - lam[1] / lam[0]) * (x[0] + x[1]) / 2, x[0], x[1])
        point = np.array([c, 0.0])
        bound = float(np.max((x - c) ** 2 + lam[1] * remaining))
    else:
        points = np.array(
            [
                [x[0], 0],
                [x[0], np.sqrt(lam[1] * remaining[0])],
                [x[1], 0],
                [x[1], np.sqrt(lam[1] * remaining[1])],
            ]
        )
        residual = np.array([lam[2] * remaining[0], 0, lam[2] * remaining[1], 0])
        y2 = np.sqrt(lam[1] * max(R * R - ((lo + hi) / 2) ** 2, 0))
        y3 = np.sqrt(lam[2] / lam[1]) * y2
        p = np.array([x.mean(), (y2 * y2 - y3 * y3) / (2 * y2) if y2 else 0.0])
        initial = np.r_[p, np.max(np.sum((points - p) ** 2, axis=1) + residual)]

        def cons(v):
            return v[2] - np.sum((points - v[:2]) ** 2, axis=1) - residual

        def jac(v):
            return np.column_stack([-2 * (v[:2] - points), np.ones(4)])

        opt = minimize(
            lambda v: v[2],
            initial,
            jac=lambda v: np.array([0, 0, 1.0]),
            constraints={"type": "ineq", "fun": cons, "jac": jac},
            bounds=[(x[0], x[1]), (0, None), (0, None)],
            method="SLSQP",
            options={"ftol": 1e-12, "maxiter": 100},
        )
        point = opt.x[:2]
        bound = float(np.max(np.sum((points - point) ** 2, axis=1) + residual))
    cache[key] = (bound, point)
    return cache[key]


cost = np.full(len(grid), 100000, dtype=int)
cost[0] = 0
previous = {}
for i in range(len(grid) - 1):
    for n in [1, 2]:
        if solve(i, i, n)[0] > limit:
            continue
        low = i
        high = len(grid) - 1
        while low < high:
            mid = (low + high + 1) // 2
            if solve(i, mid, n)[0] <= limit:
                low = mid
            else:
                high = mid - 1
        for j in range(i + 1, low + 1):
            if cost[i] + n < cost[j]:
                cost[j] = cost[i] + n
                previous[j] = (i, n)
    if i % 64 == 0:
        print("slab grid", i, "of", len(grid), flush=True)
assert len(grid) - 1 in previous
slabs = []
j = len(grid) - 1
while j:
    i, n = previous[j]
    bound, point = solve(i, j, n)
    slabs.append((i, j, n, bound, point))
    j = i
slabs.reverse()
coordinates = []
edges = [-R]
assignments = []
bounds = []
for i, j, n, bound, point in slabs:
    p = np.zeros(9)
    p[:2] = point / np.sqrt(lam[:2])
    plus = len(coordinates)
    coordinates.append(p.copy())
    minus = plus
    if n == 2:
        p[1] *= -1
        minus = len(coordinates)
        coordinates.append(p.copy())
    edges.append(grid[j])
    assignments.append([minus, plus])
    bounds.append(bound)
coordinates = np.array(coordinates)
offsets = base["mean_shift"] + coordinates @ base["physical_axes"].T
worst = max(bounds)
np.savez(
    OUT / "uniform_slabs_bank.npz",
    binary_offsets=offsets,
    coordinates=coordinates,
    physical_axes=base["physical_axes"],
    critical_q=base["critical_q"],
    mean_shift=base["mean_shift"],
    edges=edges,
    slab_assignment=assignments,
    slab_bounds=bounds,
    eigenvalues=lam,
)
rows = [
    r.split("\t")[1]
    for r in (WD / "covariance.tsv").read_text().splitlines()
    if r.startswith("parameter\t")
]
np.savetxt(
    OUT / "uniform_slabs_bank_offsets.csv",
    offsets,
    delimiter=",",
    header=",".join(rows),
)
result = {
    "templates": len(coordinates),
    "slabs": len(slabs),
    "slab_template_counts": [x[2] for x in slabs],
    "linear_exact_power_loss_upper_bound": worst,
    "with_reserved_nonlinearity_loss_upper_bound": float(
        (np.sqrt(worst) + reserve) ** 2
    ),
    "reserved_nonlinear_phase_rad": reserve,
    "uniform_times": 131072,
    "critical_q": float(base["critical_q"]),
    "placement_grid_points": len(grid),
    "scope": "Complete nine-parameter older-TOA region, uniform time weights. Analytic convex endpoint certificate bounds phase variance and hence exact coherent loss in the linearized map. Native nonlinear checks separate. Minimum within the tested interval partition family, not a globally minimal waveform bank.",
}
assert result["with_reserved_nonlinearity_loss_upper_bound"] < 0.03
(OUT / "uniform_slabs_bank_results.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
