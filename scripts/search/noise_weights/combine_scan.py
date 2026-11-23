#!/usr/bin/env python3
"""Fixed-orientation statistic at every frequency bin within +-0.3 Hz.

Reads the per-detector transforms saved by
``scripts/j0437_f_statistic.py --scan`` and checks that the
control band is flat and chi-squared distributed.  Under 1.5 GB.
"""
import json

import numpy as np

from _common import BULK, HERE
from gravfold.upper_limit import heterodyne_statistic_response

H, L = np.load(BULK / "scan_H1.npz"), np.load(BULK / "scan_L1.npz")
assert int(H["n_fft"]) == int(L["n_fft"]) and float(H["start"]) == float(L["start"])
n, kmax = int(H["n_fft"]), int(H["kmax"])
norm = dict(H1=float(H["norm"]), L1=float(L["norm"]))
f = (np.arange(2 * kmax + 1) - kmax) / n
yh, yl = H["band"].astype(np.complex128), L["band"].astype(np.complex128)
stat = dict(coherent=2 * np.abs(yh + yl) ** 2 / (norm["H1"] + norm["L1"]),
            H1=2 * np.abs(yh) ** 2 / norm["H1"], L1=2 * np.abs(yl) ** 2 / norm["L1"])
del yh, yl
response = heterodyne_statistic_response(f)
band = (np.abs(f) >= 0.03) & (np.abs(f) <= 0.2)
out = dict(bins_in_control_band=int(band.sum()), n_fft=n,
           expected_largest=float(2 * np.log(band.sum())), detectors={}, bands=[])
for name, values in stat.items():
    s = values[band] / response[band]
    out["detectors"][name] = dict(
        mean=float(s.mean()), fraction_above_1pct_point=float(np.mean(s > 9.2103)),
        largest=float(s.max()), target=float(values[kmax]))
edges = np.round(np.arange(-0.20, 0.2001, 0.01), 2)
for a, b in zip(edges[:-1], edges[1:]):
    sel = (f >= a) & (f < b) & (np.abs(f) > 1e-4)
    out["bands"].append([float(a), float(b)] + [
        float(np.mean(stat[k][sel] / response[sel])) for k in ("coherent", "H1", "L1")])
means = np.array([row[2:] for row in out["bands"]])
out["band_mean_range"] = dict(zip(("coherent", "H1", "L1"),
                                  [[float(lo), float(hi)] for lo, hi in
                                   zip(means.min(axis=0), means.max(axis=0))]))
print(json.dumps({k: v for k, v in out.items() if k != "bands"}, indent=1))
(HERE / "results/frequency_scan.json").write_text(json.dumps(out, indent=1) + "\n")
