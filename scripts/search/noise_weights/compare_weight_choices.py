#!/usr/bin/env python3
"""How the running-median length, the power cut and its margin change sensitivity.

Compares (sum w)^2 / sum(w^2 |d|^2), the signal-to-noise power a
constant-amplitude signal would receive using the realised noise, with the
same quantity for one plain variance per five-day file.  About 6
minutes.
"""
import json

import numpy as np
from scipy.ndimage import median_filter

from _common import HERE, LN2, files, load_chunk, segments
from gravid.noise_weights import stride_end_mask

table, totals = {}, {}
for det in ("H1", "L1"):
    acc, base = {}, [0.0, 0.0]
    trimmed = {"original file variance": np.zeros(2), "trimmed file variance": np.zeros(2)}
    for path in files(det):
        _, times, data, _ = load_chunk(path)
        power = np.abs(data) ** 2
        v = float(np.var(data))
        base[0] += power.size / v
        base[1] += float((power / v ** 2).sum())
        lo, hi = segments(times)
        keep = stride_end_mask(times)
        kept_power = power[keep]
        if kept_power.size:
            for name, variance in (("original file variance", v),
                                   ("trimmed file variance", np.var(data[keep]))):
                weight = 1 / variance
                trimmed[name] += (kept_power.size * weight, kept_power.sum() * weight**2)
        for window in (61, 121, 301, 601, 1201):
            variance = np.full(power.size, np.nan)
            for a, b in zip(lo, hi):
                index = a + np.flatnonzero(keep[a:b])
                if index.size >= 300:
                    variance[index] = median_filter(power[index], size=window,
                                                    mode="reflect") / LN2
            ok = np.isfinite(variance)
            ratio = np.where(ok, power / np.where(ok, variance, 1), 0)
            for cut in (np.inf, 50, 20, 12, 9):
                for margin in (0, 2):
                    use = ok & (ratio <= cut)
                    if margin and np.isfinite(cut):
                        use &= ~(np.convolve(ok & (ratio > cut),
                                             np.ones(2 * margin + 1), "same") > 0)
                    w = np.where(use, 1 / np.where(ok, variance, 1), 0)
                    e = acc.setdefault((window, cut, margin), [0.0, 0.0, 0])
                    e[0] += float(w.sum()); e[1] += float((w * w * power).sum())
                    e[2] += int(use.sum())
    totals[det] = (base, acc, trimmed)
    reference = base[0] ** 2 / base[1]
    for (window, cut, margin), e in sorted(acc.items()):
        table.setdefault(f"median {window} s, power cut {cut}, margin {margin} s", {})[det] = dict(
            kept=e[2], gain_over_file_variance=e[0] ** 2 / e[1] / reference,
            noise_power_per_unit_weight=e[1] / e[0])
reference = (sum(totals[d][0][0] for d in totals) ** 2
             / sum(totals[d][0][1] for d in totals))
for key in sorted(totals["H1"][1]):
    w = sum(totals[d][1][key][0] for d in totals)
    n = sum(totals[d][1][key][1] for d in totals)
    name = f"median {key[0]} s, power cut {key[1]}, margin {key[2]} s"
    table[name]["both_gain_over_file_variance"] = w * w / n / reference
    print(f"{name:44s} gain in signal-to-noise power {w*w/n/reference:6.2f}  "
          f"(amplitude {np.sqrt(w*w/n/reference):.2f})")
for name in totals["H1"][2]:
    sum_weights, noise = totals["H1"][2][name] + totals["L1"][2][name]
    gain = sum_weights**2 / noise / reference
    table[f"trim block ends only, {name}"] = {"both_gain_over_file_variance": float(gain)}
(HERE / "results/weight_choices.json").write_text(json.dumps(table, indent=1) + "\n")
