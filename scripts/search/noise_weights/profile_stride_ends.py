#!/usr/bin/env python3
"""Sample power against position within each segment of reduced strain.

Shows the transient before each stride end: position 3540 (3600 s after the
segment start, less the 60-s crop) and the end of the segment.  Power is
normalised by each segment's median-based variance; 1 is normal noise.
"""
import json

import numpy as np

from _common import HERE, LN2, load, segments

out = {}
for det in ("H1", "L1"):
    times, data = load(det)
    power = np.abs(data) ** 2
    del data
    lo, hi = segments(times)
    n = hi - lo
    full = lo[n == 3976]
    grid = np.stack([power[a:a + 3976] / (np.median(power[a:a + 3400]) / LN2)
                     for a in full])
    median = np.median(grid, axis=0) / LN2
    tail = {f"{k}_to_{k - 5}_s_before_end": dict(
        first_stride_end=float(median[3540 - k:3540 - k + 5].mean()),
        segment_end=float(median[3976 - k:3976 - k + 5].mean()))
        for k in range(60, 0, -5)}
    partial = [(a, b) for a, b in zip(lo, hi) if 400 <= b - a < 3976]
    ends = np.stack([power[b - 60:b] / (np.median(power[a:b - 150]) / LN2)
                     for a, b in partial])
    out[det] = dict(
        samples=int(times.size), segments=int(lo.size), full_segments=int(full.size),
        plain_over_median_variance=float(power.mean() / (np.median(power) / LN2)),
        baseline_median_500_to_3400=float(median[500:3400].mean()),
        after_first_stride_end_3540_to_3600=float(median[3540:3600].mean()),
        median_power_by_5s_before_stride_end=tail,
        shorter_segments=len(partial),
        shorter_segment_last_5s=float((np.median(ends, axis=0) / LN2)[-5:].mean()),
        shorter_segment_60_to_30s_before_end=float((np.median(ends, axis=0) / LN2)[:30].mean()))
    print(det, json.dumps(out[det], indent=1))
(HERE / "results/stride_end_profile.json").write_text(json.dumps(out, indent=1) + "\n")
