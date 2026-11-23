#!/usr/bin/env python3
"""Measure the noise autocorrelation of the reduced strain and compare the
implied equivalent noise bandwidth with the one predicted from the filter
response.  For 1-s samples the bandwidth in Hz is 1 / (1 + 2 sum of the
autocorrelation over positive lags)."""
import json

import numpy as np
from scipy.ndimage import median_filter

from _common import HERE, LN2, files, load_chunk, segments
from gravfold.noise_weights import equivalent_noise_bandwidth

out = dict(predicted_hz=equivalent_noise_bandwidth())
for det in ("H1", "L1"):
    num, den = np.zeros(41, complex), np.zeros(41)
    for path in files(det)[3::8]:
        _, times, data, _ = load_chunk(path)
        power = np.abs(data) ** 2
        for a, b in zip(*segments(times)):
            spans = [(a, a + 3480), (a + 3540, b - 60)] if b - a > 3540 else [(a, b - 60)]
            for s, e in spans:
                if e - s < 400:
                    continue
                z = data[s:e] / np.sqrt(median_filter(power[s:e], size=301, mode="reflect") / LN2)
                ok = np.abs(z) ** 2 <= 20
                ok &= ~(np.convolve(~ok, np.ones(5), "same") > 0)
                z = np.where(ok, z, 0)
                for lag in range(41):
                    num[lag] += np.sum(z[:z.size - lag] * np.conj(z[lag:]))
                    den[lag] += np.sum(ok[:z.size - lag] & ok[lag:])
    rho = num / den
    rho /= rho[0].real
    out[det] = dict(lag1=float(rho[1].real), largest_imaginary=float(np.abs(rho.imag).max()),
                    measured_hz=float(1 / (1 + 2 * np.sum(rho[1:].real))))
print(json.dumps(out, indent=1))
(HERE / "results/noise_bandwidth.json").write_text(json.dumps(out, indent=1) + "\n")
