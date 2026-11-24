# Noise weights and results of the J0437 search

How the reduced O4a strain is weighted in the search of PSR J0437−4715, the
checks on those weights, and the numbers in the article's Table 2 and
Figure 1.

## Three properties of the reduced strain

1. **A transient before the end of every stride.** CWInPy 1.5.0 reads and
   filters the strain one stride (3600 s) at a time: forward, then backward
   starting from the forward pass's final filter state. That is the wrong
   state for the reversed series, so the last ~30 s of each stride carry a
   start-up transient. Strides end 3600 s after the start of each segment and
   at the end of the segment, inside the data kept after CWInPy's 60-s crop.
   Median power relative to normal noise, in 5-s steps over the last 30 s:
   1.1, 4, 37, 730, 4,600 at Hanford and 1.8, 13, 150, 2,800, 18,000 at
   Livingston (`results/stride_end_profile.json`). These samples make the
   plain variance 25 (Hanford) and 218 (Livingston) times the median-based
   variance.
2. **Noise level changes within minutes, plus short glitches.** One variance
   per file, or per half hour, does not follow it.
3. **Filtered noise is not white.** The 1-Hz samples have passed a 0.25-Hz
   low-pass filter. The noise power spectral density at the signal frequency
   is the sample variance divided by the filter's equivalent noise bandwidth,
   0.446 Hz (predicted from the filter response; 0.452 and 0.453 Hz measured
   from the autocorrelation at the two detectors,
   `results/noise_bandwidth.json`), not by the 1-Hz sample rate.

## The weights

`gravid/noise_weights.py` (tests in `tests/test_noise_weights.py`), selected
with `scripts/j0437_upper_limit.py --noise-weights running-median`:

- drop the 60 s before each stride end (3.1% of samples);
- noise variance from a 301-s running median of sample power;
- drop samples whose power exceeds 20 times that variance, with 2 s either
  side (0.25%);
- weight by the inverse noise power spectral density (variance over
  equivalent noise bandwidth), so the statistic follows a chi-squared
  distribution in Gaussian noise.

25,478,590 of 26,366,991 samples are kept. For a constant-amplitude signal,
the estimated signal-to-noise ratio is the sum of the weights divided by
the square root of the sum of squared weights times sample power. We
combine both detectors before taking the ratio. The baseline weights are
the inverse variances of the untrimmed five-day files.

| Weights and cuts | Gain in estimated signal-to-noise ratio over baseline |
|---|---:|
| Drop the last 60 s of each filter block, keep the original file weights | 3.40 |
| Apply the same cut and recompute each file variance | 3.75 |
| Apply all cuts and use the 301-s running median | 5.89 |

These values use the measured noise in the reduced strain. Table 2 in the
article measures recovery by injecting signals with a range of amplitudes
and orientations. Across median windows of 61–1201 s and power cuts of
9–50, the gain on the constant-amplitude diagnostic ranges from 5.6 to 6.0.

The two search statistics use the same orbit. The reference filter assumes
one linear polarization (cos ι = 0, ψ = 0). It captures a median 48% of
the signal power recovered by the F-statistic, which fits both polarization
responses. At the fifth percentile the reference filter captures 16%.
Across the injected orbits the largest power loss for the F-statistic is
3.7e-9.

## Results

`results/summary.json`, from `summarize.py`. Amplitudes are recomputed exactly
from the saved sums. Spin-down limit 7.92e-27; published LVK O4a limits
7.2e-27 (Bayesian), 9.0e-27 (F-statistic), 5.8e-27 (5n-vector).

| | Reference filter (2𝒢) | F-statistic (2ℱ) |
|---|---:|---:|
| Target statistic | 2.52 | 9.85 |
| Probability in Gaussian noise | 0.28 | 0.043 |
| Controls at or above target (of 40) | 11 | 2 |
| Control mean (expected) | 2.28 (2) | 4.07 (4) |
| Loudest control | 15.58 | 16.89 |
| 95% recovery above loudest control | 3.46e-26 | 1.60e-26 |
| its 90% bootstrap interval | 2.73–3.98e-26 | 1.49–1.72e-26 |
| 95% recovery at 1% false alarm | 2.55e-26 | 1.40e-26 |
| same, expected over noise realisations | 2.64e-26 | 1.42e-26 |
| 95% recovery above the observed target | 1.42e-26 | 1.20e-26 |

`scripts/j0437_upper_limit.py` itself reports 3.59e-26 for the reference
filter because it interpolates on a 16-point amplitude grid.

## Checks

- **Statistic at every frequency bin** (`results/frequency_scan.json`): over
  the 7.2 million bins of the control band the fixed-orientation statistic
  averages 2.04 (Hanford alone 2.01, Livingston alone 2.07); 1.08% exceed the
  1% point; the largest is 32.4 against about 31.6 expected. Means over
  0.01-Hz bands stay within 2.01–2.05. The scan shows no corresponding line or band feature. Under independent
  Gaussian noise the chance that at least one of 40 controls reaches the
  observed maximum of 15.58 is 1.64%.
- The two search scripts give identical fixed-orientation target, controls
  and injection sums; `summarize.py` asserts this.
- Expected recovery over noise realisations agrees with the realised one.
- Not done: comparison with LALSuite's own F-statistic code.

## Reproduce

Run from the repository root in an environment with the `cwinpy` extra
installed.

```bash
OUT=results/j0437_search
SCAN=results/timing/noise_weights
mkdir -p "$SCAN"
# fixed orientation: under 6.5 GB, about 30 minutes
python scripts/j0437_upper_limit.py --stage upperlimit --nproc 1 --n-offsets 40 \
    --n-inj 200 --noise-weights running-median --out $OUT/result.json
# F-statistic: about 2.7 GB and 13 minutes per detector, then a light combine
python scripts/j0437_f_statistic.py --detector H1 --out-dir $OUT --scan $SCAN/scan_H1.npz
python scripts/j0437_f_statistic.py --detector L1 --out-dir $OUT --scan $SCAN/scan_L1.npz
python scripts/j0437_f_statistic.py --combine --out-dir $OUT
# diagnostics and summary (each under 1.8 GB)
cd scripts/search/noise_weights
python profile_stride_ends.py && python check_noise_bandwidth.py
python compare_weight_choices.py && python combine_scan.py && python summarize.py
python plot_results.py   # writes the article's Figure 1
```
