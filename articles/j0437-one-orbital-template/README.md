# One orbital template for a coherent continuous-wave search of PSR J0437−4715

`article.md` is the text, `article.pdf` the typeset copy and
`references.bib` its references. All paths below are from the repository root.

## Section 3: search of LIGO O4a data

| In the article | Computed by | Saved in |
|---|---|---|
| Sample counts, dropped fractions, noise weights | `gravid/noise_weights.py` | `results/j0437_search/result.json`, key `noise_weights` |
| Table 2, column 2𝒢 | `scripts/j0437_upper_limit.py` | `results/j0437_search/result.json` |
| Table 2, column 2ℱ | `scripts/j0437_f_statistic.py` | `results/j0437_search/f_statistic.json` |
| Table 2 as printed, and the power captured by each filter | `scripts/search/noise_weights/summarize.py` | `scripts/search/noise_weights/results/summary.json` |
| Scan of 7.2 million frequency bins | `combine_scan.py` in the same folder | `scripts/search/noise_weights/results/frequency_scan.json` |
| Figure 1 | `plot_results.py` in the same folder | `figures/j0437_search_results.png` here |

`scripts/search/noise_weights/README.md` explains the noise
weights, lists the checks and gives the commands.

Rebuilding Table 2 and Figure 1 from the saved outputs takes seconds:

```bash
cd scripts/search/noise_weights
python summarize.py
python plot_results.py
```

Rerunning the searches needs the reduced strain (102 files, 393 MB).
TODO: where the reduced strain is deposited, and how it was made from the GWOSC release.

## Sections 1 and 2: timing

Scripts and saved files are under `scripts/timing/` unless a full path is given.

| In the article | Computed by | Saved in |
|---|---|---|
| Refit with TEMPO2, covariance of the ten binary parameters, search model | `reproduce_j0437_full_binary_covariance.py`, `build_j0437_search_model.py` | `gravid/resources/j0437_de440_binary_covariance.json`, `gravid/resources/j0437_de440_search.par`, `binary_model/covariance.tsv` |
| Correlation 0.99999936 between the period and periastron advance | `validate_j0437_full_model.py` | `binary_model/model_validation.json`, key `PB_OMDOT_correlation` |
| Eccentricity, signal frequency, the two error relations (within 0.002%), the 0.30-ms tolerance | `later_radio_observations/summarize.py` | `later_radio_observations/results/summary.json`, keys `eccentricity`, `frequency_Hz`, `degeneracy_ratios`, `sinusoidal_tolerances` |
| Table 1, the numbers of observations, the residual scatter | same script | same file, key `predictions` |
| A separate dispersion correction each day preserves coherence | `supporting_checks/audit_external_prediction.py` | `supporting_checks/external_prediction_audit.json` |
| Orbital offset from 36 later MeerKAT observations, −66 ± 46 ns | `orbital_offset_fit/epoch_orbit_fit.py` | `orbital_offset_fit/epoch_orbit_fit.json`, key `meerkat/post` |
| Its robust standard error, 55 ns, and the InPTA result | `orbital_offset_fit/review_epoch_fit.py` | `orbital_offset_fit/orbital_fit_review.json` |
| Pulsar Portal O4b check: 543 arrival times, 17 observations, loss 0.0014% | `portal_arrival_times/validate.py`, `summarize.py` | `portal_arrival_times/results.json`, key `J0437-4715/adopted/O4b` |
| Largest loss 3.73e-8 over 256 draws and 1,024 corners; agreement of PINT, TEMPO2 and LALSuite delays | `validate_j0437_uncertainty.py` | `binary_model/uncertainty_crosscheck.json`, key `summary` |
| Independent QR solve | `export_qr_covariance.C` | `qr_solve/independent_qr_check.json` |
| Fits to the two halves of the timing data | `supporting_checks/analyze_j0437_halves.py` | `supporting_checks/j0437_free_halves.json` |
| Withheld year: 24% more scatter than predicted | `withheld_year/predictive_linear.py` | `withheld_year/heldout_linear.json`, key `chi2_per_toa` (1.53; its square root is 1.24) and `binary_derivative_subspace` |

The saved files let every number above be read directly. Rerunning the
scripts needs TEMPO2 and the radio arrival times from the cited releases;
their outputs go under `results/timing/`, which is not in the repository.
TODO: download steps for the radio data, and where the count of 11,637 arrival times is saved.
