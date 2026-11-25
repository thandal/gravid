# Reassessing binary pulsars for the O4ab narrowband search

[Read the article](article.pdf) · [Manuscript](article.md) · [References](references.bib) · [Analysis code](../../scripts/timing/o4ab_targets/)

The three MeerKAT results use two independent sets of observations. The earlier
TOAs determine how far the orbit can move within the timing errors. The O4b
TOAs test the predicted pulse phases at later dates. Both columns of Table 1
are fractions of coherent power lost. No detector strain enters this calculation.

## Reproduce Table 1 from this repository

Run from the repository root on Linux with a C++ compiler, Python, NumPy and
TEMPO2. The native programs use the TEMPO2 library and headers in `.conda/`
by default. Set `O4AB_TEMPO2_PREFIX` to use another installation of the same
version. This run used conda-forge TEMPO2 2024.12.1 (build `hddb8a8a_4`).
The optional `O4AB_WORK` sets the output directory for both commands.

```bash
bash scripts/timing/o4ab_targets/build_native.sh
OPENBLAS_NUM_THREADS=1 .conda/bin/python scripts/timing/o4ab_targets/paper_results.py
```

This calculation reads the fitted reference and limiting orbits, evaluates
fresh TEMPO2 delays at 131,072 times and finds the largest allowed phase
excursion. It also averages the saved complex phase measurements from the
O4b epochs. The reference entry lists the mean UTC MJD and channel count
for each epoch in the order used by all predictions for that pulsar. For J1231 it takes the range over 25 predictions. It compares
the recomputed values with the saved table and exits with an error on a
discrepancy. The three unrounded orbital bounds, in percent to four decimal
places, are **0.1401, 0.9733 and 1.7128**. The table rounds these upwards.
[Orbital regions](../../scripts/timing/o4ab_targets/results/limiting_orbits.json) ·
[O4b epoch measurements](../../scripts/timing/o4ab_targets/results/epoch_predictions.json) ·
[Table and supporting results](../../scripts/timing/o4ab_targets/results/paper_results.json).

## Refit from the arrival times

On 25 November 2026 the full driver rebuilt the fits in an empty working
directory from 161 verified inputs. The comparison passed for the reported
bounds, O4b predictions, LOFAR template count and historical checks. An
independent run of the Table 1 calculation used only files in this repository.

The full calculation starts with TOAs, clock corrections and the initial
timing models in `scripts/timing/o4ab_targets/models/`. Bulk timing inputs
are stored outside Git. With the archived `results/timing` tree alongside
this checkout, the following commands copy the inputs, check their SHA-256
hashes, run all fits and compare the results with the publication files.
The starter models include the fitted daily dispersion terms and the
J1737 conversion to ELL1. TEMPO2 refits the parameters flagged in each
starter model. It exports a fresh design matrix and covariance.

```bash
.conda/bin/python scripts/timing/o4ab_targets/prepare_inputs.py \
    --archive ../gravfold/results/timing \
    --partim ../gravfold/results/timing/j0437_external_validation/meerkat/partim.tar.gz
bash scripts/timing/o4ab_targets/reproduce.sh
```

`--archive` accepts any copy of the archived tree. `--partim` optionally
rebuilds the three training sets directly from the original public MeerKAT
tarball. A copy of the prepared `inputs/` directory also works. In that
case run `reproduce.sh` directly. Its first step verifies all input hashes.
The default output is `results/timing/o4ab_targets/`. Step logs go under
`logs/`. The full calculation uses Python 3.11, NumPy 1.26.4, SciPy 1.17.1
and PINT-pulsar 1.1.5 in the TEMPO2 environment. The J0437 comparison uses
the companion article's environment with CWInPy 1.5.0 and LALSuite 8.0.2.
Set `O4AB_PYTHON` and `O4AB_CW_PYTHON` to the two interpreter paths if
they differ from the project defaults. The shell driver fixes the number of
BLAS and OpenMP threads at one. [Input filenames and hashes](../../scripts/timing/o4ab_targets/input_manifest.json).

### Sources of the timing inputs

| Input | Selection |
|---|---|
| MeerKAT training TOAs | The three `.tim` files in the 4.5-year release give 3,231 J1045, 3,258 J1737 and 2,129 J1231 TOAs. The last training observation is at MJD 60157 to the nearest day. [Release tarball](https://docs.datacentral.org.au/documents/52/partim.tar.gz) · [Miles et al.](https://doi.org/10.1093/mnras/stae2572) |
| MeerKAT O4b TOAs | Retain MJD 60410.625 through 60703.70833. J1045 and J1737 use the public PTA stream with 32 channels, total intensity, one integration and no DM correction. The PTA2 stream duplicates these rows. J1231 uses the RelBin timing download. Preserve the decimal MJD strings. [Pulsar Portal](https://pulsars.org.au/) · [Data-use policy](https://pulsars.org.au/data-usage/) |
| J1400 Nançay and LOFAR | The Nançay model is the J1400 member of the LVK ephemeris bundle. The Donner release uses J1400−1438. The input is its eight-channel TOA file with author-rejected outliers removed. The appended station flags allow the fit to use one offset per station relative to the first. [LVK release](https://dcc.ligo.org/LIGO-T2600052/public) · [LOFAR release](https://doi.org/10.5281/zenodo.4290012) |
| J1400 LWA | Three reduced 49.8-MHz TOAs at MJDs 58628, 60111 and 60592, plus dispersion measured from the profiles. This replay starts at those measurements. The native residuals assume GPS-referenced station times. [LWA archive](https://lda10g.alliance.unm.edu/PulsarArchive/) · [Kumar et al.](https://doi.org/10.3847/1538-4357/adb97a) |
| J1745 | The orbit and five errors come from Table 1 of Janssen et al. The historical prediction uses the Portal folding model and 256 PTA TOAs at eight epochs. [Janssen et al.](https://doi.org/10.1051/0004-6361/200911728) · [Pulsar Portal](https://pulsars.org.au/) |
| Clocks | The snapshot includes the MeerKAT clock chain, GPS-to-UTC and TT(BIPM2024). The low-frequency station tables apply the archived GPS corrections. Constant offsets for the LOFAR stations are included in the fit. [IPTA clock tables](https://github.com/ipta/pulsar-clock-corrections) |

You can also rebuild the eight timing files for MeerKAT and LOFAR from the
original releases. Place `partim.tar.gz`, `dvt+20.tar`, the Portal timing
downloads `J1231-1411_toas.zip` and `J1745-0952_toas.zip`, and the two
public observation exports `J1045-4509_api_toas.json` and
`J1737-0811_api_toas.json` in one directory. The export script records
the original Portal fields and selection. The conversion checks that PTA
and PTA2 agree, then verifies the eight resulting files byte for byte.
The ninth input contains the three LWA measurements reduced from the archived
profiles. It ships with the article.
[Original Portal export script](https://github.com/thandal/gravfold/blob/39af6625296b268eba17c061bd6406ab83db5624/scripts/timing/gap_closure_20261017/fetch_public_toas.py).

```bash
.conda/bin/python scripts/timing/o4ab_targets/release_toas.py /path/to/release-products
.conda/bin/python scripts/timing/o4ab_targets/prepare_inputs.py --clocks /path/to/T2runtime/clock
bash scripts/timing/o4ab_targets/reproduce.sh
```

The second command uses the TEMPO2 runtime and the October 2026 IPTA clock
tables. It reconstructs the LOFAR and LWA station tables from the GPS series,
copies the three LWA measurements and verifies all 161 inputs. The clock
snapshot must match the hashes in the manifest. [Reduced LWA TOAs](../../scripts/timing/o4ab_targets/results/J1400-1431_lwa.tim).

The Portal snapshots date from October 2026. The input hashes identify the
exact selections used here. The public Table 1 calculation above works
without a Portal download or the author's archive.

## Where the numbers come from

All scripts below are in [o4ab_targets](../../scripts/timing/o4ab_targets/).
The result keys refer to [paper_results.json](../../scripts/timing/o4ab_targets/results/paper_results.json).
Power losses in JSON are fractions. Multiply by 100 for percent.

| Article result | Calculation | Result key |
|---|---|---|
| J0437 comparison over O4ab | `j0437_comparison.py` evaluates 256 Gaussian draws and 1,024 corners on 16,384 equally weighted times with the companion paper's orbital covariance | `J0437` |
| Joint errors and Table 1 bounds | `fit_orbits.py`, `phase_axes.py`, `bootstrap.py`, `profile_noise.py`, `physical_grid.py`, `bank_envelope.py`, `native_bank_check.py` | `bootstrap`, `noise_grid`, `grid_and_nonlinear_checks`, `table_1` |
| Additional fitted parameters | `fixed_parameter_sensitivity.py` refits parallax and the second spin-frequency derivative | `extra_parameters` |
| O4b epochs, sampled losses and 0.3-ms controls | `radio_predictions.py` evaluates the full earlier timing models at the later TOAs. Each complex number is an inverse-variance mean over the channels at one epoch | `table_1`, [epoch predictions and controls](../../scripts/timing/o4ab_targets/results/epoch_predictions.json) |
| O4b trend and orbital tests | `external_predictive.py` includes training-fit uncertainty and cross-epoch noise | `radio_prediction_statistics` |
| Five published errors for J1400 and J1745, LOFAR drift, LWA and historical J1745 predictions | `historical.py` and the three LOFAR fits in `fit_orbits.py` | `historical` |
| LOFAR bank | `lofar_covariance.py`, `bootstrap_lofar.py`, `lofar_metric.py`, `bank_uniform_slabs.py`, `validate_efficient_bank.py` | `LOFAR_placement`, `LOFAR_native_checks` |

The MeerKAT error regions use the 196th ordered joint-error statistic from
200 trials, giving squared radii of 19.485, 22.132 and 21.695. The noise
grid retains 247, 341 and 275 profiles within 13.277 of the minimum in
minus twice restricted log likelihood. The J1231 run crosses the 275
profiles with the reference fit and all 24 mass and inclination choices
(6,875 combinations). The dense-grid check uses the limiting combination
for each pulsar. `noise_fits.py` and `frequency_checks.py` also write the
quarter-span and radio-frequency split diagnostics in the working directory.
Their nominal $p$-values are saved under `timing_splits`.

The LOFAR offsets are in the native TEMPO2 units of the named columns.
Add each row to the free LOFAR timing fit. The count of 12 applies to
uniform weights on the O4ab grid and is a sufficient count. The linear
variance bound is 2.9961%. Allowing an additional 0.0001 rad of phase error
gives 2.9995%. Full TEMPO2 trials tested the nonlinear allowance.
[Reference fit](../../scripts/timing/o4ab_targets/results/orbits/J1400-1431_lofar.par) ·
[Twelve template offsets](../../scripts/timing/o4ab_targets/results/lofar_offsets.csv).

The five-parameter calculation uses period and ascending-node epoch in days
and projected semimajor axis in light-seconds. The earlier draft interchanged
the axis and epoch errors in its analytic bound. With the units corrected,
the J1400 bounds are 0.293% (three standard deviations) and 2.07% (eight).
The J1745 three-standard-deviation bound is 2.26%.

## Build the PDF

```bash
node scripts/build_article_pdf.cjs articles/o4ab-binary-targets/article.md
```

## Literature search, 25 November 2026

We checked the five pulsars by name, including J1400−1438, and followed
citations to the O4ab search and its orbital selection rule in INSPIRE and
to the O4ab paper in OpenAlex. We found no earlier paper reporting the O4ab
power-loss bounds and the O4b radio checks in Table 1.
[Works citing the O4ab search](https://inspirehep.net/literature?q=refersto%3Arecid%3A3136113) ·
[Works citing Mirasola et al.](https://inspirehep.net/literature?q=refersto%3Arecid%3A2774573) ·
[OpenAlex citations](https://openalex.org/works?filter=cites:W7167612366).

- LIGO reported S3/S4 targeted upper limits for J1045 and J1745. Its S5 search included J1045 and estimated power loss from correlated timing errors. The 16 filters estimated for J1045 in 2001 used the radio errors available then. [S3/S4 results](https://arxiv.org/pdf/gr-qc/0702039) · [S5, §3 and Table 1](https://arxiv.org/pdf/0909.3583) · [Dhurandhar & Vecchio, Table III](https://arxiv.org/pdf/gr-qc/0011085).
- The O4a targeted search used Nançay timing for J1231 and J1745. The public data products give upper limits. They do not include the J1745 arrival times needed to test our extrapolation. [LVK, §4.2 and Tables 2–3](https://arxiv.org/pdf/2501.01495) · [O4a data products](https://dcc.ligo.org/G2402400/public).
- Cromartie et al. give an independent J1231 orbit through MJD 60408. Their companion data release contains posterior samples for the two masses, inclination and parallax. Our 25 MeerKAT fits test sensitivity to fixed masses and inclinations. [Timing solution and Shapiro measurement](https://arxiv.org/html/2511.10529v2) · [Posterior samples](https://doi.org/10.5281/zenodo.17596801).
