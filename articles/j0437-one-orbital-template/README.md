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

Rerunning the searches needs the reduced strain; see Data below.

### What the two saved search outputs are

`result.json` and `f_statistic.json` in `results/j0437_search/` are
byte-for-byte what the two searches wrote on 20 November 2026. They were
written before this code was published, under the project's working name
"gravfold" and into a folder named `results/j0437_search_rerun_20261120/`;
both names appear inside the files. The code in this repository is the code
that ran, with names changed and comments edited. The searches have not been
rerun from this repository.

`result.json` carries a record of the run under `run_capsule`. What can be
checked from here:

```python
import hashlib, json
from pathlib import Path
from gravid.run_capsule import capsule_digest, chunk_manifest_digest

sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
result = json.loads(Path("results/j0437_search/result.json").read_text())
record = result["run_capsule"]
assert capsule_digest(record) == record["capsule_sha256"]
assert chunk_manifest_digest(record["input_hashes"]) == record["chunk_manifest_sha256"]
assert sha("gravid/resources/j0437_de440_search.par") == result["search_par_sha256"]
assert sha("gravid/resources/j0437_de440_binary_covariance.json") == result["covariance_sha256"]
```

`record["input_hashes"]` holds the SHA-256 of each of the 102 reduced strain
files. The record's digests of the source tree, the commit and the
uncommitted changes refer to the unpublished working repository and cannot
be checked from here.

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
scripts needs TEMPO2 and the radio arrival times; see Data below. Their
outputs go under `results/timing/`, which is not in the repository.

### Regenerating Table 1 from the releases

```bash
python scripts/timing/orbital_offset_fit/acquire.py        # MeerKAT and InPTA files, clock files
# save the PPTA collection's parameter file as data/J0437-4715_reardon2024_dr3.par
python scripts/timing/later_radio_observations/summarize.py --native-rerun
```

The last step builds the 14 timing models (seven per release), evaluates
their residuals with TEMPO2 without fitting anything, and prints the table.
It reproduces the saved values to nine significant figures. `acquire.py`
rewrites `download_manifest.json` with the checksums of what it fetched, so
`git diff` shows whether any upstream file has changed.

## Data

None of the input data is in this repository.

| Data | Source | What fixes the exact files |
|---|---|---|
| PPTA arrival times for J0437 (11,637), parameter file and clock files | CSIRO data collection, doi:10.25919/20rx-5f63 | `scripts/reproduce_j0437_covariance.py` holds the SHA-256 of the parameter and arrival-time files (`PAR_SHA256`, `TIM_SHA256`) and refuses any others |
| MeerKAT 4.5-year release | `https://docs.datacentral.org.au/documents/52/partim.tar.gz` | `scripts/timing/orbital_offset_fit/acquire.py` downloads it; `download_manifest.json` beside it lists the SHA-256 of every downloaded file |
| InPTA DR2 | `github.com/inpta/InPTA.DR2` and `github.com/inpta/InPTA.DR2.NA`, at the revisions named in `acquire.py` | same script and list |
| Observatory clock corrections | `github.com/ipta/pulsar-clock-corrections` | same list; the upstream files change over time |
| Pulsar Portal arrival times for the O4b check | `https://pulsars.org.au/download/J0437-4715/toas`, downloaded 16 October 2026 | 1,765,503 bytes, SHA-256 `7f632f47d40691a4067670c60693092d8f36b2777e773b01437033d340e1b13c` |
| LIGO O4a strain | GWOSC O4a release at 4 kHz, channels `H1:GWOSC-4KHZ_R1_STRAIN` and `L1:GWOSC-4KHZ_R1_STRAIN` | see below |

### Reduced strain

The searches read 102 files of reduced strain (412 MB, 26,366,991 one-second
samples), one per detector per five days, from `results/campaign_het_v2/`.
They were made from the GWOSC O4a release, GPS 1368195220 to 1389456018:

1. Times with data were taken from the GWOSC timeline (flags `H1_DATA` and
   `L1_DATA`).
2. Strain for those times was downloaded at 4 kHz with gwpy's
   `TimeSeries.fetch_open_data`.
3. CWInPy 1.5.0 heterodyned it in pieces of 432,000 s with these settings:
   stride 3,600 s; ninth-order low-pass filter at 0.25 Hz; output at 1 Hz;
   60 s cropped from both ends of every segment; Solar-system and binary
   barycentring; twice the rotation frequency; interpolation step 60 s. The
   parameter file for this step, `base_cwinpy.par`, uses DE440 and the ELL1
   binary model. The searches then apply the phase difference to the T2
   model.

All 102 files were made with the same settings. Each file stores them,
together with the parameter file and its list of segments. The SHA-256 of
every file is in `results/j0437_search/result.json` under `input_identity`,
`input_hashes`.

The files are in one archive, `j0437_o4a_reduced_strain.tar` (412 MB), in
this public folder:

https://drive.google.com/drive/folders/1GZsjyQHLniwpRwuOSjeJ_90WFVWtmRWF

Its SHA-256 is `935b116a91a0832a41e39b91b51e4ffb2c547753fc0544d2b7834568449ee2f5`.
Unpack it at the repository root with `tar -xf j0437_o4a_reduced_strain.tar`;
the files land in `results/campaign_het_v2/`. The archive also holds
`base_cwinpy.par`, the checksum record the search reads and a list of
per-file checksums. The files are derived from the GWOSC release and are
under the same CC BY 4.0 licence.
