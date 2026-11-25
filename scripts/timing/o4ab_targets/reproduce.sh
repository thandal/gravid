#!/usr/bin/env bash
set -euo pipefail
src=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$src/../../.." && pwd)
export O4AB_WORK=$(realpath -m "${O4AB_WORK:-$repo/results/timing/o4ab_targets}")
export O4AB_TEMPO2_PREFIX=$(realpath -m "${O4AB_TEMPO2_PREFIX:-$repo/.conda}")
timing_python=${O4AB_PYTHON:-$O4AB_TEMPO2_PREFIX/bin/python}
cw_python=${O4AB_CW_PYTHON:-$repo/.venv-fresh/bin/python}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$O4AB_WORK/noise" "$O4AB_WORK/logs"
run() {
    printf '%s\n' "$*"
    "$timing_python" -u "$src/$1.py" "${@:2}" > "$O4AB_WORK/logs/$1${2:+_$2}.log" 2>&1
}
run prepare_inputs
bash "$src/build_native.sh"
run fit_orbits meerkat
run phase_axes
for step in noise_fits frequency_checks bootstrap profile_noise physical_grid bank_envelope native_bank_check; do
    run "$step"
done
run native_bank_check --high-only
run fixed_parameter_sensitivity
run radio_predictions
run external_predictive
run fit_orbits lofar
for step in lofar_covariance bootstrap_lofar lofar_metric bank_uniform_slabs validate_efficient_bank historical; do
    run "$step"
done
"$cw_python" "$src/j0437_comparison.py" > "$O4AB_WORK/logs/j0437_comparison.log" 2>&1
"$timing_python" "$src/paper_results.py" --verify-run
