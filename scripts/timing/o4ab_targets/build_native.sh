#!/usr/bin/env bash
set -eu
src_dir=$(cd "$(dirname "$0")" && pwd)
repo_dir=$(cd "$src_dir/../../.." && pwd)
prefix=${O4AB_TEMPO2_PREFIX:-"$repo_dir/.conda"}
out=${O4AB_WORK:-"$repo_dir/results/timing/o4ab_targets"}/bin
mkdir -p "$out"
flags=(-O2 "-I$prefix/include" "-L$prefix/lib" "-Wl,-rpath,$prefix/lib")
g++ "${flags[@]}" -shared -fPIC "$src_dir/native/export.C" -ltempo2 -o "$out/commoncov_linux-gnu_plug.t2"
g++ "${flags[@]}" "$src_dir/native/delays.C" -ltempo2 -o "$out/delays"
g++ "${flags[@]}" "$src_dir/native/predict_design.C" -ltempo2 -o "$out/predict_design"
g++ "${flags[@]}" "$src_dir/../orbital_offset_fit/evaluate_timing.C" -ltempo2 -o "$out/evaluate"
