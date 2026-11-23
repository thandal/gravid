#!/usr/bin/env python3
"""Export all ten fitted binary parameters from native TEMPO2 QR.

The five-parameter marginal block cannot be paired with a waveform that fixes
OMDOT: its near-degeneracy with PB must be propagated jointly.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.reproduce_j0437_covariance import freeze_unmatched_jumps


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--par", type=Path, default=ROOT / "results/timing/ppta_dr3_puzzle/csiro/dr3/J0437-4715.par")
    ap.add_argument("--tim", type=Path, default=ROOT / "results/timing/ppta_dr3_puzzle/csiro/dr3/J0437-4715.tim")
    ap.add_argument("--clocks", type=Path, default=ROOT / "results/timing/ppta_dr3_puzzle/csiro/clock")
    ap.add_argument("--work-dir", type=Path, default=ROOT / "results/j0437_full_timing_covariance")
    ap.add_argument("--out", type=Path, default=ROOT / "gravfold/resources/j0437_de440_binary_covariance.json")
    ap.add_argument("--export-exists", action="store_true")
    ap.add_argument("--verify-snapshot", type=Path)
    args = ap.parse_args()
    previous = json.loads((ROOT / "gravfold/resources/j0437_reardon2024_covariance.json").read_text())
    for path, key in ((args.par, "source_par_sha256"), (args.tim, "source_tim_sha256")):
        if sha(path) != previous["provenance"][key]:
            raise ValueError(f"input hash mismatch: {path}")
    for name, digest in previous["provenance"]["clock_sha256"].items():
        if sha(args.clocks / name) != digest:
            raise ValueError(f"clock hash mismatch: {name}")
    wd = args.work_dir.resolve()
    wd.mkdir(parents=True, exist_ok=True)
    source = ROOT / "scripts/timing/export_full_binary_covariance.C"
    envroot = ROOT / ".conda"
    if not args.export_exists:
        subprocess.run(["g++", "-std=c++17", "-O2", "-fPIC", "-shared",
                        f"-I{envroot/'include'}", str(source), f"-L{envroot/'lib'}",
                        f"-Wl,-rpath,{envroot/'lib'}", "-ltempo2", "-o",
                        str(wd / "fullcov_linux-gnu_plug.t2")], check=True)
        freeze_unmatched_jumps(args.par, wd / "input.par")
        env = dict(os.environ, TEMPO2=str(envroot / "share/tempo2"),
                   TEMPO2_PLUG_PATH=str(wd), AUDIT_COV_OUTPUT=str(wd / "covariance.tsv"),
                   OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2")
        with (wd / "run.log").open("w") as log:
            subprocess.run([str(envroot / "bin/tempo2"), "-clkdir", str(args.clocks.resolve()),
                            "-f", str(wd / "input.par"), str(args.tim.resolve()),
                            "-qrfit", "-norescale", "-output", "fullcov", "-tcheck"],
                           env=env, cwd=wd, stdout=log, stderr=subprocess.STDOUT, check=True)
    lines = (wd / "covariance.tsv").read_text().splitlines()
    rows = [line.split("\t") for line in lines if line.startswith("parameter\t")]
    raw = np.loadtxt(lines[lines.index("covariance_internal") + 1:])
    if not np.array_equal(raw[:5, :5], previous["covariance_tcb_internal"]):
        raise ValueError("expanded export does not reproduce adopted five-parameter block")
    if sha(wd / "native.par") != previous["provenance"]["fit_par_sha256"]:
        raise ValueError("expanded export changed the native fitted mean")
    scale = previous["tcb_to_tdb_linear_scale"]
    factors = np.array([86400*scale, scale, 1, 1, 86400*scale,
                        1e12, 1/scale, 1, 1, scale])
    cov = raw * np.outer(factors, factors)
    np.linalg.cholesky(cov)
    result = dict(
        schema_version=1, target="J0437-4715", parameter_order=[r[1] for r in rows],
        canonical_units=["s", "light-s", "1", "rad", "s", "1e-12", "deg/yr", "rad", "rad", "solar_mass"],
        epoch_convention="gps", covariance=cov.tolist(), covariance_tcb_internal=raw.tolist(),
        covariance_kind="native TEMPO2 QR fixed-hyperparameter GLS covariance",
        normalization="unscaled (-norescale)",
        fit_values_tcb_par_units={r[1]: r[3] for r in rows},
        provenance={**previous["provenance"], "output_plugin_sha256": sha(source),
                    "fit_export_sha256": sha(wd / "covariance.tsv"),
                    "notes": "Full fitted binary block, including PB-OMDOT covariance; other timing/noise coefficients marginalized."},
    )
    if args.verify_snapshot:
        from scipy.linalg import lapack, solve_triangular
        prefix = args.verify_snapshot / "solver"
        n, p, nc = map(int, prefix.with_suffix(".meta").read_text().splitlines()[0].split())
        a = np.empty((n+nc, p), order="F")
        a[:n] = np.fromfile(str(prefix)+".white_design.bin").reshape(n, p)
        a[n:] = np.fromfile(str(prefix)+".constraints.bin").reshape(nc, p)
        norms = np.linalg.norm(a, axis=0)
        a /= norms
        work, status = lapack.dgels_lwork(n+nc, p, 1)
        qr, _, status = lapack.dgels(a, np.zeros(n+nc), lwork=int(work), overwrite_a=True)
        if status:
            raise RuntimeError(f"independent QR failed: {status}")
        indices = np.array([int(r[2]) for r in rows])
        z = solve_triangular(np.triu(qr[:p]).T, np.eye(p)[:, indices], lower=True)
        independent = (z.T @ z) / np.outer(norms[indices], norms[indices])
        corr = lambda c: c / np.sqrt(np.outer(c.diagonal(), c.diagonal()))
        check = dict(max_sigma_fractional_difference=float(np.max(abs(np.sqrt(independent.diagonal()/raw.diagonal())-1))),
                     max_correlation_difference=float(np.max(abs(corr(independent)-corr(raw)))),
                     scope="Independent column-scaled QR of exact native design/prior snapshot",
                     design_sha256=sha(str(prefix)+".white_design.bin"),
                     constraints_sha256=sha(str(prefix)+".constraints.bin"))
        if check["max_sigma_fractional_difference"] > 1e-7 or check["max_correlation_difference"] > 1e-7:
            raise ValueError(check)
        result["independent_verification"] = check
        print(json.dumps(check, indent=2))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
