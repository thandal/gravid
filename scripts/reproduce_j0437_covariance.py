#!/usr/bin/env python3
"""Reconstruct J0437's covariance with the main TEMPO2 QR fitter.

The output-only plugin exports the marginalized orbital block. Noise
hyperparameters stay fixed; no chi-squared or published-error rescaling is
applied. --verify-design also checks a post-fit relinearization using SciPy QR.
Bulk TOAs, clocks, and optional design matrices remain outside version control.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parents[1]
PAR_SHA256 = "98ecbfc528f9c3b65446b364a74c95594bf73ec2fb4a30badbb68f66d8bb3cb9"
TIM_SHA256 = "ee6a2dec40b4dc6fe73fc06c8ddc43052d2f0c24edc703bec8963b58c2541564"
TCB_TO_TDB = 0.9999999844948022
ORBITAL_ORDER = ("PB", "A1", "ECC", "OM", "T0")
UNMATCHED_JUMPS = {1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 56}
CLOCK_FILES = ("pks2gps.clk", "gps2utc.clk", "utc2tai.clk", "tai2tt_bipm2020.clk", "leap.sec")
REFERENCE = HERE / "gravfold/resources/j0437_reardon2024_covariance.json"
# TEMPO2 internal orbital covariance units: day, light-s, 1, radian, day.
CANONICAL_SCALE = np.array([86400*TCB_TO_TDB, TCB_TO_TDB, 1, 1, 86400*TCB_TO_TDB])


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_hash(path, expected):
    if sha256(path) != expected:
        raise ValueError(f"{path}: input SHA-256 does not match the pinned release")


def freeze_unmatched_jumps(source, destination):
    out, count = [], 0
    for line in Path(source).read_text().splitlines():
        fields = line.split()
        if fields and fields[0].upper() == "JUMP":
            count += 1
            if count in UNMATCHED_JUMPS:
                if fields[-1] != "1":
                    raise ValueError(f"JUMP {count}: expected fit flag 1")
                fields[-1] = "0"
                line = " ".join(fields)
        out.append(line)
    if count != 60:
        raise ValueError(f"expected 60 released JUMPs, found {count}")
    Path(destination).write_text("\n".join(out) + "\n")


def published_errors(par):
    rows = [line.split() for line in Path(par).read_text().splitlines()]
    return {key: float(next(r[3] for r in rows if r and r[0] == key))
            for key in ORBITAL_ORDER}


def read_native_export(path):
    lines = Path(path).read_text().splitlines()
    if lines[0] != "NATIVECOV1":
        raise ValueError("not a native TEMPO2 covariance export")
    pos = lines.index("covariance_internal")
    parameters = [line.split("\t") for line in lines[:pos] if line.startswith("parameter\t")]
    if [p[1] for p in parameters] != list(ORBITAL_ORDER):
        raise ValueError("incorrect orbital parameter order")
    cov = np.loadtxt(lines[pos+1:])
    if cov.shape != (5, 5) or not np.all(np.isfinite(cov)):
        raise ValueError("invalid orbital covariance")
    if not np.allclose(cov, cov.T, rtol=1e-12, atol=0):
        raise ValueError("asymmetric orbital covariance")
    cov = (cov + cov.T) / 2
    np.linalg.cholesky(cov)
    # Verify that the exporter and fitter agree on units and normalization.
    reported_sigma = np.array([float(p[4]) for p in parameters])
    reported_sigma[3] *= math.pi / 180
    if not np.allclose(reported_sigma, np.sqrt(cov.diagonal()), rtol=1e-10, atol=0):
        raise ValueError("exported covariance and reported errors disagree")
    return cov, {"total_fit": int(lines[1].split()[1]),
                 "augmented_chi2": float(lines[2].split()[1]),
                 "dof": int(lines[3].split()[1]),
                 "orbital_values_tcb_par_units": {p[1]: p[3] for p in parameters}}


def correlation(cov):
    sigma = np.sqrt(cov.diagonal())
    return cov / np.outer(sigma, sigma)


def verify_design(prefix, native):
    """Check covariance with column-scaled SciPy QR at the exported linearization.

    This shares TEMPO2 derivatives and priors; it is not an independent timing
    model. The exporter rebuilds derivatives at post-fit values, explaining
    small differences from the covariance of the preceding fit iteration.
    """
    from scipy.linalg import qr, solve_triangular
    prefix = str(prefix)
    n, p, nc = map(int, Path(prefix + ".meta").read_text().splitlines()[0].split())
    design = np.fromfile(prefix + ".design").reshape(n, p)
    constraints = np.fromfile(prefix + ".constraints").reshape(nc, p)
    toas = np.fromfile(prefix + ".toas").reshape(n, 2)
    indices = np.loadtxt(prefix + ".orbit", dtype=int)
    matrix = np.vstack([design / toas[:, 1, None], constraints])
    scales = np.linalg.norm(matrix, axis=0)
    if np.any(scales == 0):
        raise ValueError("unconstrained zero design column")
    matrix /= scales
    triangular = qr(matrix, mode="r", check_finite=False)[0][:p, :p]
    z = solve_triangular(triangular.T, np.eye(p)[:, indices], lower=True)
    cov = (z.T @ z) / np.outer(scales[indices], scales[indices])
    ds = float(np.max(abs(np.sqrt(cov.diagonal()/native.diagonal())-1)))
    dc = float(np.max(abs(correlation(cov)-correlation(native))))
    # These are explicit relinearization tolerances, not floating-point claims.
    if ds > 0.002 or dc > 0.001:
        raise ValueError(f"independent QR check failed: sigma {ds}, correlation {dc}")
    return {"method": "column-scaled SciPy QR of TEMPO2 post-fit design and priors",
            "max_relative_sigma_difference": ds, "max_absolute_correlation_difference": dc,
            "sigma_tolerance": 0.002, "correlation_tolerance": 0.001,
            "covariance_tcb_internal": cov.tolist()}


def make_search_artifact(par, fit_export, *, ephemeris, clock_dir, fit_par):
    raw, fit = read_native_export(fit_export)
    if fit["total_fit"] != 3284 or fit["dof"] != 11592:
        raise ValueError("unexpected fit dimensions for the pinned J0437 release")
    sigma = np.sqrt(raw.diagonal())
    canonical = raw * np.outer(CANONICAL_SCALE, CANONICAL_SCALE)
    np.linalg.cholesky(canonical)
    native_keys = ("PB_d", "A1_light_s", "ECC", "OM_deg", "T0_d")
    par_sigma = sigma * np.array([1, 1, 1, 180/math.pi, 1])
    rows = [line.split() for line in Path(fit_par).read_text().splitlines()]
    fit.update(tempo2_version="2024.12.1", solver="native QR", arguments=["-qrfit", "-norescale"],
               ephemeris=ephemeris, n_toas=11637, n_timing_coefficients=45,
               n_noise_coefficients=3239, released_tres_us=0.698,
               reproduced_tres_us=float(next(r[1] for r in rows if r and r[0]=="TRES")),
               frozen_unmatched_jump_numbers=sorted(UNMATCHED_JUMPS))
    return {"schema_version": 2, "target": "J0437-4715", "parameter_order": list(ORBITAL_ORDER),
            "canonical_units": ["s", "light-s", "1", "rad", "s"], "epoch_convention": "gps",
            "covariance_kind": "native TEMPO2 QR fixed-hyperparameter GLS covariance",
            "normalization": "unscaled (-norescale); no published-marginal replacement",
            "covariance": canonical.tolist(), "correlation": correlation(raw).tolist(),
            "covariance_tcb_internal": raw.tolist(),
            "reported_sigma_tcb_native": dict(zip(native_keys, published_errors(par).values())),
            "reproduced_gls_sigma_tcb_native": dict(zip(native_keys, par_sigma.tolist())),
            "tcb_to_tdb_linear_scale": TCB_TO_TDB, "fit": fit,
            "provenance": {"source_par_sha256": sha256(par), "source_tim_sha256": TIM_SHA256,
                           "fit_export_sha256": sha256(fit_export), "fit_par_sha256": sha256(fit_par),
                           "output_plugin_sha256": sha256(HERE / "scripts/timing/export_qr_covariance.C"),
                           "clock_sha256": {name: sha256(Path(clock_dir)/name) for name in CLOCK_FILES},
                           "upstream_release": "CSIRO collection 62929, DOI 10.25919/20rx-5f63",
                           "notes": "Orbital block of the full simultaneous timing/noise covariance. Noise hyperparameters fixed. Output plugin does not modify fitting. Covariance is not an author-supplied data product. Fitted central values are recorded in native TCB par units; adopting this covariance does not silently change the search's central ephemeris."}}


def verify_search_artifact(candidate, reference_path):
    reference = json.loads(Path(reference_path).read_text())
    for key in ("normalization", "parameter_order", "canonical_units", "epoch_convention"):
        if candidate[key] != reference[key]:
            raise ValueError(f"reference mismatch: {key}")
    for key in ("solver", "ephemeris", "arguments"):
        if candidate["fit"][key] != reference["fit"][key]:
            raise ValueError(f"fit convention mismatch: {key}")
    for key in ("source_par_sha256", "source_tim_sha256", "clock_sha256", "output_plugin_sha256"):
        if candidate["provenance"][key] != reference["provenance"][key]:
            raise ValueError(f"provenance mismatch: {key}")
    # Compare normalized matrices: a large T0 variance must not hide ECC errors.
    a, b = np.array(candidate["covariance"]), np.array(reference["covariance"])
    scale = np.sqrt(np.outer(b.diagonal(), b.diagonal()))
    if np.max(abs((a-b)/scale)) > 1e-7:
        raise ValueError("covariance differs from reference beyond numerical tolerance")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env", type=Path, default=HERE / ".conda")
    parser.add_argument("--par", type=Path, default=HERE / "data/J0437-4715_reardon2024_dr3.par")
    parser.add_argument("--tim", type=Path, default=HERE / "data/timing/j0437_reardon2024_dr3/J0437-4715.tim")
    parser.add_argument("--clock-dir", type=Path, default=HERE / "data/timing/j0437_reardon2024_dr3/clock")
    parser.add_argument("--work-dir", type=Path, default=HERE / "results/j0437_timing_covariance")
    parser.add_argument("--ephemeris", choices=["DE440", "DE436"], default="DE440")
    parser.add_argument("--verify-design", action="store_true")
    parser.add_argument("--adopt", action="store_true", help="replace the production covariance after checks")
    parser.add_argument("--check", action="store_true", help="compare against the production covariance")
    args = parser.parse_args()
    for key in ("env", "par", "tim", "clock_dir", "work_dir"):
        setattr(args, key, getattr(args, key).resolve())
    require_hash(args.par, PAR_SHA256); require_hash(args.tim, TIM_SHA256)
    for name in CLOCK_FILES:
        if not (args.clock_dir/name).is_file():
            raise ValueError(f"missing release clock file: {name}")
    tempo2 = args.env / "bin/tempo2"
    version = subprocess.run([str(tempo2), "-v"], capture_output=True, text=True)
    # This TEMPO2 release prints a valid version but returns status 1.
    if version.stdout.strip() != "2024.12.1":
        raise ValueError(f"expected TEMPO2 2024.12.1, got {version.stdout!r}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    source = HERE / "scripts/timing/export_qr_covariance.C"
    plugin = args.work_dir / "qrcov_linux-gnu_plug.t2"
    subprocess.run(["g++", "-std=c++17", "-O2", "-fPIC", "-shared", f"-I{args.env/'include'}",
                    str(source), f"-L{args.env/'lib'}", f"-Wl,-rpath,{args.env/'lib'}",
                    "-ltempo2", "-o", str(plugin)], check=True)
    fixed_par = args.work_dir / "input.par"
    freeze_unmatched_jumps(args.par, fixed_par)
    lines = fixed_par.read_text().splitlines()
    fixed_par.write_text("\n".join(f"EPHEM {args.ephemeris}" if l.split()[:1]==["EPHEM"] else l
                                  for l in lines) + "\n")
    export = args.work_dir / "covariance.tsv"
    env = dict(os.environ, TEMPO2=str(args.env/"share/tempo2"), AUDIT_COV_OUTPUT=str(export),
               OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2", TEMPO2_PLUG_PATH=str(args.work_dir))
    if args.verify_design:
        env["QR_DESIGN_PREFIX"] = str(args.work_dir/"model")
    else:
        env.pop("QR_DESIGN_PREFIX", None)
    command = [str(tempo2), "-clkdir", str(args.clock_dir), "-f", str(fixed_par), str(args.tim),
               "-qrfit", "-norescale", "-output", "qrcov", "-tcheck"]
    print("+", " ".join(command), flush=True)
    with (args.work_dir/"run.log").open("w") as log:
        subprocess.run(command, cwd=args.work_dir, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    candidate = make_search_artifact(args.par, export, ephemeris=args.ephemeris,
                                    clock_dir=args.clock_dir, fit_par=args.work_dir/"native.par")
    if args.verify_design:
        check = verify_design(args.work_dir/"model", np.array(candidate["covariance_tcb_internal"]))
        (args.work_dir/"independent_qr_check.json").write_text(json.dumps(check, indent=2)+"\n")
        print("Independent QR:", {k:v for k,v in check.items() if k != "covariance_tcb_internal"})
    candidate_path = args.work_dir/"j0437_reardon2024_covariance.candidate.json"
    candidate_path.write_text(json.dumps(candidate, indent=2)+"\n")
    if args.check:
        verify_search_artifact(candidate, REFERENCE)
    if args.adopt:
        REFERENCE.write_text(candidate_path.read_text())
    print(f"Covariance: {candidate_path}; adopted={args.adopt}")


if __name__ == "__main__":
    main()
