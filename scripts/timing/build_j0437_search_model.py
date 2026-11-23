#!/usr/bin/env python3
"""Convert the adopted native DE440 fit to a TDB T2 gravitational-wave model.

PINT supplies the time-scale conversion. Its internal DDK angular convention
is converted back to TEMPO2's convention for LAL's T2 implementation. Radio
propagation, instrumental offsets and stochastic timing-noise terms are not
part of the continuous-wave phase model.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np
import pint
import pint.logging
from pint.models import get_model

ROOT = Path(__file__).resolve().parents[2]
KEYS = set("PSRJ RAJ DECJ F0 F1 PEPOCH POSEPOCH PMRA PMDEC PX BINARY PB T0 "
           "A1 OM ECC PBDOT OMDOT M2 KOM KIN UNITS EPHEM CLK".split())


def build(source, covariance):
    raw = source.read_text()
    artifact = json.loads(covariance.read_text())
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != artifact["provenance"]["fit_par_sha256"]:
        raise ValueError("search mean must come from the adopted covariance's fit")
    rows = {line.split()[0]: line for line in raw.splitlines() if line.split()}
    text = "\n".join(rows[key] for key in rows if key in KEYS) + "\n"
    model = get_model(io.StringIO(text), allow_tcb=True, allow_T2=True)
    if model.BINARY.value != "DDK" or model.EPHEM.value != "DE440":
        raise ValueError("expected the adopted DE440 T2/DDK model")
    lines = ["# Native DE440 QR fit, converted TCB -> TDB with PINT.",
             "# T2 binary terms retained; KIN/KOM use TEMPO2's IAU convention.",
             "# Original heterodyne reference is a separate, unchanged file."]
    for line in model.as_parfile().splitlines():
        key = line.split()[0] if line.split() else ""
        if key not in KEYS:
            continue
        if key == "BINARY":
            line = "BINARY T2"
        elif key in ("KIN", "KOM", "RAJ", "DECJ"):
            line = rows[key]  # angles and native sexagesimal uncertainty units unchanged
        lines.append(line.rstrip())
    # TEMPO2 accepts 'SINI KIN'; LAL needs the derived numerical value.
    sini = np.sin(np.deg2rad(float(rows["KIN"].split()[1])))
    lines.extend([f"SINI {sini:.17g}", "PLANET_SHAPIRO N"])
    output = "\n".join(lines) + "\n"
    metadata = dict(
        source_fit_sha256=digest,
        covariance_sha256=hashlib.sha256(covariance.read_bytes()).hexdigest(),
        search_par_sha256=hashlib.sha256(output.encode()).hexdigest(),
        pint_version=pint.__version__, binary_model="T2", units="TDB",
        retained_binary_terms=["PB", "A1", "ECC", "OM", "T0", "PBDOT",
                               "OMDOT", "M2", "SINI", "KIN", "KOM"],
        geometry="TEMPO2 IAU KIN/KOM; PINT DDK uses transformed angles internally",
        note="No noise-hyperparameter refit or covariance rescaling.",
    )
    return output, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "scripts/timing/qr_solve/fit.par")
    parser.add_argument("--covariance", type=Path, default=ROOT / "gravfold/resources/j0437_de440_binary_covariance.json")
    parser.add_argument("--out", type=Path, default=ROOT / "gravfold/resources/j0437_de440_search.par")
    args = parser.parse_args()
    pint.logging.setup(level="ERROR")
    output, metadata = build(args.source, args.covariance)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(output)
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
