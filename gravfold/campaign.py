"""Integrity contract for the J0437 O4a heterodyne campaign.

This module deliberately has no import-time dependency on h5py.  The normal
gravfold package stays lightweight; product verification imports h5py only when
it is actually requested on a heterodyne/harvest machine.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Iterable

CHUNK_SECONDS = 5 * 86400
O4A_SLICES = (
    (1368195220, 1375282153),
    (1375282153, 1382369086),
    (1382369086, 1389456018),
)
DETECTORS = ("H1", "L1")
PRODUCT_RE = re.compile(
    r"^het_J0437-4715_(H1|L1)_2_g(\d+)\.hdf5$"
)

_TEXT_PHASE_KEYS = ("PSRJ", "EPHEM", "CLK", "UNITS", "BINARY")
_NUMERIC_PHASE_KEYS = (
    "PMRA", "PMDEC", "PX", "POSEPOCH", "F0", "F1", "PEPOCH",
    "PB", "A1", "TASC", "EPS1", "EPS2", "PBDOT",
)
_SEXAGESIMAL_PHASE_KEYS = ("RAJ", "DECJ")


class ProductVerificationError(ValueError):
    """A heterodyne product violates the campaign integrity contract."""


@dataclass(frozen=True)
class ProductRecord:
    basename: str
    detector: str
    gps_start: int
    gps_end: int
    sha256: str
    size_bytes: int
    samples: int
    time_segments: int
    livetime_samples: int
    coverage_attested: bool
    canonical_par_sha256: str

    def to_dict(self) -> dict:
        return asdict(self)


def expected_product_ranges(
    detectors: Iterable[str] = DETECTORS,
) -> dict[str, tuple[str, int, int]]:
    """Return the exact basename -> (detector, start, end) campaign plan."""
    out: dict[str, tuple[str, int, int]] = {}
    for detector in detectors:
        if detector not in DETECTORS:
            raise ValueError(f"unsupported detector {detector!r}")
        for slice_start, slice_end in O4A_SLICES:
            for start in range(slice_start, slice_end, CHUNK_SECONDS):
                end = min(start + CHUNK_SECONDS, slice_end)
                name = f"het_J0437-4715_{detector}_2_g{start}.hdf5"
                out[name] = (detector, start, end)
    return out


def expected_product_names() -> set[str]:
    return set(expected_product_ranges())


def product_name(detector: str, gps_start: int) -> str:
    name = f"het_J0437-4715_{detector}_2_g{int(gps_start)}.hdf5"
    if name not in expected_product_ranges((detector,)):
        raise ValueError(f"{detector}:g{gps_start} is not on the O4a campaign grid")
    return name


def _parse_par(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text.splitlines():
        fields = line.split()
        if fields and not fields[0].startswith("#") and len(fields) >= 2:
            values[fields[0]] = fields[1]
    return values


def _sexagesimal(value: str) -> tuple[int, Decimal, Decimal, Decimal]:
    fields = value.strip().split(":")
    if len(fields) != 3:
        raise ProductVerificationError(f"invalid sexagesimal value {value!r}")
    sign = -1 if fields[0].startswith("-") else 1
    try:
        first = Decimal(fields[0].lstrip("+-"))
        return sign, first, Decimal(fields[1]), Decimal(fields[2])
    except InvalidOperation as exc:
        raise ProductVerificationError(
            f"invalid sexagesimal value {value!r}"
        ) from exc


def _canonical_phase(values: dict[str, str]) -> dict[str, str]:
    missing = [
        key for key in (*_TEXT_PHASE_KEYS, *_NUMERIC_PHASE_KEYS,
                        *_SEXAGESIMAL_PHASE_KEYS)
        if key not in values
    ]
    if missing:
        raise ProductVerificationError(
            "phase model is missing: " + ", ".join(missing)
        )
    canonical: dict[str, str] = {}
    for key in _TEXT_PHASE_KEYS:
        canonical[key] = values[key]
    for key in _NUMERIC_PHASE_KEYS:
        try:
            canonical[key] = str(Decimal(values[key]).normalize())
        except InvalidOperation as exc:
            raise ProductVerificationError(
                f"invalid numeric phase value {key}={values[key]!r}"
            ) from exc
    for key in _SEXAGESIMAL_PHASE_KEYS:
        sign, a, b, c = _sexagesimal(values[key])
        canonical[key] = f"{sign}:{a.normalize()}:{b.normalize()}:{c.normalize()}"
    return canonical


def canonical_par_sha256(par_text: str) -> str:
    payload = json.dumps(
        _canonical_phase(_parse_par(par_text)), sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _require_matching_phase(embedded: str, expected: str) -> str:
    have = _canonical_phase(_parse_par(embedded))
    want = _canonical_phase(_parse_par(expected))
    bad = [key for key in want if have.get(key) != want[key]]
    if bad:
        details = ", ".join(
            f"{key}={have.get(key)!r} (want {want[key]!r})" for key in bad
        )
        raise ProductVerificationError("embedded phase model mismatch: " + details)
    return canonical_par_sha256(expected)


def _sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def verify_product(
    path: str | os.PathLike,
    parfile: str | os.PathLike,
    *,
    require_attested_coverage: bool = False,
) -> ProductRecord:
    """Open and fully verify one J0437 campaign HDF5 product."""
    try:
        import h5py
        import numpy as np
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("product verification requires h5py and numpy") from exc

    path = Path(path)
    match = PRODUCT_RE.fullmatch(path.name)
    if match is None:
        raise ProductVerificationError(
            f"unexpected product basename {path.name!r}"
        )
    detector, gps_text = match.groups()
    gps_start = int(gps_text)
    expected = expected_product_ranges()
    if path.name not in expected:
        raise ProductVerificationError(f"product is not on campaign grid: {path.name}")
    _, _, gps_end = expected[path.name]

    expected_par = Path(parfile).read_text()
    try:
        with h5py.File(path, "r") as h5:
            if "HeterodynedData" not in h5:
                raise ProductVerificationError("missing /HeterodynedData")
            data = h5["HeterodynedData"]
            if not isinstance(data, h5py.Dataset):
                raise ProductVerificationError(
                    "/HeterodynedData is not the expected dataset"
                )
            if data.ndim != 1 or data.size == 0:
                raise ProductVerificationError("heterodyne data is empty or non-1D")
            sample_count = int(data.size)
            # Read in bounded blocks so a large campaign product does not need a
            # second full in-memory copy merely for validation.
            for offset in range(0, data.size, 1_000_000):
                block = data[offset:offset + 1_000_000]
                if not bool(np.all(np.isfinite(block))):
                    raise ProductVerificationError(
                        f"nonfinite sample(s) in [{offset}, {offset + len(block)})"
                    )

            attrs = data.attrs
            if str(attrs.get("detector", "")) != detector:
                raise ProductVerificationError("filename/detector metadata mismatch")
            if not math.isclose(float(attrs.get("dt", math.nan)), 1.0,
                                rel_tol=0.0, abs_tol=1e-12):
                raise ProductVerificationError("dt is not 1 second")
            if not math.isclose(float(attrs.get("freqfactor", math.nan)), 2.0,
                                rel_tol=0.0, abs_tol=1e-12):
                raise ProductVerificationError("frequency factor is not 2")
            for key, want in (
                ("include_ssb", True), ("include_bsb", True),
                ("include_glitch", False),
            ):
                if bool(attrs.get(key, not want)) is not want:
                    raise ProductVerificationError(f"unexpected {key} metadata")
            if "par" not in attrs:
                raise ProductVerificationError("missing embedded par metadata")
            par_digest = _require_matching_phase(str(attrs["par"]), expected_par)

            if "HeterodynedDataTimeSegments" not in h5:
                raise ProductVerificationError("missing time-segment dataset")
            time_segments = np.asarray(h5["HeterodynedDataTimeSegments"][:],
                                       dtype=float)
            if (time_segments.ndim != 2 or time_segments.shape[1:] != (2,)
                    or len(time_segments) == 0
                    or not bool(np.all(np.isfinite(time_segments)))
                    or not bool(np.all(time_segments[:, 1] >= time_segments[:, 0]))):
                raise ProductVerificationError("invalid time-segment dataset")
            # cwinpy stores one sample at both endpoints of every output segment.
            expected_samples = int(round(float(np.sum(
                time_segments[:, 1] - time_segments[:, 0]
            )))) + len(time_segments)
            if expected_samples != sample_count:
                raise ProductVerificationError(
                    f"sample/segment mismatch: {sample_count} != {expected_samples}"
                )

            attested = bool(attrs.get("gravfold_acquisition_complete", False))
            if attested:
                for key, want in (
                    ("gravfold_requested_start", gps_start),
                    ("gravfold_requested_end", gps_end),
                ):
                    if int(attrs.get(key, -1)) != want:
                        raise ProductVerificationError(f"invalid {key} attestation")
                if "GravfoldAcquisitionSegments" not in h5:
                    raise ProductVerificationError(
                        "coverage attested but acquisition segments are absent"
                    )
                acquired = np.asarray(h5["GravfoldAcquisitionSegments"][:],
                                      dtype=float)
                if (acquired.ndim != 2 or acquired.shape[1:] != (2,)
                        or len(acquired) != len(time_segments)
                        or not bool(np.all(np.isfinite(acquired)))):
                    raise ProductVerificationError(
                        "invalid acquisition-segment attestation"
                    )
                if (not bool(np.all(time_segments[:, 0] >= acquired[:, 0]))
                        or not bool(np.all(time_segments[:, 1] <= acquired[:, 1]))):
                    raise ProductVerificationError(
                        "output time segments exceed acquired segments"
                    )
            if require_attested_coverage and not attested:
                raise ProductVerificationError(
                    "product predates fail-closed acquisition attestation"
                )
    except OSError as exc:
        raise ProductVerificationError(f"cannot open HDF5: {exc}") from exc

    return ProductRecord(
        basename=path.name,
        detector=detector,
        gps_start=gps_start,
        gps_end=gps_end,
        sha256=_sha256(path),
        size_bytes=path.stat().st_size,
        samples=sample_count,
        time_segments=len(time_segments),
        livetime_samples=expected_samples,
        coverage_attested=attested,
        canonical_par_sha256=par_digest,
    )


def load_ledger(path: str | os.PathLike) -> dict[str, dict]:
    ledger_path = Path(path)
    if not ledger_path.exists():
        return {}
    payload = json.loads(ledger_path.read_text())
    if not isinstance(payload, dict) or payload.get("schema") != 1:
        raise ValueError(f"unsupported verification ledger: {ledger_path}")
    records = payload.get("products", {})
    if not isinstance(records, dict):
        raise ValueError(f"invalid product map in {ledger_path}")
    return records


def write_ledger(path: str | os.PathLike, records: dict[str, dict]) -> None:
    """Atomically replace the digest-bound campaign verification ledger."""
    ledger_path = Path(path)
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": 1,
        "expected_products": len(expected_product_names()),
        "products": {key: records[key] for key in sorted(records)},
    }
    temporary = ledger_path.with_name(
        f".{ledger_path.name}.tmp.{os.getpid()}"
    )
    with temporary.open("w") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, ledger_path)
