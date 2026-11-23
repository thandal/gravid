"""Read and identify the frozen J0437 campaign without implicit preprocessing."""
import json
from pathlib import Path

import numpy as np

from .campaign import verify_product
from .run_capsule import sha256_file, chunk_manifest_digest


def load_chunk(path):
    import h5py
    with h5py.File(path, "r") as stream:
        dataset = stream["HeterodynedData"]
        data = dataset[()]
        detector = str(dataset.attrs["detector"])
        dt = float(dataset.attrs["dt"])
        segments = stream["HeterodynedDataTimeSegments"][()]
    # Identical endpoint convention to CWInPy's HDF5 reader. Avoid constructing
    # its noise-estimation objects: this analysis uses fixed per-file variance.
    times = np.concatenate([np.arange(a, b+dt/2, dt) for a, b in segments])
    if len(times) != len(data) or np.any(np.diff(times) <= 0) or not np.all(np.isfinite(data)):
        raise ValueError(f"invalid campaign samples: {path}")
    return detector, times, data, None


def verify_selection(files, reference_par, *, smoke=False):
    ledger_path = Path(reference_par).parent / ".verified.json"
    ledger = json.loads(ledger_path.read_text())
    expected = ledger["products"]
    names = [Path(f).name for f in files]
    if len(set(names)) != len(names):
        raise ValueError("duplicate campaign input")
    if not smoke and (len(expected) != 102 or set(names) != set(expected)):
        raise ValueError("science run requires all 102 verified products; partial inputs require smoke mode")
    hashes, counts, intervals = {}, {}, {}
    for path in files:
        name = Path(path).name
        if name not in expected:
            raise ValueError(f"unrecognized campaign product: {name}")
        record = verify_product(path, reference_par)
        if record.sha256 != expected[name]["sha256"]:
            raise ValueError(f"archived hash mismatch: {name}")
        hashes[name] = record.sha256
        counts[record.detector] = counts.get(record.detector, 0) + record.samples
        detector, times, _, _ = load_chunk(path)
        intervals.setdefault(detector, []).append((float(times[0]), float(times[-1])))
    for detector, spans in intervals.items():
        spans.sort()
        if any(left[1] >= right[0] for left, right in zip(spans, spans[1:])):
            raise ValueError(f"overlapping campaign files for {detector}")
    if not smoke and counts != {"H1": 12955581, "L1": 13411410}:
        raise ValueError(f"campaign sample counts changed: {counts}")
    return dict(scope="smoke" if smoke else "complete archived 102-file campaign",
                input_hashes=hashes, manifest_sha256=chunk_manifest_digest(hashes),
                archived_ledger_sha256=sha256_file(ledger_path), sample_counts=counts,
                reference_par_sha256=sha256_file(reference_par),
                raw_segment_coverage_attested=False)


def freeze_run(output, plan, covariance, repo_root):
    """Persist the exact plan, source snapshot, environment and run identity."""
    import io
    import shutil
    from importlib import metadata
    import subprocess
    import sys
    import tarfile
    from .run_capsule import build_capsule, SOURCE_GLOBS
    output = Path(output)
    plan_path = Path(str(output)+".plan.json")
    plan_path.write_text(json.dumps(plan, indent=2)+"\n")
    inputs = Path(str(output)+".inputs")
    inputs.mkdir(exist_ok=True)
    for name, spec in plan.get("phase_inputs", {}).items():
        source = Path(spec["path"])
        if sha256_file(source) != spec["sha256"]:
            raise ValueError(f"phase input changed during run preparation: {name}")
        shutil.copyfile(source, inputs / name)
    capsule = build_capsule(repo_root, sys.argv, posterior_path=covariance,
                            input_hashes=plan["input_identity"]["input_hashes"],
                            phase="j0437_de440_qr", phase_plan_path=plan_path)
    lock = sorted(f"{dist.metadata['Name']}=={dist.version}" for dist in metadata.distributions()
                  if dist.metadata['Name'])
    Path(str(output)+".environment.txt").write_text("\n".join(lock)+"\n")
    source = Path(str(output)+".source.tar")
    files = sorted({p.relative_to(repo_root).as_posix()
                    for pattern in SOURCE_GLOBS for p in Path(repo_root).glob(pattern) if p.is_file()})
    with tarfile.open(source, "w|") as archive:
        for name in files:
            content = (Path(repo_root)/name).read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(content)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(content))
    if sha256_file(source) != capsule["source_bundle_sha256"]:
        raise ValueError("source snapshot differs from run capsule")
    patch = subprocess.run(["git", "diff", "HEAD", "--binary"], cwd=repo_root,
                           capture_output=True, check=True).stdout
    Path(str(output)+".patch").write_bytes(patch)
    Path(str(output)+".capsule.json").write_text(json.dumps(capsule, indent=2)+"\n")
    return capsule
