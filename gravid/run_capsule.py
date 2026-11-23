"""Run-capsule identity for search executions.

A capsule binds one search run to the exact bytes that define it: a
canonical source archive digest, git commit + dirty-patch digest, the
execution environment (container image digest or an explicit uncontainerized
identity, plus a package lock), the command line, the orbital posterior, and
a content-hash manifest of every strain input.  Its digest is embedded in
every shard a run produces; merge and adoption refuse to combine shards
whose capsules differ.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import tarfile
from pathlib import Path

SOURCE_GLOBS = ("gravid/**/*.py", "scripts/**/*.py", "scripts/**/*.C", "compute/**/*.sh",
                "pyproject.toml")
# Fields hashed into the capsule digest; created_utc and the digest itself
# are deliberately outside so rebuilding an identical run reproduces the
# digest bit-exactly.  The capsule is the RUN CONTRACT (what code, inputs,
# posterior, and phase plan a run must execute); the environment that
# actually executed a shard is attested separately by the WORKER at run time
# (worker_attestation) — a controller-built capsule must never masquerade as
# an execution-environment record (Codex capsule review, finding 1).
IDENTITY_FIELDS = ("source_bundle_sha256", "git_commit", "dirty_patch_sha256",
                   "builder_image", "builder_environment_lock_sha256",
                   "builder_environment", "builder_invocation",
                   "phase", "phase_plan_sha256",
                   "posterior_sha256", "chunk_manifest_sha256", "input_count")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_bundle_sha256(repo_root):
    """Digest of a canonical tar of the source tree.

    Deterministic by construction: sorted paths, zeroed timestamps and
    ownership, normalized modes.  Includes untracked files, so that code
    which ran without being committed is still identified.
    """
    repo_root = Path(repo_root)
    files = sorted(
        {p.relative_to(repo_root).as_posix()
         for pattern in SOURCE_GLOBS for p in repo_root.glob(pattern)
         if p.is_file()})
    digest = hashlib.sha256()

    class _Hash(io.RawIOBase):
        def writable(self):
            return True

        def write(self, data):
            digest.update(data)
            return len(data)

    # Stream mode ("w|"): sequential writes only, no seeks — required for a
    # hashing sink, and it keeps the archive format canonical.
    with tarfile.open(fileobj=_Hash(), mode="w|") as archive:
        for name in files:
            path = repo_root / name
            data = path.read_bytes()
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            info.mtime = 0
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(data))
    return digest.hexdigest(), len(files)


def git_identity(repo_root):
    """Read-only git identity: HEAD commit and a digest of the dirty diff."""
    def run(*argv):
        try:
            proc = subprocess.run(
                ["git", "-C", str(repo_root), *argv], capture_output=True,
                timeout=60)
            return proc.stdout if proc.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            return None
    head = run("rev-parse", "HEAD")
    diff = run("diff", "HEAD", "--binary")
    return (
        head.decode().strip() if head else None,
        hashlib.sha256(diff).hexdigest() if diff is not None else None,
    )


def environment_identity():
    """Package lock digest plus the interpreter/platform summary."""
    from importlib import metadata
    entries = sorted(
        f"{dist.metadata['Name']}=={dist.version}"
        for dist in metadata.distributions()
        if dist.metadata["Name"])
    lock = hashlib.sha256("\n".join(entries).encode()).hexdigest()
    summary = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
    for name in ("numpy", "scipy", "cupy"):
        try:
            summary[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    return lock, summary


def chunk_manifest_digest(input_hashes):
    """Digest of the canonical basename -> sha256 input mapping."""
    canonical = json.dumps(dict(sorted(input_hashes.items())),
                           sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def capsule_digest(capsule):
    identity = {key: capsule.get(key) for key in IDENTITY_FIELDS}
    canonical = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def build_capsule(repo_root, builder_invocation, posterior_path=None,
                  input_hashes=None, created_utc=None, phase="clean_search",
                  phase_plan_path=None):
    """Assemble a run-contract capsule.

    ``input_hashes`` maps basename -> sha256.  ``phase`` names the campaign
    stage the capsule authorizes; any phase other than ``clean_search``
    requires a phase plan (the jobs file), whose content hash becomes part of
    the capsule identity so a job lookup can never silently fall back to a
    different search than the one this capsule contracted.
    """
    bundle_sha, bundle_files = source_bundle_sha256(repo_root)
    commit, dirty_sha = git_identity(repo_root)
    lock_sha, environment = environment_identity()
    input_hashes = dict(input_hashes or {})
    if phase != "clean_search" and phase_plan_path is None:
        raise ValueError(f"phase {phase!r} requires a phase plan (jobs file)")
    capsule = {
        "schema": "gravid.run_capsule.v2",
        "source_bundle_sha256": bundle_sha,
        "source_bundle_file_count": bundle_files,
        "source_bundle_globs": list(SOURCE_GLOBS),
        "git_commit": commit,
        "dirty_patch_sha256": dirty_sha,
        # Builder context, honestly labeled: this is where the CAPSULE was
        # made, never where a shard ran.
        "builder_image": (os.environ.get("GRAVID_IMAGE_DIGEST")
                          or f"uncontainerized:{platform.node()}"),
        "builder_environment_lock_sha256": lock_sha,
        "builder_environment": environment,
        "builder_invocation": list(builder_invocation),
        "phase": phase,
        "phase_plan_sha256": (sha256_file(phase_plan_path)
                              if phase_plan_path else None),
        "posterior_sha256": (sha256_file(posterior_path)
                             if posterior_path else None),
        "chunk_manifest_sha256": (chunk_manifest_digest(input_hashes)
                                  if input_hashes else None),
        "input_count": len(input_hashes),
        "input_hashes": input_hashes,
    }
    if created_utc is not None:
        capsule["created_utc"] = created_utc
    capsule["capsule_sha256"] = capsule_digest(capsule)
    return capsule


def verify_capsule_dict(capsule, origin="capsule"):
    """Full integrity check of a capsule dict; returns it or raises.

    Verifies the self-declared digest AND that the embedded input_hashes
    reproduce the digested chunk_manifest_sha256 — the hashes ride outside
    the identity fields for size, so without this recomputation a recorded
    input hash could be edited while the capsule SHA stayed valid (Codex
    capsule review, finding 4).
    """
    expected = capsule.get("capsule_sha256")
    actual = capsule_digest(capsule)
    if expected != actual:
        raise ValueError(
            f"{origin} digest mismatch: recorded {expected}, recomputed "
            f"{actual} — the capsule was edited after creation")
    hashes = capsule.get("input_hashes") or {}
    declared = capsule.get("chunk_manifest_sha256")
    if hashes or declared:
        recomputed = chunk_manifest_digest(hashes) if hashes else None
        if recomputed != declared:
            raise ValueError(
                f"{origin} input-hash manifest mismatch: declared "
                f"{declared}, recomputed {recomputed} — an embedded input "
                "hash was altered")
    if int(capsule.get("input_count", 0)) != len(hashes):
        raise ValueError(f"{origin} input_count does not match input_hashes")
    return capsule


def load_capsule(path):
    """Load a capsule file with full integrity verification."""
    return verify_capsule_dict(json.loads(Path(path).read_text()),
                               origin=f"capsule {path}")


def verify_shard_capsule(capsule_json_str, expected_sha):
    """Verify a shard's embedded capsule JSON against its declared digest.

    The standalone capsule_sha256 string and the capsule_json payload are
    separate npz members; consumers must never trust the JSON without
    proving it reproduces the declared digest (Codex capsule review,
    finding 4: a forged capsule_json was accepted while the declared SHA
    stayed unchanged).
    """
    capsule = verify_capsule_dict(json.loads(capsule_json_str),
                                  origin="shard capsule_json")
    if capsule["capsule_sha256"] != expected_sha:
        raise ValueError(
            f"shard capsule_json digest {capsule['capsule_sha256']} does "
            f"not match the shard's declared capsule_sha256 {expected_sha}")
    return capsule


def shard_capsule_fields(capsule):
    """The two arrays embedded in every shard npz."""
    return {
        "capsule_sha256": capsule["capsule_sha256"],
        "capsule_json": json.dumps(capsule, sort_keys=True),
    }


def worker_attestation(argv):
    """Execution-environment record captured BY the process that ran a shard.

    This — not the capsule — is what says which environment executed the
    search: the worker's own package lock, platform, image reference (a
    mutable tag is a reference, not a digest; the environment lock hash is
    the content-based identity), hostname, and full search argv.
    """
    lock_sha, environment = environment_identity()
    return {
        "schema": "gravid.worker_attestation.v1",
        "environment_lock_sha256": lock_sha,
        "environment": environment,
        "image_reference": (os.environ.get("GRAVID_IMAGE_DIGEST")
                            or f"uncontainerized:{platform.node()}"),
        # What the machine was launched to do.  A product can then say what kind of machine produced it rather
        # than leaving it to be inferred from which flags happened to be set.
        "role": os.environ.get("GRAVID_ROLE") or "unrecorded",
        "hostname": platform.node(),
        "search_invocation": list(argv),
    }


def require_equal_capsules(labels_and_digests, allow_missing=False):
    """Refuse to combine shards whose capsule digests differ or are absent.

    Accepts an iterable of (label, digest-or-None).  A fully capsule-less
    set is refused unless ``allow_missing`` is explicitly set (products
    that predate capsules only); a MIXED set is always an error.
    """
    rows = list(labels_and_digests)
    digests = {digest for _label, digest in rows}
    if digests == {None}:
        if allow_missing:
            return None
        raise ValueError(
            "capsule-less shard set refused: only legacy clean-search "
            "products may opt in explicitly (B3SF_ALLOW_NO_CAPSULE=1 / "
            "--allow-uncapsuled)")
    if None in digests or len(digests) != 1:
        detail = ", ".join(f"{label}={digest or 'MISSING'}"
                           for label, digest in rows)
        raise ValueError(f"shard run-capsule mismatch: {detail}")
    return digests.pop()


def manifest_provenance(capsule_json_str, expected_sha):
    """The manifest['provenance'] block from a VERIFIED shard capsule."""
    capsule = verify_shard_capsule(capsule_json_str, expected_sha)
    return {
        "capsule_sha256": capsule["capsule_sha256"],
        "source_bundle_sha256": capsule["source_bundle_sha256"],
        "git_commit": capsule["git_commit"],
        "dirty_patch_sha256": capsule["dirty_patch_sha256"],
        "phase": capsule["phase"],
        "phase_plan_sha256": capsule["phase_plan_sha256"],
        "builder_image": capsule["builder_image"],
        "posterior_sha256": capsule["posterior_sha256"],
        "chunk_manifest_sha256": capsule["chunk_manifest_sha256"],
    }
