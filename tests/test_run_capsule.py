import json
import os

import pytest

from gravfold.run_capsule import (build_capsule, capsule_digest,
                                  chunk_manifest_digest, load_capsule,
                                  require_equal_capsules,
                                  shard_capsule_fields, source_bundle_sha256,
                                  verify_shard_capsule, worker_attestation)


def _tree(tmp_path):
    (tmp_path / "gravfold").mkdir()
    (tmp_path / "scripts").mkdir()
    (tmp_path / "gravfold" / "a.py").write_text("x = 1\n")
    (tmp_path / "scripts" / "b.py").write_text("y = 2\n")
    return tmp_path


def test_source_bundle_is_deterministic_and_mtime_invariant(tmp_path):
    root = _tree(tmp_path)
    first, count = source_bundle_sha256(root)
    assert count == 2
    os.utime(root / "gravfold" / "a.py", (0, 0))
    again, _ = source_bundle_sha256(root)
    assert again == first
    (root / "gravfold" / "a.py").write_text("x = 3\n")
    changed, _ = source_bundle_sha256(root)
    assert changed != first


def test_capsule_digest_excludes_timestamp_and_detects_tampering(tmp_path):
    root = _tree(tmp_path)
    inputs = {"chunk_a.npz": "aa", "chunk_b.npz": "bb"}
    one = build_capsule(root, ["prog", "--x"], input_hashes=inputs)
    two = build_capsule(root, ["prog", "--x"], input_hashes=inputs,
                        created_utc="2026-07-23T00:00:00+00:00")
    assert one["capsule_sha256"] == two["capsule_sha256"]
    other = build_capsule(root, ["prog", "--y"], input_hashes=inputs)
    assert other["capsule_sha256"] != one["capsule_sha256"]
    path = tmp_path / "capsule.json"
    path.write_text(json.dumps(two))
    assert load_capsule(path)["capsule_sha256"] == two["capsule_sha256"]
    two["builder_invocation"] = ["prog", "--tampered"]
    path.write_text(json.dumps(two))
    with pytest.raises(ValueError, match="digest mismatch"):
        load_capsule(path)


def test_load_capsule_detects_edited_input_hash(tmp_path):
    # input_hashes ride outside the digested identity fields for size; the
    # loader must recompute the manifest digest so an edited entry cannot
    # hide behind a still-valid capsule SHA (Codex finding 4).
    capsule = build_capsule(_tree(tmp_path), ["prog"],
                            input_hashes={"a.npz": "aa", "b.npz": "bb"})
    capsule["input_hashes"]["a.npz"] = "ee"
    path = tmp_path / "capsule.json"
    path.write_text(json.dumps(capsule))
    with pytest.raises(ValueError, match="input-hash manifest mismatch"):
        load_capsule(path)


def test_verify_shard_capsule_rejects_forged_json(tmp_path):
    capsule = build_capsule(_tree(tmp_path), ["prog"],
                            input_hashes={"a.npz": "aa"})
    fields = shard_capsule_fields(capsule)
    assert verify_shard_capsule(fields["capsule_json"],
                                fields["capsule_sha256"])
    forged = json.loads(fields["capsule_json"])
    forged["builder_image"] = "forged:image"
    with pytest.raises(ValueError, match="digest mismatch"):
        verify_shard_capsule(json.dumps(forged), fields["capsule_sha256"])
    # A self-consistent capsule from a DIFFERENT run must not pass either.
    other_root = tmp_path / "o"
    other_root.mkdir()
    other = build_capsule(_tree(other_root), ["other"])
    with pytest.raises(ValueError, match="does not match"):
        verify_shard_capsule(json.dumps(other), fields["capsule_sha256"])


def test_non_clean_phase_requires_a_phase_plan(tmp_path):
    root = _tree(tmp_path)
    with pytest.raises(ValueError, match="requires a phase plan"):
        build_capsule(root, ["prog"], phase="injection")
    plan = tmp_path / "jobs.json"
    plan.write_text('{"tiles": {}}')
    capsule = build_capsule(root, ["prog"], phase="injection",
                            phase_plan_path=plan)
    assert capsule["phase"] == "injection"
    assert capsule["phase_plan_sha256"]


def test_chunk_manifest_digest_is_order_independent():
    assert chunk_manifest_digest({"a": "1", "b": "2"}) == \
        chunk_manifest_digest({"b": "2", "a": "1"})
    assert chunk_manifest_digest({"a": "1"}) != \
        chunk_manifest_digest({"a": "2"})


def test_require_equal_capsules_fails_closed():
    # Capsule-less sets need an explicit legacy opt-in (Codex finding 5).
    with pytest.raises(ValueError, match="capsule-less"):
        require_equal_capsules([("s0", None), ("s1", None)])
    assert require_equal_capsules([("s0", None), ("s1", None)],
                                  allow_missing=True) is None
    assert require_equal_capsules([("s0", "d"), ("s1", "d")]) == "d"
    with pytest.raises(ValueError, match="mismatch"):
        require_equal_capsules([("s0", "d"), ("s1", None)],
                               allow_missing=True)
    with pytest.raises(ValueError, match="mismatch"):
        require_equal_capsules([("s0", "d"), ("s1", "e")])


def test_worker_attestation_records_execution_identity():
    record = worker_attestation(["run_gpu_search.py", "--shard", "0"])
    assert record["environment_lock_sha256"]
    assert record["search_invocation"][-1] == "0"
    assert record["image_reference"]


def test_capsule_digest_matches_identity_fields(tmp_path):
    capsule = build_capsule(_tree(tmp_path), ["prog"])
    assert capsule["capsule_sha256"] == capsule_digest(capsule)
