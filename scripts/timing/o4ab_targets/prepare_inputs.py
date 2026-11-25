"""Copy and verify the exact TOAs and clock tables used for the full timing run."""

import argparse
import hashlib
import json
import shutil
import tarfile
from pathlib import Path

from paths import HERE, INPUT, NAMES, TEMPO


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize_release(archive):
    with tarfile.open(archive) as bundle:
        for name in NAMES:
            member = bundle.extractfile(f"partim/{name}.tim")
            rows = [
                line
                for line in member.read().decode().splitlines()
                if len(words := line.split()) > 4
                and words[0] not in {"C", "#", "FORMAT", "MODE"}
            ]
            yield (
                name + "_training.tim",
                ("FORMAT 1\nMODE 1\n" + "\n".join(rows) + "\n").encode(),
            )


def assemble_clocks(clock_dir, manifest):
    runtime = TEMPO / "share/tempo2/clock"
    # Build the common clock directory before deriving the station aliases.
    for entry in sorted(manifest):
        record = manifest[entry]
        if not entry.endswith(".clk"):
            continue
        name = Path(entry).name
        # The low-frequency station aliases use the same GPS-to-UTC series
        # with the relevant TEMPO2 station name in the clock header.
        station = name.removesuffix("2utc.clk")
        if entry.startswith("low_frequency_clock/") and station in {
            "lwa1",
            "lofar",
            "EFlfrhba",
            "JUlfrhba",
            "NDlfrhba",
            "POlfrhba",
            "TBlfrhba",
            "UWlfrhba",
        }:
            if station == "lwa1":
                header = "# UTC(LWA1) UTC\n# GPS-derived archive time; global GPS correction, not a measured LWA clock record."
            else:
                header = f"# UTC({station}) UTC\n# GPS-derived station time; global correction plus fitted constant station offset."
            gps = (INPUT / "clock/gps2utc.clk").read_text().splitlines()
            data = (
                "\n".join(
                    [header]
                    + [line for line in gps if line and not line.startswith("#")]
                )
                + "\n"
            ).encode()
        else:
            candidates = (clock_dir / name, runtime / name)
            source = next(
                (
                    p
                    for p in candidates
                    if p.is_file() and digest(p) == record["sha256"]
                ),
                None,
            )
            if source is None:
                raise ValueError(
                    f"Missing clock table {name} with SHA-256 {record['sha256']}"
                )
            data = source.read_bytes()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(f"Unexpected generated clock table {entry}")
        destination = INPUT / entry
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive", type=Path, help="Root of the archived results/timing tree"
    )
    parser.add_argument(
        "--partim",
        type=Path,
        help="Use the original MeerKAT partim.tar.gz for the three training sets",
    )
    parser.add_argument(
        "--clocks",
        type=Path,
        help="Rebuild the clock directories from the IPTA snapshot and TEMPO2 runtime",
    )
    args = parser.parse_args()
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    INPUT.mkdir(parents=True, exist_ok=True)
    if args.partim:
        for name, data in normalize_release(args.partim):
            if hashlib.sha256(data).hexdigest() != manifest[name]["sha256"]:
                raise ValueError(f"Unexpected MeerKAT release for {name}")
            (INPUT / name).write_bytes(data)
    if args.clocks:
        assemble_clocks(args.clocks, manifest)
    lwa = INPUT / "J1400-1431_lwa.tim"
    if not args.archive and not lwa.exists():
        shutil.copyfile(HERE / "results" / lwa.name, lwa)
    for name, item in manifest.items():
        destination = INPUT / name
        if (
            args.archive
            and not (args.clocks and name.endswith(".clk"))
            and not (args.partim and name in {n + "_training.tim" for n in NAMES})
        ):
            source = args.archive / item["archive_path"]
            if digest(source) != item["sha256"]:
                raise ValueError(f"Input checksum mismatch for {source}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.resolve() != destination.resolve():
                shutil.copyfile(source, destination)
        if not destination.is_file() or digest(destination) != item["sha256"]:
            raise ValueError(f"Missing or incorrect input: {destination}")
    print(f"Verified {len(manifest)} input files in {INPUT}")


if __name__ == "__main__":
    main()
