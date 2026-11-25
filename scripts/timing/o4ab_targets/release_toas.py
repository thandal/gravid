"""Reconstruct the eight original TOA selections from release files and Portal exports."""

import argparse
import hashlib
import json
import tarfile
import zipfile
from pathlib import Path

from paths import HERE, INPUT
from prepare_inputs import normalize_release

O4B = (60410.625, 60703 + 17 / 24)


def timing_file(rows):
    return ("FORMAT 1\n" + "\n".join(rows) + "\n").encode()


def portal_zip(path, stream, dates=None):
    rows = []
    with zipfile.ZipFile(path) as archive:
        for name in sorted(archive.namelist()):
            if not name.endswith(".tim") or Path(name).parent.name != stream:
                continue
            for line in archive.read(name).decode().splitlines():
                words = line.split()
                if len(words) < 5 or words[0] in {
                    "#",
                    "C",
                    "FORMAT",
                    "MODE",
                    "INCLUDE",
                }:
                    continue
                if dates is None or dates[0] <= float(words[2]) <= dates[1]:
                    rows.append(line)
    return timing_file(dict.fromkeys(rows))


def portal_rows(path):
    rows = json.loads(path.read_text())
    usable = [
        r
        for r in rows
        if not r["observation"]["restricted"]
        and not r["dmCorrected"]
        and r["nsubType"] == "1"
        and r["obsNchan"] == 32
        and r["obsNpol"] == 1
        and r["mjdErr"] > 0
        and O4B[0] <= float(r["mjd"]) <= O4B[1]
    ]

    def signature(stream):
        return {
            (r["mjd"], r["freqMhz"], r["mjdErr"], r["tmplt"])
            for r in usable
            if r["project"]["short"] == stream
        }

    if signature("PTA") != signature("PTA2"):
        raise ValueError(f"PTA and PTA2 differ in {path}")
    selected = sorted(
        (r for r in usable if r["project"]["short"] == "PTA"),
        key=lambda r: (float(r["mjd"]), r["freqMhz"]),
    )
    flags = "fe be f bw tobs tmplt gof nbin nch chan rcvr snr length subint".split()
    lines = []
    for r in selected:
        observation = r["observation"]
        day = observation["utcStart"].replace("+00:00", "").replace(":", "-")
        line = f"api_{day}_b{observation['beam']} {r['freqMhz']:.12g} {r['mjd']} {r['mjdErr']:.12g} {r['telescope']}"
        lines.append(
            line
            + "".join(f" -{key} {r[key]}" for key in flags if r.get(key) is not None)
        )
    return timing_file(lines)


def lofar(path):
    with tarfile.open(path) as bundle:
        tim = bundle.extractfile(
            "dvt+20/J1400-1438/J1400-1438_LOFAR_8chan_outliers_rejected.tim"
        )
        rows = [
            words
            for line in tim.read().decode().splitlines()
            if len(words := line.split()) > 5
        ]
    return timing_file(" ".join(words) + f" -tel {words[4]}" for words in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "products",
        type=Path,
        help="Directory containing the original releases and Portal exports",
    )
    args = parser.parse_args()
    src = args.products
    selections = dict(normalize_release(src / "partim.tar.gz"))
    for name in ["J1045-4509", "J1737-0811"]:
        selections[name + "_o4b.tim"] = portal_rows(src / (name + "_api_toas.json"))
    selections["J1231-1411_o4b.tim"] = portal_zip(
        src / "J1231-1411_toas.zip", "RelBin", O4B
    )
    selections["J1745-0952_historical.tim"] = portal_zip(
        src / "J1745-0952_toas.zip", "PTA"
    )
    selections["J1400-1431_training.tim"] = lofar(src / "dvt+20.tar")
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    for name, data in selections.items():
        if hashlib.sha256(data).hexdigest() != manifest[name]["sha256"]:
            raise ValueError(
                f"The release products give a different selection for {name}"
            )
    INPUT.mkdir(parents=True, exist_ok=True)
    for name, data in selections.items():
        (INPUT / name).write_bytes(data)
    print(f"Reconstructed and verified {len(selections)} TOA files in {INPUT}")


if __name__ == "__main__":
    main()
