"""File layout for the O4ab binary-target analysis."""

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WORK = Path(os.environ.get("O4AB_WORK", ROOT / "results/timing/o4ab_targets")).resolve()
TEMPO = Path(os.environ.get("O4AB_TEMPO2_PREFIX", ROOT / ".conda")).resolve()
INPUT = WORK / "inputs"
TOOLS = WORK / "bin"
BASE = WORK / "fits"
OUT = WORK / "noise"
LOFAR = WORK / "lofar"
NAMES = ["J1045-4509", "J1737-0811", "J1231-1411"]
