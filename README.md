# gravfold

Code and results for coherent continuous gravitational-wave searches of
binary pulsars, published article by article.

## Articles

| Article | Folder |
|---|---|
| One orbital template for a coherent continuous-wave search of PSR J0437−4715 | [`articles/j0437-one-orbital-template/`](articles/j0437-one-orbital-template/) |

Each article folder holds the text, the PDF, the figures and a README that
maps every number in the article to the script that computes it and the
saved file that holds it.

## Layout

| Path | Contents |
|---|---|
| `articles/` | one folder per article |
| `gravfold/` | the Python library the analyses import |
| `scripts/` | search scripts, and under `scripts/timing/` one folder per timing analysis with its code and small saved results |
| `results/` | saved search outputs |
| `tests/` | tests of the library |

Bulk outputs (`results/timing/`), downloaded radio data and the reduced
strain are not in the repository. Each article's README says where they come
from.

## Install

Python 3.10 or later.

```bash
pip install -e ".[dev]"            # library and tests
pip install -e ".[dev,cwinpy]"     # also what the strain searches need
pytest
```

The timing scripts also need TEMPO2. They expect it in a conda environment
at `.conda/` in the repository root.

## Licence

The code is under the BSD 3-Clause licence (`LICENSE`). The article text and
figures under `articles/` are under the Creative Commons Attribution 4.0
International licence (`articles/LICENSE`).
