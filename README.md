# dna-storage-simulator

[![tests](https://github.com/opaulovina/dna-storage-simulator/actions/workflows/tests.yml/badge.svg)](https://github.com/opaulovina/dna-storage-simulator/actions/workflows/tests.yml)

Interactive DNA Storage Error Simulator: A comparative analysis between MAFFT and Reed-Solomon algorithms for data recovery in synthetic DNA.

## Requirements

- Python 3.12+
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) 7.505 (on Debian/Ubuntu: `sudo apt-get install mafft`)

## Development

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
ruff check .
```

## Validation against Xie et al. (2023)

Condition A (MAFFT) is checked against the behaviour reported by Xie et al. (2023), using the
same sequence layout (120 nt: 8 nt index + 112 nt payload), I:D:S ratios, error rates and
sequencing depths. The input is the English abstract of the monograph (`data/abstract.txt`):

```bash
python experiments/validate_xie.py data/abstract.txt
```

The raw results, the run environment and a heatmap are written to
`results/validation_xie.csv`, `results/validation_xie.json` and `results/validation_xie.png`.

## Design decisions

Every scope and design decision, with its source, is listed (in Portuguese) in
[`docs/ESCOPO.txt`](docs/ESCOPO.txt).
