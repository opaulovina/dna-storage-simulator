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
same file size (140 sequences of 120 nt), I:D:S ratios, error rates and sequencing depths:

```bash
python experiments/validate_xie.py data/the_grandmother.txt
```

The raw results, the run environment and a heatmap are written to
`results/validation_xie.csv`, `results/validation_xie.json` and `results/validation_xie.png`.
