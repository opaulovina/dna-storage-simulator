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
