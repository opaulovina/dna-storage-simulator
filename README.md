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
sequencing depths. The input is the English abstract of the monograph (`data/abstract.txt`),
repeated up to 3,920 bytes so that it gives the same 140 sequences used by Xie et al.
(`data/abstract_3920.txt`):

```bash
python experiments/validate_xie.py data/abstract_3920.txt
```

The raw results, the run environment and a heatmap are written to
`results/validation_xie.csv`, `results/validation_xie.json` and `results/validation_xie.png`.

## Comparison with Xie et al. figures

```bash
python experiments/compare_xie.py
```

Compares `results/validation_xie.csv` with values transcribed by hand from Figures 1 and 2a of
Xie et al. (2023), stored in `data/xie_figuras.csv`.

## Pilot run

Measures, on the current machine, the cost of each condition, whether running several
experiments in parallel changes the measured time, and how long the final grid would take:

```bash
python experiments/pilot.py data/abstract_3920.txt
```

## Running the experiments

Each run is executed in a fresh process so that time and memory are measured in isolation.
Results are appended to a CSV as they finish, so an interrupted run can be resumed by repeating
the command. The final grids were run one at a time (`--workers 1`), because running in parallel
inflated the measured time by 65-148% on the experiment machine (see `results/pilot.json`).

Main comparison (text and image):

```bash
python experiments/run_experiments.py data/abstract_3920.txt \
    --error-rates 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.25 0.30 \
    --depths 5 10 20 30 50 100 --ratios 1:1:1 1:1:2 --repetitions 3 \
    --workers 1 --output results/final_texto

python experiments/run_experiments.py data/imagem.pgm \
    --error-rates 0.04 0.08 0.12 0.16 0.20 --depths 10 30 50 --repetitions 3 \
    --workers 1 --output results/final_imagem
```

Ablation of the outer Reed-Solomon code (Condition B only, same grid):

```bash
python experiments/run_experiments.py data/abstract_3920.txt \
    --error-rates 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.25 0.30 \
    --depths 5 10 20 30 50 100 --ratios 1:1:1 1:1:2 --repetitions 3 \
    --conditions "RS híbrido" --redundancy 0 0.05 0.10 0.15 0.25 0.50 \
    --workers 1 --output results/ablacao
```

Reconstruction variants of Condition B (Double-sided BMA and NW-POA):

```bash
python experiments/run_experiments.py data/abstract_3920.txt \
    --error-rates 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.25 0.30 \
    --depths 5 10 20 30 50 100 --ratios 1:1:1 1:1:2 --repetitions 3 \
    --conditions "RS híbrido" --redundancy 0 0.15 --reconstructions DBMA NW-POA \
    --workers 1 --output results/reconstrucoes
```

Each CSV has one row per condition and run; the matching JSON stores the input checksum, the
hardware and the software versions.

## Failure sites

Re-runs the grid and classifies every failed sequence as an index error (first 8 nt), a payload
error or both. Time is not measured here, so it can run in parallel:

```bash
python experiments/failure_sites.py data/abstract_3920.txt \
    --error-rates 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.25 0.30 \
    --depths 5 10 20 30 50 100 --ratios 1:1:1 1:1:2 --repetitions 3 \
    --redundancy 0 0.05 0.15 0.50 --output results/falhas
```

## Analysis

```bash
python experiments/analyze.py
```

Reads the text grids (`final_texto`, `ablacao`, `reconstrucoes`), the image grid and the failure
sites, and writes tables (CSV) and figures (PDF and PNG) to `results/analise/`: recovery per
depth, full-recovery limits, Reed-Solomon gain, reconstruction variants, minimum sequencing cost
and time per sequence. The findings are summarised (in Portuguese) in section 15 of
[`docs/ESCOPO.txt`](docs/ESCOPO.txt).

## Interactive app

```bash
streamlit run app/app.py
```

Upload any file up to 1.75 MB (the 8-base index addresses at most 4^8 = 65,536 sequences of
28 bytes) and choose the error rate, sequencing depth, I:D:S weights, seed, the reconstruction
used by Condition B (BMA, Double-sided BMA or NW-POA) and the Reed-Solomon redundancy. The app
runs both conditions on the same noisy reads, shows the metrics side by side, maps which
sequences failed, and previews and lets you download the original and both recovered files. A
second tab plots the results of the batch experiments, and a third shows the validation against
Xie et al. (2023).

## Design decisions

Every scope and design decision, with its source, is listed (in Portuguese) in
[`docs/ESCOPO.txt`](docs/ESCOPO.txt).
