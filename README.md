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

## Validation of Condition B

```bash
python experiments/validate_organick.py data/abstract_3920.txt
python experiments/positional_errors.py data/abstract_3920.txt
```

`validate_organick.py` runs BMA + RS 15% at the Illumina error regime measured by Organick et
al. (2018) (0.64% per base, I:D:S 1:5:10) with 1 to 14 reads per sequence and 10 seeds, and
records the "used error resilience" of every matrix row, (2 x errors + erasures) / parity. Reads
are drawn either with a fixed depth (the channel of the main experiments) or with a Poisson
number of reads per molecule, so that some molecules get no reads and become erasures. It also
checks the decoder directly in GF(2^16): for every combination of 0 to 11 errors and 0 to 22
erasures in a 161-symbol codeword with 21 parity symbols, it counts how often the message comes
back right, how often the decoder reports a failure and how often it returns a wrong message.
Output: `results/organick_cobertura.csv` and `results/organick_fronteira.csv`.

`positional_errors.py` counts, without Reed-Solomon, how often each of the 120 positions of the
reconstructed strand is wrong for BMA, Double-sided BMA and NW-POA (12% and 16% error, 10 and
20 reads, 10 seeds), to compare with the error propagation described by Sharma et al. (2024).
Output: `results/perfil_posicional.csv`.

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

## MAFFT with a fixed consensus length

Condition A follows Xie et al. (2023): the consensus drops every column where the gap wins, so it
can end up shorter or longer than 120 nt, while the BMA always returns 120 nt. This diagnostic
aligns each cluster once with MAFFT and scores both the free consensus and a consensus forced to
120 columns (the columns with the fewest gaps, the same rule as NW-POA):

```bash
python experiments/mafft_length.py data/abstract_3920.txt \
    --error-rates 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.25 0.30 \
    --depths 5 10 20 30 50 100 --ratios 1:1:1 1:1:2 --repetitions 3 \
    --output results/mafft_120
```

## Analysis

```bash
python experiments/analyze.py
```

Reads the text grids (`final_texto`, `ablacao`, `reconstrucoes`), the image grid, the failure
sites, the fixed-length MAFFT diagnostic and the Condition B validation, and writes tables (CSV) and figures (PDF and PNG) to `results/analise/`: recovery per
depth, full-recovery limits, Reed-Solomon gain, reconstruction variants, minimum sequencing cost
and time per sequence. The findings are summarised (in Portuguese) in section 15 of
[`docs/ESCOPO.txt`](docs/ESCOPO.txt).

## Interactive app

```bash
streamlit run app/app.py
```

Upload any file up to 1.75 MiB (the 8-base index addresses at most 4^8 = 65,536 sequences of
28 bytes) and choose the error rate, sequencing depth, I:D:S weights, seed, the reconstruction
used by Condition B (BMA, Double-sided BMA or NW-POA) and the Reed-Solomon redundancy. The app
runs both conditions on the same noisy reads, shows the metrics side by side, maps which
sequences failed, and previews and lets you download the original and both recovered files. A
second tab plots the results of the batch experiments. With parity, Condition B needs more
molecules than data sequences; files that would exceed the 65,536 addresses for the chosen
redundancy are refused with an explanation.

## Condition B as in Organick et al. (2018)

Since version 0.2.0 the outer code uses 16-bit symbols (GF(2^16)) and one long codeword per file,
and the decoder places each reconstructed strand by the address it reads: strands with an
address out of range or claimed twice are dropped, and every empty column is passed to the
Reed-Solomon decoder as an erasure; a molecule with no reads leaves its column empty. Results produced with version 0.1.0 (8-bit symbols, strands
placed by the simulator) are not comparable for Condition B.

## Design decisions

Every scope and design decision, with its source, is listed (in Portuguese) in
[`docs/ESCOPO.txt`](docs/ESCOPO.txt).
