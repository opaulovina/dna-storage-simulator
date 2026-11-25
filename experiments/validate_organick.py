import argparse
import itertools
import json
import random
from pathlib import Path

import pandas as pd

from dna_storage.codec.payload import encode
from dna_storage.correction.reed_solomon import layout
from dna_storage.experiment import Trial, environment, label
from dna_storage.validation import (
    COVERAGE_MODES,
    ORGANICK_ERROR_RATE,
    ORGANICK_IDS_RATIO,
    coverage_run,
    frontier_cell,
)


def main():
    parser = argparse.ArgumentParser(description="Validação da Condição B frente a Organick")
    parser.add_argument("file", type=Path)
    parser.add_argument(
        "--coverages", type=int, nargs="+", default=[1, 2, 3, 4, 5, 6, 8, 10, 12, 14]
    )
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument("--redundancy", type=float, default=0.15)
    parser.add_argument("--reconstruction", default="BMA")
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--output", type=Path, default=Path("results/organick"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    if args.redundancy <= 0:
        raise SystemExit("a validação precisa de redundância maior que 0")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "file": args.file.name,
        **environment(data),
        "error_rate": ORGANICK_ERROR_RATE,
        "ids_ratio": label(ORGANICK_IDS_RATIO),
        "redundancy": args.redundancy,
        "reconstruction": args.reconstruction,
    }

    runs = []
    for mode, coverage, seed in itertools.product(
        COVERAGE_MODES, args.coverages, range(args.repetitions)
    ):
        trial = Trial(
            ORGANICK_ERROR_RATE,
            coverage,
            seed,
            ORGANICK_IDS_RATIO,
            args.redundancy,
            args.reconstruction,
        )
        runs.append(coverage_run(data, trial, mode))
        print(f"\rcobertura {len(runs)}", end="", flush=True)
    print()
    coverage_path = args.output.parent / f"{args.output.name}_cobertura"
    pd.DataFrame(runs).to_csv(coverage_path.with_suffix(".csv"), index=False)
    coverage_path.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False)
    )

    data_count = len(encode(data))
    parity = len(layout(data_count, args.redundancy)[0][1])
    rng = random.Random(0)
    cells = []
    for errors in range(parity // 2 + 2):
        for erasures in range(parity + 2):
            cells.append(frontier_cell(data_count, parity, errors, erasures, args.trials, rng))
        print(f"\rfronteira {errors + 1}/{parity // 2 + 2}", end="", flush=True)
    print()
    frontier_path = args.output.parent / f"{args.output.name}_fronteira"
    pd.DataFrame(cells).to_csv(frontier_path.with_suffix(".csv"), index=False)
    frontier_path.with_suffix(".json").write_text(
        json.dumps({**metadata, "trials": args.trials, "seed": 0}, indent=2, ensure_ascii=False)
    )
    print(f"resultados em {coverage_path}.csv e {frontier_path}.csv")


if __name__ == "__main__":
    main()
