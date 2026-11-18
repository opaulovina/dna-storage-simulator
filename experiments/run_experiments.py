import argparse
import csv
import itertools
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from dna_storage.experiment import CONDITIONS, Trial, environment, label, run_condition

FIELDS = (
    "condition",
    "ids_ratio",
    "error_rate",
    "depth",
    "seed",
    "redundancy",
    "sequences",
    "molecules",
    "density",
    "recovery_rate",
    "hamming",
    "seconds",
    "memory_kib",
)


def ratio(text):
    parts = tuple(float(part) for part in text.split(":"))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"proporção I:D:S inválida: {text}")
    return parts


def key(condition, ids_ratio, error_rate, depth, seed, redundancy):
    return (condition, ids_ratio, float(error_rate), int(depth), int(seed), float(redundancy))


def done_keys(path):
    if not path.exists():
        return set()
    with open(path, newline="") as handle:
        return {key(*(row[field] for field in FIELDS[:6])) for row in csv.DictReader(handle)}


def check_input(path, metadata):
    if path.exists():
        previous = json.loads(path.read_text())
        if previous["sha256"] != metadata["sha256"]:
            raise SystemExit(f"{path} foi gerado com outro arquivo de entrada; use outro --output")
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False))


def summarize(path):
    table = pd.read_csv(path)
    summary = table.pivot_table(
        index=["ids_ratio", "error_rate"],
        columns=["condition", "depth"],
        values="recovery_rate",
        aggfunc="mean",
    )
    print(summary.round(2).to_string())


def main():
    parser = argparse.ArgumentParser(description="Condições A e B sobre o mesmo ruído")
    parser.add_argument("file", type=Path)
    parser.add_argument("--error-rates", type=float, nargs="+", required=True)
    parser.add_argument("--depths", type=int, nargs="+", required=True)
    parser.add_argument("--ratios", type=ratio, nargs="+", default=[(1.0, 1.0, 1.0)])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--redundancy", type=float, default=0.15)
    parser.add_argument("--conditions", nargs="+", choices=CONDITIONS, default=list(CONDITIONS))
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--output", type=Path, default=Path("results/experiments"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    metadata = {
        "file": args.file.name,
        **environment(data),
        "workers": args.workers,
        "redundancy": args.redundancy,
    }
    csv_path = args.output.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    check_input(args.output.with_suffix(".json"), metadata)

    finished = done_keys(csv_path)
    jobs = [
        (Trial(error_rate, depth, seed, ids_ratio, args.redundancy), condition)
        for ids_ratio, error_rate, depth, seed, condition in itertools.product(
            args.ratios, args.error_rates, args.depths, range(args.repetitions), args.conditions
        )
    ]
    jobs = [
        (trial, condition)
        for trial, condition in jobs
        if key(
            condition,
            label(trial.ids_ratio),
            trial.error_rate,
            trial.depth,
            trial.seed,
            trial.redundancy,
        )
        not in finished
    ]
    print(f"{len(finished)} execuções já feitas; {len(jobs)} a fazer")

    new_file = not csv_path.exists()
    with (
        open(csv_path, "a", newline="") as handle,
        ProcessPoolExecutor(max_workers=args.workers, max_tasks_per_child=1) as pool,
    ):
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
        futures = [pool.submit(run_condition, data, trial, condition) for trial, condition in jobs]
        for count, future in enumerate(as_completed(futures), start=1):
            writer.writerow(future.result())
            handle.flush()
            print(f"\r{count}/{len(jobs)}", end="", flush=True)
    print()
    summarize(csv_path)


if __name__ == "__main__":
    main()
