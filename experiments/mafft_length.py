import argparse
import csv
import itertools
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import psutil

from dna_storage.experiment import Trial, environment, label, mafft_length_counts

FIELDS = (
    "ids_ratio",
    "error_rate",
    "depth",
    "seed",
    "sequences",
    "recovery_rate",
    "recovery_rate_120",
)


def ratio(text):
    parts = tuple(float(part) for part in text.split(":"))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"proporção I:D:S inválida: {text}")
    return parts


def key(ids_ratio, error_rate, depth, seed):
    return ids_ratio, float(error_rate), int(depth), int(seed)


def done_keys(path):
    if not path.exists():
        return set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise SystemExit(f"{path} tem colunas de uma versão anterior; use outro --output")
        return {key(*(row[field] for field in FIELDS[:4])) for row in reader}


def main():
    parser = argparse.ArgumentParser(description="MAFFT com o consenso livre e com 120 colunas")
    parser.add_argument("file", type=Path)
    parser.add_argument("--error-rates", type=float, nargs="+", required=True)
    parser.add_argument("--depths", type=int, nargs="+", required=True)
    parser.add_argument("--ratios", type=ratio, nargs="+", default=[(1.0, 1.0, 1.0)])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--workers", type=int, default=psutil.cpu_count(logical=False))
    parser.add_argument("--output", type=Path, default=Path("results/mafft_120"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    csv_path = args.output.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"file": args.file.name, **environment(data)}
    json_path = args.output.with_suffix(".json")
    if json_path.exists() and json.loads(json_path.read_text())["sha256"] != metadata["sha256"]:
        raise SystemExit(f"{json_path} foi gerado com outro arquivo de entrada; use outro --output")
    json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False))

    finished = done_keys(csv_path)
    trials = [
        Trial(error_rate, depth, seed, ids_ratio)
        for ids_ratio, error_rate, depth, seed in itertools.product(
            args.ratios, args.error_rates, args.depths, range(args.repetitions)
        )
    ]
    trials = [
        trial
        for trial in trials
        if key(label(trial.ids_ratio), trial.error_rate, trial.depth, trial.seed) not in finished
    ]
    print(f"{len(finished)} execuções já feitas; {len(trials)} a fazer")

    new_file = not csv_path.exists()
    with (
        open(csv_path, "a", newline="") as handle,
        ProcessPoolExecutor(max_workers=args.workers) as pool,
    ):
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
        futures = [pool.submit(mafft_length_counts, data, trial) for trial in trials]
        for count, future in enumerate(as_completed(futures), start=1):
            writer.writerow(future.result())
            handle.flush()
            print(f"\r{count}/{len(trials)}", end="", flush=True)
    print()


if __name__ == "__main__":
    main()
