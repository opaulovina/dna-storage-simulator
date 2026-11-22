import argparse
import csv
import itertools
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import psutil

from dna_storage.correction.mafft import MafftCorrector
from dna_storage.correction.reed_solomon import RsHybridCorrector
from dna_storage.experiment import CONDITIONS, SITES, Trial, environment, failure_counts, label

FIELDS = (
    "condition",
    "ids_ratio",
    "error_rate",
    "depth",
    "seed",
    "redundancy",
    "reconstruction",
    "sequences",
    "recovery_rate",
    *SITES,
)


def ratio(text):
    parts = tuple(float(part) for part in text.split(":"))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"proporção I:D:S inválida: {text}")
    return parts


def key(condition, ids_ratio, error_rate, depth, seed, redundancy):
    return condition, ids_ratio, float(error_rate), int(depth), int(seed), float(redundancy)


def build_jobs(args):
    jobs = []
    for ids_ratio, error_rate, depth, seed in itertools.product(
        args.ratios, args.error_rates, args.depths, range(args.repetitions)
    ):
        if MafftCorrector.name in args.conditions:
            jobs.append((Trial(error_rate, depth, seed, ids_ratio, 0.0), MafftCorrector.name))
        if RsHybridCorrector.name not in args.conditions:
            continue
        for redundancy in args.redundancy:
            trial = Trial(error_rate, depth, seed, ids_ratio, redundancy, "BMA")
            jobs.append((trial, RsHybridCorrector.name))
    return jobs


def job_key(trial, condition):
    redundancy = 0.0 if condition == MafftCorrector.name else trial.redundancy
    return key(
        condition, label(trial.ids_ratio), trial.error_rate, trial.depth, trial.seed, redundancy
    )


def done_keys(path):
    if not path.exists():
        return set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise SystemExit(f"{path} tem colunas de uma versão anterior; use outro --output")
        return {key(*(row[field] for field in FIELDS[:6])) for row in reader}


def main():
    parser = argparse.ArgumentParser(description="Onde ficam as falhas: índice, payload ou ambos")
    parser.add_argument("file", type=Path)
    parser.add_argument("--error-rates", type=float, nargs="+", required=True)
    parser.add_argument("--depths", type=int, nargs="+", required=True)
    parser.add_argument("--ratios", type=ratio, nargs="+", default=[(1.0, 1.0, 1.0)])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--redundancy", type=float, nargs="+", default=[0.0, 0.15])
    parser.add_argument("--conditions", nargs="+", choices=CONDITIONS, default=list(CONDITIONS))
    parser.add_argument("--workers", type=int, default=psutil.cpu_count(logical=False))
    parser.add_argument("--output", type=Path, default=Path("results/falhas"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    csv_path = args.output.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"file": args.file.name, **environment(data), "redundancy": args.redundancy}
    json_path = args.output.with_suffix(".json")
    if json_path.exists() and json.loads(json_path.read_text())["sha256"] != metadata["sha256"]:
        raise SystemExit(f"{json_path} foi gerado com outro arquivo de entrada; use outro --output")
    json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False))

    finished = done_keys(csv_path)
    jobs = [job for job in build_jobs(args) if job_key(*job) not in finished]
    print(f"{len(finished)} execuções já feitas; {len(jobs)} a fazer")

    new_file = not csv_path.exists()
    with (
        open(csv_path, "a", newline="") as handle,
        ProcessPoolExecutor(max_workers=args.workers) as pool,
    ):
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if new_file:
            writer.writeheader()
        futures = [pool.submit(failure_counts, data, *job) for job in jobs]
        for count, future in enumerate(as_completed(futures), start=1):
            writer.writerow(future.result())
            handle.flush()
            print(f"\r{count}/{len(jobs)}", end="", flush=True)
    print()


if __name__ == "__main__":
    main()
