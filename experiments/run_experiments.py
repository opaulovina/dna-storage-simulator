import argparse
import csv
import itertools
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import psutil

from dna_storage.correction.mafft import MafftCorrector
from dna_storage.correction.reed_solomon import RECONSTRUCTIONS
from dna_storage.experiment import (
    CONDITIONS,
    Trial,
    environment,
    label,
    run_condition,
    series_name,
)

FIELDS = (
    "condition",
    "ids_ratio",
    "error_rate",
    "depth",
    "seed",
    "redundancy",
    "reconstruction",
    "sequences",
    "molecules",
    "density",
    "recovery_rate",
    "hamming",
    "seconds",
    "memory_kib",
    "workers",
)


def ratio(text):
    parts = tuple(float(part) for part in text.split(":"))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"proporção I:D:S inválida: {text}")
    return parts


def key(condition, ids_ratio, error_rate, depth, seed, redundancy, reconstruction):
    return (
        condition,
        ids_ratio,
        float(error_rate),
        int(depth),
        int(seed),
        float(redundancy),
        reconstruction,
    )


def job_key(trial, condition):
    reconstruction = "-" if condition == MafftCorrector.name else trial.reconstruction
    return key(
        condition,
        label(trial.ids_ratio),
        trial.error_rate,
        trial.depth,
        trial.seed,
        trial.redundancy,
        reconstruction,
    )


def done_keys(path):
    if not path.exists():
        return set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise SystemExit(f"{path} tem colunas de uma versão anterior; use outro --output")
        return {key(*(row[field] for field in FIELDS[:7])) for row in reader}


def build_jobs(args):
    jobs = []
    for ids_ratio, error_rate, depth, seed in itertools.product(
        args.ratios, args.error_rates, args.depths, range(args.repetitions)
    ):
        if MafftCorrector.name in args.conditions:
            trial = Trial(error_rate, depth, seed, ids_ratio, args.redundancy[0])
            jobs.append((trial, MafftCorrector.name))
        for condition in args.conditions:
            if condition == MafftCorrector.name:
                continue
            for redundancy, reconstruction in itertools.product(
                args.redundancy, args.reconstructions
            ):
                trial = Trial(error_rate, depth, seed, ids_ratio, redundancy, reconstruction)
                jobs.append((trial, condition))
    return jobs


def check_input(path, metadata):
    if path.exists():
        previous = json.loads(path.read_text())
        if previous["sha256"] != metadata["sha256"]:
            raise SystemExit(f"{path} foi gerado com outro arquivo de entrada; use outro --output")
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False))


def summarize(path):
    table = pd.read_csv(path)
    table["series"] = [
        series_name(*values)
        for values in zip(
            table["condition"], table["reconstruction"], table["redundancy"], strict=True
        )
    ]
    summary = table.pivot_table(
        index=["ids_ratio", "error_rate"],
        columns=["series", "depth"],
        values="recovery_rate",
        aggfunc="mean",
    )
    print(summary.round(2).to_string())


# roda a grade pedida nas condições escolhidas, com "--workers" processos ao mesmo tempo (o padrão
# é um por núcleo físico), grava cada execução no CSV assim que termina e pula as que já foram
# feitas
def main():
    parser = argparse.ArgumentParser(description="Condições A e B sobre o mesmo ruído")
    parser.add_argument("file", type=Path)
    parser.add_argument("--error-rates", type=float, nargs="+", required=True)
    parser.add_argument("--depths", type=int, nargs="+", required=True)
    parser.add_argument("--ratios", type=ratio, nargs="+", default=[(1.0, 1.0, 1.0)])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--redundancy", type=float, nargs="+", default=[0.15])
    parser.add_argument(
        "--reconstructions", nargs="+", choices=list(RECONSTRUCTIONS), default=["BMA"]
    )
    parser.add_argument("--conditions", nargs="+", choices=CONDITIONS, default=list(CONDITIONS))
    parser.add_argument("--workers", type=int, default=psutil.cpu_count(logical=False))
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
        "reconstructions": args.reconstructions,
    }
    csv_path = args.output.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    check_input(args.output.with_suffix(".json"), metadata)

    finished = done_keys(csv_path)
    jobs = [job for job in build_jobs(args) if job_key(*job) not in finished]
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
            writer.writerow({**future.result(), "workers": args.workers})
            handle.flush()
            print(f"\r{count}/{len(jobs)}", end="", flush=True)
    print()
    summarize(csv_path)


if __name__ == "__main__":
    main()
