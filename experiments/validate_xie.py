import argparse
import csv
import itertools
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import BYTES_PER_SEQUENCE, encode
from dna_storage.correction.mafft import MafftCorrector
from dna_storage.experiment import environment, label
from dna_storage.metrics import mean_hamming, recovery_rate

SEQUENCES = 140
FILE_SIZE = SEQUENCES * BYTES_PER_SEQUENCE
REPETITIONS = 3
ERROR_RATES = (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30)
DEPTHS = (5, 10, 20, 30, 50, 70, 100, 200)
IDS_RATIOS = ((1, 1, 1), (1, 1, 2))
FIELDS = ("ids_ratio", "error_rate", "depth", "seed", "sequences", "recovery_rate", "hamming")


def run(sequences, ids_ratio, error_rate, depth, seed):
    config = ChannelConfig(error_rate=error_rate, depth=depth, seed=seed, ids_ratio=ids_ratio)
    recovered = MafftCorrector().correct(simulate_reads(sequences, config))
    return {
        "ids_ratio": label(ids_ratio),
        "error_rate": error_rate,
        "depth": depth,
        "seed": seed,
        "sequences": len(sequences),
        "recovery_rate": recovery_rate(recovered, sequences),
        "hamming": mean_hamming(recovered, sequences),
    }


def load(path):
    text = path.read_bytes()
    if not text:
        raise SystemExit(f"{path} está vazio")
    copies = FILE_SIZE // len(text) + 1
    return (text * copies)[:FILE_SIZE]


def mean_rates(rows, field="recovery_rate"):
    groups = {}
    for row in rows:
        group = (row["ids_ratio"], row["error_rate"], row["depth"])
        groups.setdefault(group, []).append(row[field])
    return {group: mean(values) for group, values in groups.items()}


def summarize(rows, error_rates, depths):
    rates = mean_rates(rows)
    distances = mean_rates(rows, "hamming")
    for ratio in dict.fromkeys(row["ids_ratio"] for row in rows):
        print(f"\nI:D:S {ratio}")
        for error_rate in error_rates:
            full = [depth for depth in depths if rates[(ratio, error_rate, depth)] == 1]
            first = full[0] if full else "-"
            last = rates[(ratio, error_rate, depths[-1])]
            accuracy = 1 - distances[(ratio, error_rate, depths[-1])]
            print(f"  {error_rate:>4.0%}: 100% a partir de {first!s:>4}; ", end="")
            print(f"em {depths[-1]}: {last:.1%} (bases certas: {accuracy:.1%})")


def plot(rows, error_rates, depths, path):
    rates = mean_rates(rows)
    ratios = list(dict.fromkeys(row["ids_ratio"] for row in rows))
    figure, axes = plt.subplots(1, len(ratios), figsize=(6 * len(ratios), 6), squeeze=False)
    for axis, ratio in zip(axes[0], ratios, strict=True):
        grid = [[rates[(ratio, rate, depth)] for depth in depths] for rate in error_rates]
        axis.imshow(grid, cmap="RdBu", vmin=0, vmax=1, aspect="auto")
        for row, values in enumerate(grid):
            for column, value in enumerate(values):
                axis.text(column, row, f"{value:.2f}", ha="center", va="center")
        axis.set_xticks(range(len(depths)), [str(depth) for depth in depths])
        axis.set_yticks(range(len(error_rates)), [f"{rate:.0%}" for rate in error_rates])
        axis.set_xlabel("profundidade")
        axis.set_ylabel("taxa de erro por base")
        axis.set_title(f"MAFFT, I:D:S {ratio}")
    figure.tight_layout()
    figure.savefig(path, dpi=150)


def main():
    parser = argparse.ArgumentParser(description="Condição A nos parâmetros de Xie et al. (2023)")
    parser.add_argument("file", type=Path)
    parser.add_argument("--repetitions", type=int, default=REPETITIONS)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--error-rates", type=float, nargs="+", default=list(ERROR_RATES))
    parser.add_argument("--depths", type=int, nargs="+", default=list(DEPTHS))
    parser.add_argument("--output", type=Path, default=Path("results/validation_xie"))
    args = parser.parse_args()

    data = load(args.file)
    sequences = encode(data)
    jobs = list(
        itertools.product(IDS_RATIOS, args.error_rates, args.depths, range(args.repetitions))
    )

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run, sequences, *job) for job in jobs]
        for done, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            print(f"\r{done}/{len(jobs)}", end="", flush=True)
    rows.sort(key=lambda row: (row["ids_ratio"], row["error_rate"], row["depth"], row["seed"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output.with_suffix(".csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "file": args.file.name,
        **environment(data),
        "error_rates": args.error_rates,
        "depths": args.depths,
        "ids_ratios": [list(ratio) for ratio in IDS_RATIOS],
        "repetitions": args.repetitions,
    }
    with open(args.output.with_suffix(".json"), "w") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
    plot(rows, args.error_rates, args.depths, args.output.with_suffix(".png"))
    summarize(rows, args.error_rates, args.depths)


if __name__ == "__main__":
    main()
