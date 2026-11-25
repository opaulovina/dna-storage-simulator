import argparse
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import psutil

from dna_storage.correction.reed_solomon import RECONSTRUCTIONS
from dna_storage.experiment import Trial, environment
from dna_storage.validation import position_errors


def ratio(text):
    parts = tuple(float(part) for part in text.split(":"))
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"proporção I:D:S inválida: {text}")
    return parts


def main():
    parser = argparse.ArgumentParser(description="Erro por posição da fita reconstruída (sem RS)")
    parser.add_argument("file", type=Path)
    parser.add_argument("--error-rates", type=float, nargs="+", default=[0.12, 0.16])
    parser.add_argument("--depths", type=int, nargs="+", default=[10, 20])
    parser.add_argument("--ratio", type=ratio, default=(1.0, 1.0, 1.0))
    parser.add_argument("--repetitions", type=int, default=10)
    parser.add_argument(
        "--reconstructions", nargs="+", choices=list(RECONSTRUCTIONS), default=list(RECONSTRUCTIONS)
    )
    parser.add_argument("--workers", type=int, default=psutil.cpu_count(logical=False))
    parser.add_argument("--output", type=Path, default=Path("results/perfil_posicional"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    trials = [
        Trial(error_rate, depth, seed, args.ratio, 0.0, reconstruction)
        for error_rate, depth, reconstruction, seed in itertools.product(
            args.error_rates, args.depths, args.reconstructions, range(args.repetitions)
        )
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for count, result in enumerate(
            pool.map(position_errors, itertools.repeat(data), trials), start=1
        ):
            counts = result.pop("errors")
            rows.extend(
                {**result, "position": position, "errors": errors}
                for position, errors in enumerate(counts)
            )
            print(f"\r{count}/{len(trials)}", end="", flush=True)
    print()
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".csv"), index=False)
    args.output.with_suffix(".json").write_text(
        json.dumps({"file": args.file.name, **environment(data)}, indent=2, ensure_ascii=False)
    )
    print(f"resultados em {args.output}.csv")


if __name__ == "__main__":
    main()
