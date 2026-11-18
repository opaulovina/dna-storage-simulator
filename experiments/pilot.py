import argparse
import csv
import itertools
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import mean

import psutil

from dna_storage.correction.mafft import MafftCorrector
from dna_storage.experiment import CONDITIONS, Trial, environment, run_condition

FINAL_ERROR_RATES = (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30)
FINAL_DEPTHS = (5, 10, 20, 30, 50, 100)
FINAL_RATIOS = 2
FINAL_REPETITIONS = 3
STARTUP_SAMPLES = 20


#abre o MAFFT 20 vezes com uma entrada minima ("A", "A") e tira a média
# isso é o custo de só abrir o programa
def mafft_startup(samples=STARTUP_SAMPLES):
    corrector = MafftCorrector()
    start = time.perf_counter()
    for _ in range(samples):
        corrector.align(["A", "A"])
    return (time.perf_counter() - start) / samples


# roda todas as execuções com "workers" processos ao mesmo tempo e devolve os resultados e o tempo
# total; o main chama duas vezes: primeiro isolado (workers=1) e depois em paralelo
def run_all(data, jobs, workers):
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=1) as pool:
        rows = list(
            pool.map(
                run_condition,
                itertools.repeat(data),
                [trial for trial, _ in jobs],
                [condition for _, condition in jobs],
            )
        )
    wall = time.perf_counter() - start
    return [{**row, "workers": workers} for row in rows], wall


def outcome(row):
    return row["condition"], row["recovery_rate"], row["hamming"]

# outcome e report comparam linha a linha se os dois modos deram o mesmo resultado
# tem que dar, porque o ruído vem da semente
#isso é um teste de reprodutibilidade.

def report(isolated, parallel, walls, workers, final_runs, startup):
    same = all(
        outcome(alone) == outcome(together)
        for alone, together in zip(isolated, parallel, strict=True)
    )
    print(f"\nresultados idênticos isolado x paralelo: {'sim' if same else 'NÃO'}")
    speedup = walls[1] / walls[workers]
    print(f"tempo total: isolado {walls[1]:.0f} s; {workers} em paralelo {walls[workers]:.0f} s")
    print(f"aceleração: {speedup:.1f}x")

    projection = {}
    # imprime, para cada condição: tempo isolado contra paralelo (a inflação), memória, a fração do
    # tempo gasta na partida do MAFFT e a projeção da grade final
    for condition in dict.fromkeys(row["condition"] for row in isolated):
        alone = mean(row["seconds"] for row in isolated if row["condition"] == condition)
        together = mean(row["seconds"] for row in parallel if row["condition"] == condition)
        memory = mean(row["memory_kib"] for row in isolated if row["condition"] == condition)
        inflation = together / alone - 1
        hours_isolated = final_runs * alone / 3600
        hours_parallel = hours_isolated / speedup
        projection[condition] = {
            "seconds_isolated": alone,
            "seconds_parallel": together,
            "inflation": inflation,
            "memory_kib": memory,
            "projected_hours_isolated": hours_isolated,
            "projected_hours_parallel": hours_parallel,
        }
        print(f"\n{condition}")
        print(f"  tempo por execução: isolado {alone:.2f} s; em paralelo {together:.2f} s")
        print(f"  inflação do tempo em paralelo: {inflation:+.0%}")
        print(f"  memória média: {memory / 1024:.1f} MiB")
        if condition == MafftCorrector.name:
            share = startup * isolated[0]["sequences"] / alone
            projection[condition]["startup_share"] = share
            print(f"  partida do MAFFT: {startup:.3f} s por grupo, {share:.0%} do tempo")
        print(f"  grade final ({final_runs} execuções):")
        print(f"    uma por vez: {hours_isolated:.1f} h")
        print(f"    {workers} em paralelo: {hours_parallel:.1f} h")
    return same, speedup, projection


def main():
    parser = argparse.ArgumentParser(description="Experimento piloto: custo real e paralelismo")
    parser.add_argument("file", type=Path)
    parser.add_argument("--workers", type=int, default=psutil.cpu_count(logical=False))
    parser.add_argument("--error-rates", type=float, nargs="+", default=[0.08, 0.16])
    parser.add_argument("--depths", type=int, nargs="+", default=[10, 50])
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--output", type=Path, default=Path("results/pilot"))
    args = parser.parse_args()

    data = args.file.read_bytes()
    if not data:
        raise SystemExit(f"{args.file} está vazio")
    jobs = [
        (Trial(error_rate, depth, seed), condition)
        for error_rate, depth, seed, condition in itertools.product(
            args.error_rates, args.depths, args.seeds, CONDITIONS
        )
    ]
    startup = mafft_startup()
    print(f"{len(jobs)} execuções, primeiro uma por vez, depois {args.workers} em paralelo")
    isolated, alone = run_all(data, jobs, 1)
    parallel, together = run_all(data, jobs, args.workers)
    walls = {1: alone, args.workers: together}

    final_runs = len(FINAL_ERROR_RATES) * len(FINAL_DEPTHS) * FINAL_RATIOS * FINAL_REPETITIONS
    same, speedup, projection = report(
        isolated, parallel, walls, args.workers, final_runs, startup
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output.with_suffix(".csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(isolated[0]))
        writer.writeheader()
        writer.writerows(isolated + parallel)
    metadata = {
        "file": args.file.name,
        **environment(data),
        "workers": args.workers,
        "wall_seconds": {str(key): value for key, value in walls.items()},
        "speedup": speedup,
        "mafft_startup_seconds": startup,
        "identical_results": same,
        "final_grid_runs_per_condition": final_runs,
        "projection": projection,
    }
    args.output.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False)
    )


if __name__ == "__main__":
    main()
