import json
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from dna_storage.analysis import (
    cheapest_depth,
    full_recovery_limit,
    load_runs,
    paired_difference,
)

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "experiments" / "analyze.py"


def run(condition, error_rate, depth, seed, recovery, redundancy=0.15, density=1.6):
    return {
        "condition": condition,
        "ids_ratio": "1:1:1",
        "error_rate": error_rate,
        "depth": depth,
        "seed": seed,
        "redundancy": redundancy,
        "density": density,
        "recovery_rate": recovery,
    }


def write(folder, name, rows, sha="abc", old=False):
    table = pd.DataFrame(rows)
    if not old:
        table["reconstruction"] = ["-" if c == "MAFFT" else "BMA" for c in table["condition"]]
    table.to_csv(folder / f"{name}.csv", index=False)
    (folder / f"{name}.json").write_text(json.dumps({"sha256": sha}))
    return folder / f"{name}.csv"


def test_load_runs_names_series_fills_old_files_and_drops_repeated_runs(tmp_path):
    old = write(
        tmp_path,
        "final",
        [run("MAFFT", 0.1, 5, 0, 0.5), run("RS híbrido", 0.1, 5, 0, 0.9)],
        old=True,
    )
    new = write(
        tmp_path,
        "ablacao",
        [run("RS híbrido", 0.1, 5, 0, 0.9), run("RS híbrido", 0.1, 5, 0, 0.8, 0.0)],
    )
    runs = load_runs([old, new])
    assert sorted(runs["series"]) == ["BMA (sem RS)", "BMA + RS 15%", "MAFFT"]
    assert runs.loc[runs["series"] == "MAFFT", "reconstruction"].item() == "-"
    assert runs.loc[runs["series"] == "BMA + RS 15%", "read_cost"].item() == pytest.approx(5 / 1.6)


def test_load_runs_refuses_different_input_files(tmp_path):
    first = write(tmp_path, "a", [run("MAFFT", 0.1, 5, 0, 1.0)], sha="x")
    second = write(tmp_path, "b", [run("MAFFT", 0.1, 5, 0, 1.0)], sha="y")
    with pytest.raises(ValueError, match="diferentes"):
        load_runs([first, second])


def test_limit_stops_at_the_first_rate_that_fails_in_any_seed(tmp_path):
    rows = [run("MAFFT", rate, 10, seed, 1.0) for rate in (0.02, 0.04, 0.08) for seed in (0, 1)]
    rows += [run("MAFFT", 0.06, 10, 0, 1.0), run("MAFFT", 0.06, 10, 1, 0.99)]
    rows += [run("MAFFT", 0.02, 5, 0, 0.9)]
    runs = load_runs([write(tmp_path, "r", rows)])
    limits = full_recovery_limit(runs).set_index("depth")["limit"]
    assert limits[10] == 0.04
    assert math.isnan(limits[5])


def test_paired_difference_compares_the_same_reads(tmp_path):
    rows = [run("MAFFT", 0.1, 5, seed, 0.5) for seed in (0, 1)]
    rows += [run("RS híbrido", 0.1, 5, 0, 0.75), run("RS híbrido", 0.1, 5, 1, 0.5)]
    runs = load_runs([write(tmp_path, "r", rows)])
    difference = paired_difference(runs, "BMA + RS 15%", "MAFFT")
    assert list(difference["difference"]) == [0.25, 0.0]


def test_cheapest_depth_needs_full_recovery_in_every_seed(tmp_path):
    rows = [run("MAFFT", 0.1, 5, 0, 1.0, density=2.0), run("MAFFT", 0.1, 5, 1, 0.9, density=2.0)]
    rows += [
        run("MAFFT", 0.1, depth, seed, 1.0, density=2.0) for depth in (10, 20) for seed in (0, 1)
    ]
    runs = load_runs([write(tmp_path, "r", rows)])
    cheapest = cheapest_depth(runs)
    assert cheapest[["depth", "read_cost"]].values.tolist() == [[10, 5.0]]


def test_analysis_script_writes_tables_and_figures_from_the_results(tmp_path):
    subprocess.run([sys.executable, str(SCRIPT), "--output", str(tmp_path)], cwd=ROOT, check=True)
    for name in ("limites.csv", "contra_mafft.csv", "ganho_rs.csv", "custo_minimo.csv"):
        assert (tmp_path / name).exists()
    assert (tmp_path / "recuperacao_111.pdf").exists()
    assert (tmp_path / "custo_leitura.png").exists()
