import csv
import subprocess
import sys
from pathlib import Path

import pytest

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import encode
from dna_storage.correction.reed_solomon import add_parity
from dna_storage.experiment import (
    CONDITIONS,
    SITES,
    Trial,
    corrector_for,
    failure_counts,
    label,
    run_condition,
    series_name,
)

DATA = bytes(range(56))
SCRIPT = Path(__file__).parents[1] / "experiments" / "run_experiments.py"
SITES_SCRIPT = Path(__file__).parents[1] / "experiments" / "failure_sites.py"


def test_parity_at_the_end_keeps_the_data_reads_identical():
    sequences = encode(bytes(range(256)) * 4)
    config = ChannelConfig(error_rate=0.1, depth=5, seed=3)
    with_parity = simulate_reads(add_parity(sequences), config)
    without_parity = simulate_reads(sequences, config)
    assert with_parity[: len(sequences)] == without_parity


def test_label_writes_the_ratio_as_in_the_monograph():
    assert label((1, 1, 2)) == "1:1:2"
    assert label((1.0, 1.0, 1.0)) == "1:1:1"


def test_unknown_condition_is_rejected():
    with pytest.raises(ValueError, match="desconhecida"):
        corrector_for("RS puro", 2, 0.15)


@pytest.mark.parametrize("condition", CONDITIONS)
def test_clean_channel_recovers_everything(condition):
    row = run_condition(DATA, Trial(error_rate=0, depth=3, seed=0), condition)
    assert row["condition"] == condition
    assert row["recovery_rate"] == 1
    assert row["hamming"] == 0
    assert row["seconds"] >= 0
    assert row["memory_kib"] >= 0


def test_condition_b_pays_for_the_parity_molecules():
    trial = Trial(error_rate=0, depth=1, seed=0)
    mafft = run_condition(DATA, trial, CONDITIONS[0])
    hybrid = run_condition(DATA, trial, CONDITIONS[1])
    assert (mafft["sequences"], mafft["molecules"]) == (2, 2)
    assert (hybrid["sequences"], hybrid["molecules"]) == (2, 3)
    assert mafft["density"] == pytest.approx(8 * 56 / (2 * 120))
    assert hybrid["density"] == pytest.approx(8 * 56 / (3 * 120))


def test_same_trial_gives_the_same_result():
    trial = Trial(error_rate=0.15, depth=5, seed=7)
    first = run_condition(DATA, trial, CONDITIONS[1])
    second = run_condition(DATA, trial, CONDITIONS[1])
    assert first["recovery_rate"] == second["recovery_rate"]
    assert first["hamming"] == second["hamming"]


def test_batch_script_writes_one_row_per_run_and_resumes(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "out"
    command = [
        sys.executable, str(SCRIPT), str(source),
        "--error-rates", "0", "0.05",
        "--depths", "3",
        "--repetitions", "2",
        "--workers", "2",
        "--output", str(output),
    ]
    subprocess.run(command, check=True, capture_output=True)
    with open(output.with_suffix(".csv"), newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2 * 2 * len(CONDITIONS)
    assert output.with_suffix(".json").exists()

    rerun = subprocess.run(command, check=True, capture_output=True, text=True)
    assert "8 execuções já feitas; 0 a fazer" in rerun.stdout
    with open(output.with_suffix(".csv"), newline="") as handle:
        assert len(list(csv.DictReader(handle))) == len(rows)


def test_batch_script_refuses_to_mix_input_files(tmp_path):
    output = tmp_path / "out"
    for content in (DATA, DATA[::-1]):
        source = tmp_path / "input.bin"
        source.write_bytes(content)
        result = subprocess.run(
            [sys.executable, str(SCRIPT), str(source), "--error-rates", "0",
             "--depths", "1", "--repetitions", "1", "--output", str(output)],
            capture_output=True, text=True,
        )
    assert result.returncode != 0
    assert "outro arquivo" in result.stderr


def test_series_names_separate_reconstruction_and_redundancy():
    assert series_name("MAFFT", "-", 0.15) == "MAFFT"
    assert series_name("RS híbrido", "BMA", 0.15) == "BMA + RS 15%"
    assert series_name("RS híbrido", "NW-POA", 0.0) == "NW-POA (sem RS)"


def test_rows_record_the_reconstruction_of_condition_b_only():
    trial = Trial(error_rate=0, depth=1, seed=0, reconstruction="DBMA")
    assert run_condition(DATA, trial, CONDITIONS[0])["reconstruction"] == "-"
    assert run_condition(DATA, trial, CONDITIONS[1])["reconstruction"] == "DBMA"


def test_zero_redundancy_has_the_density_of_condition_a():
    trial = Trial(error_rate=0, depth=1, seed=0, redundancy=0)
    hybrid = run_condition(DATA, trial, CONDITIONS[1])
    assert hybrid["molecules"] == 2


def test_ablation_runs_only_condition_b_for_each_redundancy_and_reconstruction(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "ablation"
    subprocess.run(
        [
            sys.executable, str(SCRIPT), str(source),
            "--error-rates", "0.05",
            "--depths", "3",
            "--repetitions", "1",
            "--conditions", "RS híbrido",
            "--redundancy", "0", "0.15",
            "--reconstructions", "BMA", "NW-POA",
            "--output", str(output),
        ],
        check=True,
        capture_output=True,
    )
    with open(output.with_suffix(".csv"), newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {(row["redundancy"], row["reconstruction"]) for row in rows} == {
        ("0.0", "BMA"), ("0.0", "NW-POA"), ("0.15", "BMA"), ("0.15", "NW-POA")
    }


def test_batch_script_refuses_a_csv_with_old_columns(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "old"
    output.with_suffix(".csv").write_text("condition,ids_ratio\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(source), "--error-rates", "0",
         "--depths", "1", "--repetitions", "1", "--output", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "versão anterior" in result.stderr


@pytest.mark.parametrize("condition", CONDITIONS)
def test_failure_counts_add_up_to_the_failed_sequences(condition):
    trial = Trial(error_rate=0.15, depth=3, seed=1)
    row = failure_counts(DATA, trial, condition)
    failed = round((1 - row["recovery_rate"]) * row["sequences"])
    assert sum(row[site] for site in SITES) == failed


def test_failure_sites_script_runs_both_conditions_and_resumes(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "falhas"
    command = [
        sys.executable,
        str(SITES_SCRIPT),
        str(source),
        "--error-rates",
        "0.05",
        "--depths",
        "3",
        "--repetitions",
        "1",
        "--workers",
        "2",
        "--output",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True)
    with open(output.with_suffix(".csv"), newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {(row["condition"], row["redundancy"]) for row in rows} == {
        ("MAFFT", "0.0"),
        ("RS híbrido", "0.0"),
        ("RS híbrido", "0.15"),
    }
    again = subprocess.run(command, check=True, capture_output=True, text=True)
    assert "0 a fazer" in again.stdout
