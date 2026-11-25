import csv
import random
import subprocess
import sys
from pathlib import Path

import pytest

from dna_storage.codec.payload import SEQUENCE_LENGTH, encode
from dna_storage.correction.reed_solomon import RsHybridCorrector, add_parity
from dna_storage.experiment import Trial
from dna_storage.validation import (
    ORGANICK_ERROR_RATE,
    ORGANICK_IDS_RATIO,
    coverage_run,
    frontier_cell,
    poisson_reads,
    position_errors,
    used_resilience,
)

DATA = bytes(range(256)) * 15 + bytes(range(80))
ROOT = Path(__file__).parents[1] / "experiments"


def test_clean_reads_use_none_of_the_resilience():
    for mode in ("fixa", "poisson"):
        row = coverage_run(DATA, Trial(0.0, 20, 0), mode)
        assert row["max_resilience"] == 0
        assert row["recovery_rate"] == 1


def test_an_error_costs_two_coordinates_and_an_erasure_one():
    stored = add_parity(encode(DATA))
    corrector = RsHybridCorrector(140)
    columns = corrector.receive([[sequence] for sequence in stored])
    truth = {address: list(symbols) for address, symbols in columns.items()}
    columns[3] = [symbol ^ 1 for symbol in columns[3]]
    del columns[4]
    ratios = used_resilience(corrector, columns, [truth[address] for address in range(161)])
    assert ratios == [3 / 21] * 14


def test_poisson_reads_are_reproducible_and_lose_molecules_at_low_coverage():
    stored = add_parity(encode(DATA))
    trial = Trial(ORGANICK_ERROR_RATE, 1, 0, ORGANICK_IDS_RATIO)
    reads = poisson_reads(stored, trial)
    assert reads == poisson_reads(stored, trial)
    assert any(not cluster for cluster in reads)


def test_unknown_mode_and_zero_redundancy_are_rejected():
    with pytest.raises(ValueError, match="desconhecido"):
        coverage_run(DATA, Trial(0.0, 5, 0), "gama")
    with pytest.raises(ValueError, match="redundância"):
        coverage_run(DATA, Trial(0.0, 5, 0, redundancy=0.0), "fixa")


@pytest.mark.parametrize(("errors", "erasures"), [(10, 1), (5, 11), (0, 21)])
def test_the_decoder_corrects_everything_inside_the_bound(errors, erasures):
    cell = frontier_cell(140, 21, errors, erasures, 3, random.Random(0))
    assert cell["certo"] == 3


@pytest.mark.parametrize(("errors", "erasures"), [(11, 0), (0, 22), (6, 10)])
def test_the_decoder_never_recovers_beyond_the_bound(errors, erasures):
    cell = frontier_cell(140, 21, errors, erasures, 3, random.Random(0))
    assert cell["certo"] == 0


def test_position_errors_count_one_value_per_position():
    row = position_errors(DATA, Trial(0.0, 3, 0, redundancy=0.0, reconstruction="DBMA"))
    assert row["errors"] == [0] * SEQUENCE_LENGTH


def test_organick_script_writes_coverage_and_frontier(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "organick"
    command = [
        sys.executable,
        str(ROOT / "validate_organick.py"),
        str(source),
        "--coverages",
        "3",
        "--repetitions",
        "1",
        "--trials",
        "1",
        "--output",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True)
    with open(tmp_path / "organick_cobertura.csv", newline="") as handle:
        assert {row["mode"] for row in csv.DictReader(handle)} == {"fixa", "poisson"}
    with open(tmp_path / "organick_fronteira.csv", newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 12 * 23


def test_positional_script_writes_one_row_per_position(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(DATA)
    output = tmp_path / "perfil"
    command = [
        sys.executable,
        str(ROOT / "positional_errors.py"),
        str(source),
        "--error-rates",
        "0.05",
        "--depths",
        "3",
        "--repetitions",
        "1",
        "--reconstructions",
        "BMA",
        "--workers",
        "1",
        "--output",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True)
    with open(output.with_suffix(".csv"), newline="") as handle:
        assert len(list(csv.DictReader(handle))) == SEQUENCE_LENGTH
