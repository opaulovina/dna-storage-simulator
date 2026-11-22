import subprocess
import time

import pytest

from dna_storage.codec.payload import encode
from dna_storage.correction.base import Corrector
from dna_storage.metrics import (
    PeakMemory,
    alignment,
    failure_site,
    hamming,
    mean_hamming,
    measure,
    recovery_rate,
)


def test_alignment_counts_differences_and_columns():
    assert alignment("ATCG", "ATCG") == (0, 4)
    assert alignment("ATGG", "ATCG") == (1, 4)
    assert alignment("ATG", "ATCG") == (1, 4)
    assert alignment("ATACG", "ATCG") == (1, 5)
    assert alignment("", "ATCG") == (4, 4)


def test_alignment_prefers_the_shortest_among_the_cheapest():
    assert alignment("AT", "TA") == (2, 2)


def test_hamming_matches_the_examples_of_section_2_2():
    assert hamming("ATGG", "ATCG") == 0.25
    assert hamming("ATG", "ATCG") == 0.25
    assert hamming("ATACG", "ATCG") == 0.2


def test_hamming_equals_plain_hamming_when_there_are_only_substitutions():
    original = encode(bytes(range(28)))[0]
    changed = "".join(
        "A" if position % 10 == 0 else base for position, base in enumerate(original)
    )
    plain = sum(new != old for new, old in zip(changed, original, strict=True))
    assert hamming(changed, original) == plain / len(original)


def test_hamming_stays_between_zero_and_one():
    assert hamming("AAAA", "TTTT") == 1.0
    assert hamming("", "") == 0.0


def test_recovery_rate_counts_identical_sequences_only():
    assert recovery_rate(["ATCG", "ATGG", "ATCG"], ["ATCG", "ATCG", "ATCG"]) == 2 / 3


def test_metrics_reject_an_empty_comparison():
    with pytest.raises(ValueError, match="nenhuma"):
        recovery_rate([], [])
    with pytest.raises(ValueError, match="nenhuma"):
        mean_hamming([], [])


def test_mean_hamming_averages_over_all_sequences():
    assert mean_hamming(["ATCG", "ATGG"], ["ATCG", "ATCG"]) == 0.125


class Echo(Corrector):
    name = "eco"

    def correct(self, clusters: list[list[str]]) -> list[str]:
        return [reads[0] for reads in clusters]


def test_measure_returns_the_correction_with_time_and_memory():
    result = measure(Echo(), [["ATCG"], ["GGCC"]])
    assert result.recovered == ["ATCG", "GGCC"]
    assert result.seconds >= 0
    assert result.memory_kib >= 0


def test_peak_memory_sees_an_allocation_made_inside_the_block():
    with PeakMemory(interval=0.001) as memory:
        block = bytearray(50 * 1024 * 1024)
        time.sleep(0.05)
        del block
    assert memory.kib >= 40 * 1024


def test_peak_memory_counts_external_programs():
    with PeakMemory(interval=0.001) as memory:
        subprocess.run(
            ["python3", "-c", "import time; block = bytearray(60 * 1024 * 1024); time.sleep(0.3)"],
            check=True,
        )
    assert memory.kib >= 40 * 1024


def test_failure_site_separates_index_and_payload():
    original = encode(bytes(range(28)))[0]
    wrong_index = "C" + original[1:]
    wrong_payload = original[:-1] + ("A" if original[-1] != "A" else "C")
    assert failure_site(original, original) == "ok"
    assert failure_site(wrong_index, original) == "index"
    assert failure_site(wrong_payload, original) == "payload"
    assert failure_site("C" + wrong_payload[1:], original) == "both"
    assert failure_site(original[:100], original) == "payload"
