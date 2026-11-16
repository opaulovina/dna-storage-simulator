import random

import pytest

from dna_storage.channel import ChannelConfig, corrupt, simulate_reads
from dna_storage.codec.payload import encode

SEQUENCE = encode(bytes(range(28)))[0]


def test_one_group_of_reads_per_sequence_with_depth_copies():
    sequences = encode(bytes(range(100)))
    reads = simulate_reads(sequences, ChannelConfig(error_rate=0.1, depth=7, seed=1))
    assert len(reads) == len(sequences)
    assert all(len(copies) == 7 for copies in reads)


def test_zero_error_rate_copies_exactly():
    reads = simulate_reads([SEQUENCE], ChannelConfig(error_rate=0, depth=5, seed=1))
    assert reads == [[SEQUENCE] * 5]


def test_same_seed_gives_same_reads():
    config = ChannelConfig(error_rate=0.2, depth=10, seed=42)
    assert simulate_reads([SEQUENCE], config) == simulate_reads([SEQUENCE], config)


def test_different_seeds_give_different_reads():
    a = simulate_reads([SEQUENCE], ChannelConfig(error_rate=0.2, depth=10, seed=1))
    b = simulate_reads([SEQUENCE], ChannelConfig(error_rate=0.2, depth=10, seed=2))
    assert a != b


def test_substitution_always_changes_the_base():
    read = corrupt(SEQUENCE, 1.0, (0, 0, 1), random.Random(0))
    assert len(read) == len(SEQUENCE)
    assert all(new != old for new, old in zip(read, SEQUENCE, strict=True))


def test_deletion_at_full_rate_removes_every_base():
    assert corrupt(SEQUENCE, 1.0, (0, 1, 0), random.Random(0)) == ""


def test_insertion_keeps_the_original_base_right_after_it():
    read = corrupt(SEQUENCE, 1.0, (1, 0, 0), random.Random(0))
    assert len(read) == 2 * len(SEQUENCE)
    assert read[1::2] == SEQUENCE


def test_observed_substitution_rate_matches_nominal():
    rng = random.Random(3)
    reads = [corrupt(SEQUENCE, 0.1, (0, 0, 1), rng) for _ in range(1000)]
    changed = sum(new != old for read in reads for new, old in zip(read, SEQUENCE, strict=True))
    assert abs(changed / (1000 * len(SEQUENCE)) - 0.1) < 0.005


def test_observed_deletion_rate_matches_nominal():
    rng = random.Random(4)
    reads = [corrupt(SEQUENCE, 0.1, (0, 1, 0), rng) for _ in range(1000)]
    mean_length = sum(len(read) for read in reads) / 1000
    assert abs(mean_length - 108) < 0.6


def test_equal_insertions_and_deletions_keep_mean_length():
    rng = random.Random(5)
    reads = [corrupt(SEQUENCE, 0.2, (1, 1, 0), rng) for _ in range(2000)]
    mean_length = sum(len(read) for read in reads) / 2000
    assert abs(mean_length - 120) < 0.6


@pytest.mark.parametrize(
    "kwargs",
    [
        {"error_rate": -0.1, "depth": 1, "seed": 0},
        {"error_rate": 1.5, "depth": 1, "seed": 0},
        {"error_rate": 0.1, "depth": 0, "seed": 0},
        {"error_rate": 0.1, "depth": 1, "seed": 0, "ids_ratio": (0, 0, 0)},
        {"error_rate": 0.1, "depth": 1, "seed": 0, "ids_ratio": (1, -1, 1)},
    ],
)
def test_invalid_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        ChannelConfig(**kwargs)
