import random

import pytest

from dna_storage.codec.payload import (
    BYTES_PER_SEQUENCE,
    INDEX_LENGTH,
    PAYLOAD_LENGTH,
    SEQUENCE_LENGTH,
    decode,
    encode,
)


def random_bytes(size: int) -> bytes:
    rng = random.Random(size)
    return bytes(rng.randrange(256) for _ in range(size))


def test_parameters_follow_xie_2023():
    assert SEQUENCE_LENGTH == 120
    assert INDEX_LENGTH == 8
    assert PAYLOAD_LENGTH == 112
    assert BYTES_PER_SEQUENCE == 28


def test_one_kilobyte_needs_37_sequences():
    assert len(encode(random_bytes(1024))) == 37


def test_140_sequences_hold_3920_bytes():
    assert len(encode(random_bytes(3920))) == 140


def test_every_sequence_has_120_bases():
    for sequence in encode(random_bytes(1000)):
        assert len(sequence) == SEQUENCE_LENGTH


def test_index_counts_from_zero():
    sequences = encode(random_bytes(100))
    assert sequences[0][:INDEX_LENGTH] == "AAAAAAAA"
    assert sequences[1][:INDEX_LENGTH] == "AAAAAAAT"
    assert sequences[3][:INDEX_LENGTH] == "AAAAAAAC"


@pytest.mark.parametrize("size", [0, 1, 27, 28, 29, 1024, 3920])
def test_round_trip(size):
    data = random_bytes(size)
    assert decode(encode(data), size) == data


def test_round_trip_with_shuffled_sequences():
    data = random_bytes(1024)
    sequences = encode(data)
    random.Random(1).shuffle(sequences)
    assert decode(sequences, len(data)) == data


def test_missing_sequence_is_rejected():
    data = random_bytes(1024)
    sequences = encode(data)
    del sequences[5]
    with pytest.raises(ValueError, match="ausente"):
        decode(sequences, len(data))


def test_sequence_with_wrong_length_is_rejected():
    data = random_bytes(1024)
    sequences = encode(data)
    sequences[5] = sequences[5][:-1]
    with pytest.raises(ValueError, match="119 bases"):
        decode(sequences, len(data))


def test_real_zero_bytes_at_the_end_are_kept():
    data = b"abc" + b"\x00" * 30
    assert decode(encode(data), len(data)) == data
