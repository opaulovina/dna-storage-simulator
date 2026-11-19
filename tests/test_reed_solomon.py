import pytest

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.mapping import BASES, bases_to_int
from dna_storage.codec.payload import INDEX_LENGTH, SEQUENCE_LENGTH, encode
from dna_storage.correction.base import Corrector
from dna_storage.correction.reed_solomon import (
    RsHybridCorrector,
    add_parity,
    layout,
    parity_count,
    unit_size,
)

DATA = bytes(range(256)) * 15 + bytes(range(80))
SEQUENCES = encode(DATA)


def substitute(sequence: str, position: int) -> str:
    other = next(base for base in BASES if base != sequence[position])
    return sequence[:position] + other + sequence[position + 1 :]


def test_xie_file_size_gives_140_molecules_and_21_parity_molecules():
    assert len(SEQUENCES) == 140
    assert parity_count(140, 0.15) == 21
    assert layout(140, 0.15) == [(range(0, 140), range(140, 161))]


def test_parity_count_ignores_floating_point_noise():
    assert parity_count(20, 0.15) == 3


def test_units_never_exceed_255_symbols():
    assert unit_size(0.15) == 221
    for data, extra in layout(1000, 0.15):
        assert len(data) + len(extra) <= 255
    assert sum(len(data) for data, _ in layout(1000, 0.15)) == 1000


def test_zero_redundancy_adds_nothing():
    assert add_parity(SEQUENCES, 0) == SEQUENCES


def test_negative_redundancy_is_rejected():
    with pytest.raises(ValueError, match="negativa"):
        unit_size(-0.1)


def test_parity_molecules_come_after_the_data_with_the_next_indices():
    stored = add_parity(SEQUENCES)
    assert stored[:140] == SEQUENCES
    assert len(stored) == 161
    assert all(len(sequence) == SEQUENCE_LENGTH for sequence in stored)
    assert [bases_to_int(sequence[:INDEX_LENGTH]) for sequence in stored[140:]] == list(
        range(140, 161)
    )


def test_rs_hybrid_is_a_corrector():
    assert isinstance(RsHybridCorrector(140), Corrector)


def test_clean_reads_return_the_data_molecules():
    clusters = [[sequence] for sequence in add_parity(SEQUENCES)]
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES


def test_wrong_number_of_clusters_is_rejected():
    with pytest.raises(ValueError, match="esperado 161"):
        RsHybridCorrector(140).correct([[sequence] for sequence in SEQUENCES])


def test_up_to_ten_wrong_bytes_per_row_are_corrected():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    for column in range(10):
        clusters[column] = [substitute(stored[column], INDEX_LENGTH)]
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES


def test_eleven_wrong_bytes_in_one_row_are_beyond_capacity():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    for column in range(11):
        clusters[column] = [substitute(stored[column], INDEX_LENGTH)]
    assert RsHybridCorrector(140).correct(clusters) != SEQUENCES


def test_the_index_is_not_protected_by_the_outer_code():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    clusters[5] = [substitute(stored[5], 0)]
    recovered = RsHybridCorrector(140).correct(clusters)
    assert recovered[5] != SEQUENCES[5]
    assert recovered[5][INDEX_LENGTH:] == SEQUENCES[5][INDEX_LENGTH:]


def test_moderate_noise_is_fully_recovered():
    stored = add_parity(SEQUENCES)
    clusters = simulate_reads(stored, ChannelConfig(error_rate=0.06, depth=10, seed=1))
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES


@pytest.mark.parametrize("reconstruction", ["BMA", "DBMA", "NW-POA"])
def test_every_reconstruction_recovers_moderate_noise(reconstruction):
    stored = add_parity(SEQUENCES)
    clusters = simulate_reads(stored, ChannelConfig(error_rate=0.02, depth=10, seed=1))
    assert RsHybridCorrector(140, reconstruction=reconstruction).correct(clusters) == SEQUENCES


def test_unknown_reconstruction_is_rejected():
    with pytest.raises(ValueError, match="desconhecida"):
        RsHybridCorrector(140, reconstruction="Trellis")
