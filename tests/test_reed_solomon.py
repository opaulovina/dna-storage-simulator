import pytest

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.mapping import BASES, bases_to_int, int_to_bases
from dna_storage.codec.payload import INDEX_LENGTH, SEQUENCE_LENGTH, encode
from dna_storage.correction.base import Corrector
from dna_storage.correction.reed_solomon import (
    CODEWORD_LIMIT,
    RsHybridCorrector,
    add_parity,
    layout,
    molecule_count,
    parity_count,
    to_payload,
    to_symbols,
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


def test_symbols_have_16_bits_and_one_long_codeword_covers_a_typical_file():
    assert CODEWORD_LIMIT == 65_535
    assert to_symbols(bytes([1, 2, 3, 4])) == [0x0102, 0x0304]
    assert to_payload([0x0102, 0x0304]) == bytes([1, 2, 3, 4])
    assert layout(586, 0.15) == [(range(0, 586), range(586, 674))]


def test_units_never_exceed_the_codeword_limit():
    for data, extra in layout(60_000, 0.15):
        assert len(data) + len(extra) <= CODEWORD_LIMIT
    assert sum(len(data) for data, _ in layout(60_000, 0.15)) == 60_000


def test_data_and_parity_must_fit_in_the_index():
    assert molecule_count(140, 0.15) == 161
    with pytest.raises(ValueError, match="65536"):
        add_parity(["A" * SEQUENCE_LENGTH] * 60_000, 0.15)


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


def readdress(sequence: str, address: int) -> str:
    return int_to_bases(address, INDEX_LENGTH) + sequence[INDEX_LENGTH:]


def test_an_address_out_of_range_becomes_an_erasure_and_is_recovered():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    clusters[5] = [substitute(stored[5], 0)]
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES


def test_a_repeated_address_discards_every_strand_that_claims_it():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    clusters[5] = [readdress(stored[5], 6)]
    corrector = RsHybridCorrector(140)
    reconstructed = [cluster[0] for cluster in clusters]
    assert 5 not in corrector.place(reconstructed)
    assert 6 not in corrector.place(reconstructed)
    assert corrector.correct(clusters) == SEQUENCES


def test_up_to_21_erasures_are_recovered():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    for column in range(21):
        clusters[column] = [readdress(stored[column], 1000 + column)]
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES
    clusters[21] = [readdress(stored[21], 2000)]
    assert RsHybridCorrector(140).correct(clusters) != SEQUENCES


def test_two_errors_cost_as_much_as_four_erasures():
    stored = add_parity(SEQUENCES)
    clusters = [[sequence] for sequence in stored]
    for column in range(10):
        clusters[column] = [substitute(stored[column], INDEX_LENGTH)]
    clusters[20] = [readdress(stored[20], 1000)]
    assert RsHybridCorrector(140).correct(clusters) == SEQUENCES
    clusters[21] = [readdress(stored[21], 1001)]
    assert RsHybridCorrector(140).correct(clusters) != SEQUENCES


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
