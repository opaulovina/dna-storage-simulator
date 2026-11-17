import pytest

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import encode
from dna_storage.correction.base import Corrector
from dna_storage.correction.mafft import MafftCorrector, consensus, parse_fasta, to_fasta

SEQUENCE = encode(bytes(range(28)))[0]


def test_consensus_takes_the_most_frequent_symbol_of_each_column():
    assert consensus(["ATCG", "ATCG", "AGCG"]) == "ATCG"


def test_consensus_drops_a_column_where_the_gap_wins():
    assert consensus(["AT-CG", "AT-CG", "ATACG"]) == "ATCG"


def test_consensus_keeps_the_base_where_the_gap_loses():
    assert consensus(["ATCG", "ATCG", "AT-G"]) == "ATCG"


def test_consensus_tie_goes_to_the_first_read():
    assert consensus(["A", "T"]) == "A"
    assert consensus(["T", "A"]) == "T"


def test_consensus_of_nothing_is_empty():
    assert consensus([]) == ""


def test_consensus_rejects_rows_of_different_lengths():
    with pytest.raises(ValueError):
        consensus(["ATCG", "ATC"])


def test_fasta_round_trip_returns_uppercase():
    assert parse_fasta(to_fasta(["atcg", "ATG"])) == ["ATCG", "ATG"]


def test_parse_fasta_joins_wrapped_lines():
    assert parse_fasta(">0\nATCG\nATCG\n>1\nAT\n") == ["ATCGATCG", "AT"]


def test_mafft_corrector_is_a_corrector():
    assert isinstance(MafftCorrector(), Corrector)


def test_missing_executable_is_reported():
    with pytest.raises(RuntimeError, match="não encontrado"):
        MafftCorrector("mafft-inexistente")


def test_alignment_rows_have_equal_length_and_valid_symbols():
    alignment = MafftCorrector().align(["ATCG", "ATG", "ATACG"])
    assert len(alignment) == 3
    assert len({len(row) for row in alignment}) == 1
    assert all(set(row) <= set("ATGC-") for row in alignment)


def test_reads_of_the_monograph_example_recover_the_original():
    reads = ["ATCG", "ATCG", "ATCG", "ATG", "ATACG"]
    assert MafftCorrector().correct([reads]) == ["ATCG"]


def test_one_error_of_each_type_is_corrected():
    other = "A" if SEQUENCE[50] != "A" else "T"
    substitution = SEQUENCE[:50] + other + SEQUENCE[51:]
    deletion = SEQUENCE[:30] + SEQUENCE[31:]
    insertion = SEQUENCE[:90] + "G" + SEQUENCE[90:]
    reads = [SEQUENCE, SEQUENCE, substitution, deletion, insertion]
    assert MafftCorrector().correct([reads]) == [SEQUENCE]


def test_one_recovered_sequence_per_cluster_in_the_same_order():
    sequences = encode(bytes(range(84)))
    clusters = [[sequence] * 3 for sequence in sequences]
    assert MafftCorrector().correct(clusters) == sequences


def test_a_single_read_is_returned_as_is():
    assert MafftCorrector().correct([[SEQUENCE]]) == [SEQUENCE]


def test_low_error_rate_is_fully_recovered():
    sequences = encode(bytes(range(256)))
    clusters = simulate_reads(sequences, ChannelConfig(error_rate=0.03, depth=30, seed=1))
    assert any(read not in sequences for reads in clusters for read in reads)
    assert MafftCorrector().correct(clusters) == sequences
