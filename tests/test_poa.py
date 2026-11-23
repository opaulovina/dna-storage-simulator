from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import SEQUENCE_LENGTH, encode
from dna_storage.correction.poa import fixed_length_consensus, nw_poa

SEQUENCE = encode(bytes(range(28)))[0]


def test_identical_reads_give_the_sequence_back():
    assert nw_poa([SEQUENCE] * 3, SEQUENCE_LENGTH) == SEQUENCE


def test_reads_of_the_monograph_example_recover_the_original():
    assert nw_poa(["ATCG", "ATCG", "ATCG", "ATG", "ATACG"], 4) == "ATCG"


def test_one_error_of_each_type_is_corrected():
    other = "A" if SEQUENCE[50] != "A" else "T"
    substitution = SEQUENCE[:50] + other + SEQUENCE[51:]
    deletion = SEQUENCE[:30] + SEQUENCE[31:]
    insertion = SEQUENCE[:90] + "G" + SEQUENCE[90:]
    reads = [SEQUENCE, SEQUENCE, substitution, deletion, insertion]
    assert nw_poa(reads, SEQUENCE_LENGTH) == SEQUENCE


def test_extra_columns_with_most_gaps_are_dropped():
    longer = SEQUENCE[:60] + "G" + SEQUENCE[60:]
    assert nw_poa([SEQUENCE, SEQUENCE, longer], SEQUENCE_LENGTH) == SEQUENCE


def test_without_gaps_the_leftmost_extra_column_is_dropped():
    longer = SEQUENCE[:60] + "G" + SEQUENCE[60:]
    assert nw_poa([longer, longer], SEQUENCE_LENGTH) == longer[1:]


def test_output_always_has_the_requested_length():
    assert len(nw_poa([SEQUENCE[:100], SEQUENCE[:90]], SEQUENCE_LENGTH)) == SEQUENCE_LENGTH
    assert nw_poa(["", ""], 5) == "AAAAA"


def test_low_error_rate_is_fully_recovered():
    sequences = encode(bytes(range(256)))
    clusters = simulate_reads(sequences, ChannelConfig(error_rate=0.03, depth=30, seed=1))
    assert [nw_poa(reads, SEQUENCE_LENGTH) for reads in clusters] == sequences


def test_fixed_length_consensus_drops_gap_columns_and_pads_short_alignments():
    alignment = ["AC-GT", "AC-GT", "ACTGT"]
    assert fixed_length_consensus(alignment, 4) == "ACGT"
    assert fixed_length_consensus(["AC", "AC"], 4) == "ACAA"
