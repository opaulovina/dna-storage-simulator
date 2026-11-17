from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import SEQUENCE_LENGTH, encode
from dna_storage.correction.bma import bma, plurality, realign

SEQUENCE = encode(bytes(range(28)))[0]


def test_plurality_tie_goes_to_the_first_symbol():
    assert plurality(["A", "T"]) == "A"
    assert plurality(["T", "A", "A", "T"]) == "T"


def test_plurality_of_nothing_is_none():
    assert plurality([]) is None


def test_realign_hypotheses_move_the_pointer_as_in_the_monograph():
    ahead = ["G", "A", "T", "C"]
    assert realign("CGATC", 0, "T", ahead) == 1
    assert realign("GATC", 0, "T", ahead) == 0
    assert realign("CTGATC", 0, "T", ahead) == 2


def test_identical_reads_give_the_sequence_back():
    assert bma([SEQUENCE] * 3, SEQUENCE_LENGTH) == SEQUENCE


def test_reads_of_the_monograph_example_recover_the_original():
    assert bma(["ATCG", "ATCG", "ATCG", "ATG", "ATACG"], 4) == "ATCG"


def test_one_error_of_each_type_is_corrected():
    other = "A" if SEQUENCE[50] != "A" else "T"
    substitution = SEQUENCE[:50] + other + SEQUENCE[51:]
    deletion = SEQUENCE[:30] + SEQUENCE[31:]
    insertion = SEQUENCE[:90] + "G" + SEQUENCE[90:]
    reads = [SEQUENCE, SEQUENCE, substitution, deletion, insertion]
    assert bma(reads, SEQUENCE_LENGTH) == SEQUENCE


def test_output_always_has_the_requested_length():
    assert len(bma([SEQUENCE[:100], SEQUENCE[:90]], SEQUENCE_LENGTH)) == SEQUENCE_LENGTH
    assert len(bma([SEQUENCE + "ACGT" * 5], SEQUENCE_LENGTH)) == SEQUENCE_LENGTH
    assert bma(["", ""], 5) == "AAAAA"


def test_low_error_rate_is_fully_recovered():
    sequences = encode(bytes(range(256)))
    clusters = simulate_reads(sequences, ChannelConfig(error_rate=0.03, depth=30, seed=1))
    assert any(read not in sequences for reads in clusters for read in reads)
    assert [bma(reads, SEQUENCE_LENGTH) for reads in clusters] == sequences
