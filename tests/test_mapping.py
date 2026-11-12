"""Testes da camada de conversão entre bits e bases."""

import random

import pytest

from dna_storage.codec.mapping import (
    BASES,
    BASES_PER_BYTE,
    bases_to_bytes,
    bases_to_int,
    bytes_to_bases,
    int_to_bases,
)

# --- a constante --------------------------------------------------------------


def test_alphabet_follows_the_pfc1_table():
    # 00→A, 01→T, 10→G, 11→C, lido como valores 0, 1, 2, 3.
    assert BASES == "ATGC"


def test_alphabet_has_four_distinct_bases():
    # Pega erro de digitação numa constante que nada mais no sistema sabe conferir.
    assert len(BASES) == 4
    assert len(set(BASES)) == 4


# --- payload: bytes para bases ------------------------------------------------


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        (b"", ""),
        (b"\x00", "AAAA"),  # 00 00 00 00
        (b"\xff", "CCCC"),  # 11 11 11 11
        (b"\x1b", "ATGC"),  # 00 01 10 11
        (b"\x1b\x00", "ATGCAAAA"),
    ],
)
def test_known_vectors(data, expected):
    assert bytes_to_bases(data) == expected
    assert bases_to_bytes(expected) == data


def test_length_is_always_four_bases_per_byte():
    data = bytes(range(256))
    assert len(bytes_to_bases(data)) == BASES_PER_BYTE * len(data)


def test_round_trip_over_every_possible_byte():
    # 256 valores: dá para cobrir o espaço inteiro de entrada, o que é raro.
    data = bytes(range(256))
    assert bases_to_bytes(bytes_to_bases(data)) == data


def test_round_trip_over_random_data():
    rng = random.Random(42)
    for size in (1, 7, 64, 1000):
        data = bytes(rng.randrange(256) for _ in range(size))
        assert bases_to_bytes(bytes_to_bases(data)) == data


def test_output_uses_only_the_alphabet():
    assert set(bytes_to_bases(bytes(range(256)))) <= set(BASES)


# --- payload: erros -----------------------------------------------------------


def test_desynchronized_length_is_rejected():
    # Sintoma de deleção não corrigida: a sequência não fecha em bytes.
    with pytest.raises(ValueError, match="sincronizada"):
        bases_to_bytes("ATG")


def test_invalid_base_is_rejected():
    with pytest.raises(ValueError, match="inválida"):
        bases_to_bytes("ATGN")


def test_lowercase_is_rejected():
    # MAFFT pode devolver o alinhamento em minúsculas; a normalização é
    # responsabilidade de quem chama, não desta camada.
    with pytest.raises(ValueError, match="inválida"):
        bases_to_bytes("atgc")


# --- índice: inteiro para bases -----------------------------------------------


@pytest.mark.parametrize(
    ("value", "length", "expected"),
    [
        (0, 4, "AAAA"),
        (1, 4, "AAAT"),
        (27, 4, "ATGC"),  # 0*64 + 1*16 + 2*4 + 3
        (255, 4, "CCCC"),
        (5, 8, "AAAAAATT"),  # zeros à esquerda são bases A
    ],
)
def test_index_known_values(value, length, expected):
    assert int_to_bases(value, length) == expected
    assert bases_to_int(expected) == value


def test_index_has_fixed_width():
    # O comprimento fixo é o que permite fatiar índice e payload por posição.
    for value in (0, 1, 1000, 65535):
        assert len(int_to_bases(value, 8)) == 8


def test_index_round_trip():
    rng = random.Random(7)
    for length in (4, 8):
        capacity = 4**length
        values = [0, 1, capacity - 1] + [rng.randrange(capacity) for _ in range(50)]
        for value in values:
            assert bases_to_int(int_to_bases(value, length)) == value


def test_index_capacity_for_eight_bases():
    # 8 bases de índice endereçam 65.536 sequências (D-09).
    assert int_to_bases(65535, 8) == "CCCCCCCC"
    with pytest.raises(ValueError, match="não cabe"):
        int_to_bases(65536, 8)


def test_index_negative_is_rejected():
    with pytest.raises(ValueError, match="negativo"):
        int_to_bases(-1, 4)


def test_index_invalid_base_is_rejected():
    with pytest.raises(ValueError, match="inválida"):
        bases_to_int("ATGN")


# --- a invariante que sustenta a tabela única ---------------------------------


def test_both_readings_of_the_table_agree():
    """O par de bits de uma base é a representação binária do seu dígito.

    Este teste tranca a propriedade que justifica usar uma constante só para o
    payload e para o índice. Se alguém reordenar BASES, ele quebra.
    """
    for byte in range(256):
        assert bytes_to_bases(bytes([byte])) == int_to_bases(byte, BASES_PER_BYTE)
