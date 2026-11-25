import math
from functools import cache

from reedsolo import ReedSolomonError, RSCodec

from dna_storage.codec.mapping import (
    RADIX,
    bases_to_bytes,
    bases_to_int,
    bytes_to_bases,
    int_to_bases,
)
from dna_storage.codec.payload import BYTES_PER_SEQUENCE, INDEX_LENGTH, SEQUENCE_LENGTH
from dna_storage.correction.base import Corrector
from dna_storage.correction.bma import bma, double_sided_bma
from dna_storage.correction.poa import nw_poa

SYMBOL_BITS = 16
SYMBOL_BYTES = SYMBOL_BITS // 8
ROWS = BYTES_PER_SEQUENCE // SYMBOL_BYTES
CODEWORD_LIMIT = 2**SYMBOL_BITS - 1
ADDRESSES = RADIX**INDEX_LENGTH
REDUNDANCY = 0.15
RECONSTRUCTIONS = {"BMA": bma, "DBMA": double_sided_bma, "NW-POA": nw_poa}


# calcula o número de moléculas de paridade necessárias para uma dada quantidade de dados e
# redundância
def parity_count(data_count: int, redundancy: float) -> int:
    # o teto de 15% x dados. O round(..., 9) evita que um erro de ponto flutuante transforme
    # 3,0000000001 em 4.
    return math.ceil(round(redundancy * data_count, 9))


# maior número de moléculas de dados por bloco
def unit_size(redundancy: float) -> int:
    if redundancy < 0:
        raise ValueError(f"redundância negativa: {redundancy}")
    size = CODEWORD_LIMIT
    while size > 0 and size + parity_count(size, redundancy) > CODEWORD_LIMIT:
        size -= 1
    if size == 0:
        raise ValueError(f"redundância grande demais para a palavra-código: {redundancy}")
    return size


# divide as moléculas em blocos e numera a paridade depois dos dados
def layout(data_count: int, redundancy: float) -> list[tuple[range, range]]:
    size = unit_size(redundancy)
    units = []
    parity_start = data_count
    for start in range(0, data_count, size):
        data = range(start, min(start + size, data_count))
        count = parity_count(len(data), redundancy)
        units.append((data, range(parity_start, parity_start + count)))
        parity_start += count
    return units


def molecule_count(data_count: int, redundancy: float) -> int:
    return data_count + sum(len(extra) for _, extra in layout(data_count, redundancy))


@cache
def codec(parity: int) -> RSCodec:
    return RSCodec(parity, nsize=CODEWORD_LIMIT, c_exp=SYMBOL_BITS)


def to_symbols(payload: bytes) -> list[int]:
    return [
        int.from_bytes(payload[start : start + SYMBOL_BYTES])
        for start in range(0, len(payload), SYMBOL_BYTES)
    ]


def to_payload(symbols: list[int]) -> bytes:
    return b"".join(symbol.to_bytes(SYMBOL_BYTES) for symbol in symbols)


# adiciona moléculas de paridade às sequências de dados, usando o layout definido acima
def add_parity(sequences: list[str], redundancy: float = REDUNDANCY) -> list[str]:
    total = molecule_count(len(sequences), redundancy)
    if total > ADDRESSES:
        raise ValueError(
            f"{total} moléculas com a paridade; o índice de {INDEX_LENGTH} bases numera no "
            f"máximo {ADDRESSES}"
        )
    # para cada linha do bloco, codifica os bytes dos dados e guarda os bytes de paridade.
    # Cada coluna de paridade vira uma molécula com índice próprio.
    symbols = [to_symbols(bases_to_bytes(sequence[INDEX_LENGTH:])) for sequence in sequences]
    parity = []
    for data, extra in layout(len(sequences), redundancy):
        if not extra:
            continue
        rows = [
            codec(len(extra)).encode([symbols[column][row] for column in data])[-len(extra) :]
            for row in range(ROWS)
        ]
        for position, index in enumerate(extra):
            payload = to_payload([row[position] for row in rows])
            parity.append(int_to_bases(index, INDEX_LENGTH) + bytes_to_bases(payload))
    return sequences + parity


# 1. aplica o BMA em cada grupo
# 2. monta a matriz de bytes
# 3. decodifica linha por linha. se uma linha tem erros demais, ela fica como veio, e a métrica
# registra a falha
# 4. devolve as 140 sequências de dados,
class RsHybridCorrector(Corrector):
    name = "RS híbrido"

    # inicializa o corretor híbrido de Reed-Solomon, que combina a reconstrução pelo BMA com a
    # decodificação Reed-Solomon das linhas da matriz
    def __init__(
        self, data_count: int, redundancy: float = REDUNDANCY, reconstruction: str = "BMA"
    ) -> None:
        if reconstruction not in RECONSTRUCTIONS:
            raise ValueError(f"reconstrução desconhecida: {reconstruction}")
        self.data_count = data_count
        self.redundancy = redundancy
        self.reconstruction = reconstruction
        self.units = layout(data_count, redundancy)
        self.total = molecule_count(data_count, redundancy)

    def place(self, strands: list[str]) -> dict[int, list[int]]:
        claims = {}
        for strand in strands:
            address = bases_to_int(strand[:INDEX_LENGTH])
            if address < self.total:
                claims.setdefault(address, []).append(strand)
        return {
            address: to_symbols(bases_to_bytes(found[0][INDEX_LENGTH:]))
            for address, found in claims.items()
            if len(found) == 1
        }

    def receive(self, clusters: list[list[str]]) -> dict[int, list[int]]:
        if len(clusters) != self.total:
            raise ValueError(f"{len(clusters)} grupos recebidos, esperado {self.total}")
        reconstruct = RECONSTRUCTIONS[self.reconstruction]
        return self.place([reconstruct(reads, SEQUENCE_LENGTH) for reads in clusters if reads])

    def decode(self, columns: dict[int, list[int]]) -> list[str]:
        empty = [0] * ROWS
        matrix = [list(columns.get(address, empty)) for address in range(self.total)]
        for data, extra in self.units:
            order = (*data, *extra)
            erased = [position for position, address in enumerate(order) if address not in columns]
            if not extra or len(erased) > len(extra):
                continue
            for row in range(ROWS):
                word = [matrix[address][row] for address in order]
                try:
                    message = codec(len(extra)).decode(word, erase_pos=erased or None)[0]
                except ReedSolomonError:
                    continue
                for address, value in zip(data, message, strict=True):
                    matrix[address][row] = value
        return [
            int_to_bases(address, INDEX_LENGTH) + bytes_to_bases(to_payload(matrix[address]))
            for address in range(self.data_count)
        ]

    def correct(self, clusters: list[list[str]]) -> list[str]:
        return self.decode(self.receive(clusters))
