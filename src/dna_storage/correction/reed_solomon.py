import math

from reedsolo import ReedSolomonError, RSCodec

from dna_storage.codec.mapping import bases_to_bytes, bytes_to_bases, int_to_bases
from dna_storage.codec.payload import BYTES_PER_SEQUENCE, INDEX_LENGTH, SEQUENCE_LENGTH
from dna_storage.correction.base import Corrector
from dna_storage.correction.bma import bma

CODEWORD_LIMIT = 255
REDUNDANCY = 0.15


# calcula o número de moléculas de paridade necessárias para uma dada quantidade de dados e
# redundância
def parity_count(data_count: int, redundancy: float) -> int:
    # o teto de 15% x dados. O round(..., 9) evita que um erro de ponto flutuante transforme
    # 3,0000000001 em 4.
    return math.ceil(round(redundancy * data_count, 9))


# maior número de moléculas de dados por bloco sem que dados + paridade passem de 255 símbolos
def unit_size(redundancy: float) -> int:
    if redundancy < 0:
        raise ValueError(f"redundância negativa: {redundancy}")
    size = CODEWORD_LIMIT
    while size > 0 and size + parity_count(size, redundancy) > CODEWORD_LIMIT:
        size -= 1
    if size == 0:
        raise ValueError(f"redundância grande demais para palavras de 255 bytes: {redundancy}")
    return size


# divide as moléculas em blocos de no máximo 255 símbolos e numera a paridade depois dos dados
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


# adiciona moléculas de paridade às sequências de dados, usando o layout definido acima
def add_parity(sequences: list[str], redundancy: float = REDUNDANCY) -> list[str]:
    # para cada linha (byte) do bloco, codifica os bytes dos dados e guarda os bytes de paridade.
    # Cada coluna de paridade vira uma molécula com índice próprio.
    payloads = [bases_to_bytes(sequence[INDEX_LENGTH:]) for sequence in sequences]
    parity = []
    for data, extra in layout(len(sequences), redundancy):
        if not extra:
            continue
        codec = RSCodec(len(extra))
        rows = [
            codec.encode(bytes(payloads[column][row] for column in data))[-len(extra) :]
            for row in range(BYTES_PER_SEQUENCE)
        ]
        for position, index in enumerate(extra):
            payload = bytes(row[position] for row in rows)
            parity.append(int_to_bases(index, INDEX_LENGTH) + bytes_to_bases(payload))
    return sequences + parity


# 1. aplica o BMA em cada grupo
# 2. monta a matriz de bytes
# 3. decodifica linha por linha. se uma linha tem erros demais, ela fica como veio, e a métrica
# registra a falha
# 4. devolve as 140 sequências de dados, cada uma com um índice reconstruído pelo BMA e o payload
# corrigido.
class RsHybridCorrector(Corrector):
    name = "RS híbrido"

    # inicializa o corretor híbrido de Reed-Solomon, que combina a reconstrução pelo BMA com a
    # decodificação Reed-Solomon das linhas da matriz
    def __init__(self, data_count: int, redundancy: float = REDUNDANCY) -> None:
        self.data_count = data_count
        self.redundancy = redundancy
        self.units = layout(data_count, redundancy)

    def correct(self, clusters: list[list[str]]) -> list[str]:
        expected = self.data_count + sum(len(extra) for _, extra in self.units)
        if len(clusters) != expected:
            raise ValueError(f"{len(clusters)} grupos recebidos, esperado {expected}")
        strands = [bma(reads, SEQUENCE_LENGTH) for reads in clusters]
        payloads = [bytearray(bases_to_bytes(strand[INDEX_LENGTH:])) for strand in strands]
        for data, extra in self.units:
            if not extra:
                continue
            codec = RSCodec(len(extra))
            for row in range(BYTES_PER_SEQUENCE):
                word = bytearray(payloads[column][row] for column in (*data, *extra))
                try:
                    message = codec.decode(word)[0]
                except ReedSolomonError:
                    continue
                for column, value in zip(data, message, strict=True):
                    payloads[column][row] = value
        return [
            strands[column][:INDEX_LENGTH] + bytes_to_bases(payloads[column])
            for column in range(self.data_count)
        ]
