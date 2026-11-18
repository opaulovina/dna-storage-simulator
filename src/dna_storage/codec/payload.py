import math

from dna_storage.codec.mapping import (
    BASES_PER_BYTE,
    bases_to_bytes,
    bases_to_int,
    bytes_to_bases,
    int_to_bases,
)

SEQUENCE_LENGTH = 120
INDEX_LENGTH = 8
PAYLOAD_LENGTH = SEQUENCE_LENGTH - INDEX_LENGTH
BYTES_PER_SEQUENCE = PAYLOAD_LENGTH // BASES_PER_BYTE


# divide os dados em blocos de BYTES_PER_SEQUENCE bytes, adiciona zeros se necessário (A), e
# converte cada bloco em uma sequência de bases com índice
def encode(data: bytes) -> list[str]:
    sequences = []
    for start in range(0, len(data), BYTES_PER_SEQUENCE):
        chunk = data[start : start + BYTES_PER_SEQUENCE]
        missing = BYTES_PER_SEQUENCE - len(chunk)
        chunk = chunk + b"\x00" * missing
        index = start // BYTES_PER_SEQUENCE
        sequence = int_to_bases(index, INDEX_LENGTH) + bytes_to_bases(chunk)
        sequences.append(sequence)
    return sequences


# converte as sequências de bases em bytes, verifica se todas as sequências estão presentes e
# retorna os primeiros 'size' bytes
def decode(sequences: list[str], size: int) -> bytes:
    chunks = {}
    for sequence in sequences:
        if len(sequence) != SEQUENCE_LENGTH:
            raise ValueError(f"sequência com {len(sequence)} bases, esperado {SEQUENCE_LENGTH}")
        index = bases_to_int(sequence[:INDEX_LENGTH])
        chunks[index] = bases_to_bytes(sequence[INDEX_LENGTH:])

    count = math.ceil(size / BYTES_PER_SEQUENCE)
    data = bytearray()
    # verifica se todas as sequências estão presentes, caso contrário levanta um erro
    for index in range(count):
        if index not in chunks:
            raise ValueError(f"sequência de índice {index} ausente")
        data.extend(chunks[index])
    return bytes(data[:size])


# rebuild: remonta o arquivo para download na ordem das sequências.
# sequência curta é completada com A e sequência longa é cortada.
# o arquivo sai sempre do tamanho certo, com a corrupção visível (numa imagem aparecem os blocos
# estragados)
# não entra nas métricas (5.6).
def rebuild(sequences: list[str], size: int) -> bytes:
    data = bytearray()
    for sequence in sequences:
        payload = sequence[INDEX_LENGTH:SEQUENCE_LENGTH].ljust(PAYLOAD_LENGTH, "A")
        data.extend(bases_to_bytes(payload))
    return bytes(data[:size])
