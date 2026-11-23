from spoa import poa

from dna_storage.codec.mapping import BASES
from dna_storage.correction.bma import plurality

GAP = "-"
GLOBAL_ALIGNMENT = 1


# alinha as leituras com o NW global sobre grafo de ordem parcial
# se o alinhamento passar de 120 colunas, tira as que têm mais lacunas
# cada coluna vota só entre bases
# se o alinhamento tiver menos de 120 colunas, completa com A's à direita
def nw_poa(reads: list[str], length: int) -> str:
    _, alignment = poa(reads, algorithm=GLOBAL_ALIGNMENT)
    return fixed_length_consensus(alignment, length)


def fixed_length_consensus(alignment: list[str], length: int) -> str:
    columns = list(zip(*alignment, strict=True))
    extra = len(columns) - length
    if extra > 0:
        by_gaps = sorted(range(len(columns)), key=lambda i: columns[i].count(GAP), reverse=True)
        dropped = set(by_gaps[:extra])
        columns = [column for i, column in enumerate(columns) if i not in dropped]
    consensus = "".join(
        plurality([symbol for symbol in column if symbol != GAP]) or BASES[0]
        for column in columns
    )
    return consensus.ljust(length, BASES[0])
