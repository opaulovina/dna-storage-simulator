from collections import Counter

from dna_storage.codec.mapping import BASES

LOOKAHEAD = 4


# conta os símbolos e devolve o símbolo mais frequente na lista, ou None se a lista estiver vazia
def plurality(symbols: list[str]) -> str | None:
    if not symbols:
        return None
    return Counter(symbols).most_common(1)[0][0]


# conta quantas posições de read, a partir de start, batem com target
def matches(read: str, start: int, target: list[str | None]) -> int:
    return sum(
        1
        # percorre os elementos de target com seus índices
        for offset, base in enumerate(target)
        # verifica se a posição atual de read é igual à base de target, e se não ultrapassa o
        # tamanho de read
        if start + offset < len(read) and read[start + offset] == base
    )


# a tentativa de realinhar leitura acontece aqui.
def realign(read: str, pointer: int, base: str, ahead: list[str | None]) -> int:
    # em caso de substituição, a leitura continua, em pointer+1, segue.
    # em caso de deleção, a leitura continua, em pointer.
    # em caso de inserção, a leitura continua, em pointer+2, em pointer+1 deveria estar a base.
    substitution = matches(read, pointer + 1, ahead)
    deletion = matches(read, pointer, ahead)
    insertion = matches(read, pointer + 1, [base, *ahead[:-1]])
    best = max(substitution, deletion, insertion)
    if substitution == best:
        return pointer + 1
    if deletion == best:
        return pointer
    return pointer + 2


# aqui é a função que implementa o algoritmo BMA, que recebe uma lista de leituras, o tamanho da
# sequência original e um parâmetro lookahead (quantas posições à frente considerar para
# realinhamento).
def bma(reads: list[str], length: int, lookahead: int = LOOKAHEAD) -> str:
    # para cada uma das 120 posições:
    # active: as leituras que ainda têm base sob o ponteiro
    # base: o voto de pluralidade, se ninguém vota é A (0)
    # agreeing: quem concordou, avança 1
    # ahead: as próximas 4 bases, votadas só entre quem concordou, porque essas leituras estão
    # sincronizadas
    # quem discordou passa por realign
    pointers = [0] * len(reads)
    consensus = []
    for _ in range(length):
        active = [k for k, read in enumerate(reads) if pointers[k] < len(read)]
        base = plurality([reads[k][pointers[k]] for k in active]) or BASES[0]
        agreeing = {k for k in active if reads[k][pointers[k]] == base}
        for k in agreeing:
            pointers[k] += 1
        ahead = [
            plurality(
                [
                    reads[k][pointers[k] + offset]
                    for k in agreeing
                    if pointers[k] + offset < len(reads[k])
                ]
            )
            for offset in range(lookahead)
        ]
        for k in active:
            if k not in agreeing:
                pointers[k] = realign(reads[k], pointers[k], base, ahead)
        consensus.append(base)
    return "".join(consensus)
