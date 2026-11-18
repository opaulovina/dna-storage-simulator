import resource
import time
from dataclasses import dataclass

from dna_storage.correction.base import Corrector


# Needleman-Wunsch com custo 1 para substituição, inserção e deleção.
def alignment(recovered: str, original: str) -> tuple[int, int]:
    # para cada base da recuperada e cada other da original, cada célula guarda o par (custo,
    # tamanho) e escolhe o menor entre três caminhos:
    # diagonal: casa as duas bases, custo +1 se forem diferentes;
    # de cima: lacuna numa sequência;
    # da esquerda: lacuna na outra.
    # o par compara primeiro o custo e depois o tamanho, assim, entre alinhamentos igualmente
    # baratos, fica o mais curto
    # devolve (diferenças, colunas) do alinhamento ótimo
    previous = [(column, column) for column in range(len(original) + 1)]
    for row, base in enumerate(recovered, start=1):
        current = [(row, row)]
        for column, other in enumerate(original, start=1):
            cost, size = previous[column - 1]
            diagonal = (cost + (base != other), size + 1)
            up = (previous[column][0] + 1, previous[column][1] + 1)
            left = (current[column - 1][0] + 1, current[column - 1][1] + 1)
            current.append(min(diagonal, up, left))
        previous = current
    return previous[-1]


def hamming(recovered: str, original: str) -> float: # diferenças / colunas. sempre entre 0 e 1.
    # com só substituições o resultado é igual à Hamming comum.
    differences, size = alignment(recovered, original)
    return differences / size if size else 0.0


def recovery_rate(recovered: list[str], original: list[str]) -> float:
    if not original:
        raise ValueError("nenhuma sequência para comparar")
    hits = sum(new == old for new, old in zip(recovered, original, strict=True))
    return hits / len(original)


#a média entre as sequências recuperadas e as originais, considerando a distância de Hamming. Se
#não houver sequências originais, levanta um ValueError.
def mean_hamming(recovered: list[str], original: list[str]) -> float:
    if not original:
        raise ValueError("nenhuma sequência para comparar")
    distances = [
        0.0 if new == old else hamming(new, old)
        for new, old in zip(recovered, original, strict=True)
    ]
    return sum(distances) / len(distances)


@dataclass(frozen=True)
class Measurement:
    recovered: list[str]
    seconds: float
    memory_kib: int


# mede o tempo de parede da correção com o MAFFT incluído;
def measure(corrector: Corrector, clusters: list[list[str]]) -> Measurement:
    # a memória que a correção acrescenta ao processo python; o pico dos processos externos.
    # só vale com cada execução num processo novo.
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    start = time.perf_counter()
    recovered = corrector.correct(clusters)
    seconds = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return Measurement(recovered, seconds, after - before + children)
