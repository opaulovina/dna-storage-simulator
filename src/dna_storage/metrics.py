import threading
import time
from dataclasses import dataclass

import psutil

from dna_storage.codec.payload import INDEX_LENGTH
from dna_storage.correction.base import Corrector

SAMPLING_INTERVAL = 0.01


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


# se as 8 primeiras bases diferem da original
# é o index
# se o resto difere da original - payload
# devolve: ok, index, payload ou both
def failure_site(recovered: str, original: str) -> str:
    index = recovered[:INDEX_LENGTH] != original[:INDEX_LENGTH]
    payload = recovered[INDEX_LENGTH:] != original[INDEX_LENGTH:]
    if index and payload:
        return "both"
    if index:
        return "index"
    if payload:
        return "payload"
    return "ok"


@dataclass(frozen=True)
class Measurement:
    recovered: list[str]
    seconds: float
    memory_kib: int


def resident_memory(process: psutil.Process) -> int:
    total = process.memory_info().rss
    for child in process.children(recursive=True):
        try:
            if child.cmdline() != process.cmdline():
                total += child.memory_info().rss
        except psutil.Error:
            continue
    return total


# um vigia de memória residente que mede o pico de memória do processo e filhos. O uso é com um
# gerenciador de contexto:
class PeakMemory:
    def __init__(self, interval: float = SAMPLING_INTERVAL) -> None:
        self.interval = interval
        self.process = psutil.Process()
        self.stop = threading.Event()
        self.baseline = 0
        self.peak = 0

    def watch(self) -> None:
        while not self.stop.wait(self.interval):
            self.peak = max(self.peak, resident_memory(self.process))

    def __enter__(self) -> "PeakMemory":
        self.baseline = self.peak = resident_memory(self.process)
        self.thread = threading.Thread(target=self.watch, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *exception) -> None:
        self.stop.set()
        self.thread.join()
        self.peak = max(self.peak, resident_memory(self.process))

    @property
    def kib(self) -> int:
        return (self.peak - self.baseline) // 1024


# mede o tempo de parede da correção com o MAFFT incluído;
def measure(corrector: Corrector, clusters: list[list[str]]) -> Measurement:
    # a memória que a correção acrescenta: o processo python e os programas do MAFFT, somados.
    with PeakMemory() as memory:
        start = time.perf_counter()
        recovered = corrector.correct(clusters)
        seconds = time.perf_counter() - start
    return Measurement(recovered, seconds, memory.kib)
