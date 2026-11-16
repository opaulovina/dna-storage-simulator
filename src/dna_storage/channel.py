import random
from dataclasses import dataclass

from dna_storage.codec.mapping import BASES

ERROR_TYPES = ("I", "D", "S")


@dataclass(frozen=True)
class ChannelConfig:
    error_rate: float
    depth: int
    seed: int
    ids_ratio: tuple[float, float, float] = (1, 1, 1)

    def __post_init__(self) -> None:
        if not 0 <= self.error_rate <= 1:
            raise ValueError(f"taxa de erro fora de [0, 1]: {self.error_rate}")
        if self.depth < 1:
            raise ValueError(f"profundidade precisa ser pelo menos 1: {self.depth}")
        if len(self.ids_ratio) != 3 or min(self.ids_ratio) < 0 or sum(self.ids_ratio) == 0:
            raise ValueError(f"proporção I:D:S inválida: {self.ids_ratio}")


def corrupt(
    sequence: str,
    error_rate: float,
    ids_ratio: tuple[float, float, float],
    rng: random.Random,
) -> str:
    read = []
    for base in sequence:
        if rng.random() >= error_rate:
            read.append(base)
            continue
        error = rng.choices(ERROR_TYPES, weights=ids_ratio)[0]
        if error == "S":
            read.append(rng.choice([other for other in BASES if other != base]))
        elif error == "I":
            read.append(rng.choice(BASES))
            read.append(base)
        elif error == "D":
            continue
    return "".join(read)


def simulate_reads(sequences: list[str], config: ChannelConfig) -> list[list[str]]:
    rng = random.Random(config.seed)
    reads = []
    for sequence in sequences:
        copies = []
        for _ in range(config.depth):
            copies.append(corrupt(sequence, config.error_rate, config.ids_ratio, rng))
        reads.append(copies)
    return reads
