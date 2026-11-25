import random
from collections import Counter

import numpy as np
from reedsolo import ReedSolomonError

from dna_storage.channel import ChannelConfig, corrupt, simulate_reads
from dna_storage.codec.mapping import bases_to_bytes
from dna_storage.codec.payload import INDEX_LENGTH, SEQUENCE_LENGTH, encode
from dna_storage.correction.reed_solomon import (
    RECONSTRUCTIONS,
    ROWS,
    SYMBOL_BITS,
    RsHybridCorrector,
    add_parity,
    codec,
    to_symbols,
)
from dna_storage.experiment import Trial, label
from dna_storage.metrics import recovery_rate

ORGANICK_ERROR_RATE = 0.0064
ORGANICK_IDS_RATIO = (1.0, 5.0, 10.0)
COVERAGE_MODES = ("fixa", "poisson")
OUTCOMES = ("certo", "falha", "errado")


def poisson_reads(stored: list[str], trial: Trial) -> list[list[str]]:
    depths = np.random.default_rng(trial.seed).poisson(trial.depth, len(stored))
    rng = random.Random(trial.seed)
    return [
        [corrupt(sequence, trial.error_rate, trial.ids_ratio, rng) for _ in range(depth)]
        for sequence, depth in zip(stored, depths, strict=True)
    ]


def used_resilience(
    corrector: RsHybridCorrector, columns: dict[int, list[int]], truth: list[list[int]]
) -> list[float]:
    ratios = []
    for data, extra in corrector.units:
        if not extra:
            continue
        order = (*data, *extra)
        erased = sum(address not in columns for address in order)
        for row in range(ROWS):
            errors = sum(
                columns[address][row] != truth[address][row]
                for address in order
                if address in columns
            )
            ratios.append((2 * errors + erased) / len(extra))
    return ratios


def coverage_run(data: bytes, trial: Trial, mode: str) -> dict:
    if mode not in COVERAGE_MODES:
        raise ValueError(f"modo de cobertura desconhecido: {mode}")
    if trial.redundancy <= 0:
        raise ValueError("a used error resilience precisa de redundância")
    sequences = encode(data)
    stored = add_parity(sequences, trial.redundancy)
    if mode == "fixa":
        config = ChannelConfig(trial.error_rate, trial.depth, trial.seed, trial.ids_ratio)
        clusters = simulate_reads(stored, config)
    else:
        clusters = poisson_reads(stored, trial)
    corrector = RsHybridCorrector(len(sequences), trial.redundancy, trial.reconstruction)
    columns = corrector.receive(clusters)
    truth = [to_symbols(bases_to_bytes(sequence[INDEX_LENGTH:])) for sequence in stored]
    ratios = used_resilience(corrector, columns, truth)
    return {
        "mode": mode,
        "ids_ratio": label(trial.ids_ratio),
        "error_rate": trial.error_rate,
        "coverage": trial.depth,
        "seed": trial.seed,
        "redundancy": trial.redundancy,
        "reconstruction": trial.reconstruction,
        "molecules": len(stored),
        "dropouts": sum(not reads for reads in clusters),
        "erasures": len(stored) - len(columns),
        "max_resilience": round(max(ratios), 4),
        "mean_resilience": round(sum(ratios) / len(ratios), 4),
        "recovery_rate": recovery_rate(corrector.decode(columns), sequences),
    }


def frontier_cell(
    data_count: int, parity: int, errors: int, erasures: int, trials: int, rng: random.Random
) -> dict:
    length = data_count + parity
    outcomes = Counter()
    for _ in range(trials):
        message = [rng.randrange(2**SYMBOL_BITS) for _ in range(data_count)]
        word = list(codec(parity).encode(message))
        positions = rng.sample(range(length), errors + erasures)
        erased = positions[:erasures]
        for position in positions[erasures:]:
            word[position] = (word[position] + rng.randrange(1, 2**SYMBOL_BITS)) % 2**SYMBOL_BITS
        for position in erased:
            word[position] = 0
        try:
            decoded = list(codec(parity).decode(word, erase_pos=erased or None)[0])
        except ReedSolomonError:
            outcomes["falha"] += 1
            continue
        outcomes["certo" if decoded == message else "errado"] += 1
    return {
        "errors": errors,
        "erasures": erasures,
        "budget": 2 * errors + erasures,
        "parity": parity,
        "trials": trials,
        **{outcome: outcomes[outcome] for outcome in OUTCOMES},
    }


def position_errors(data: bytes, trial: Trial) -> dict:
    sequences = encode(data)
    config = ChannelConfig(trial.error_rate, trial.depth, trial.seed, trial.ids_ratio)
    reconstruct = RECONSTRUCTIONS[trial.reconstruction]
    counts = [0] * SEQUENCE_LENGTH
    for reads, original in zip(simulate_reads(sequences, config), sequences, strict=True):
        strand = reconstruct(reads, SEQUENCE_LENGTH)
        for position, (new, old) in enumerate(zip(strand, original, strict=True)):
            counts[position] += new != old
    return {
        "ids_ratio": label(trial.ids_ratio),
        "error_rate": trial.error_rate,
        "depth": trial.depth,
        "seed": trial.seed,
        "reconstruction": trial.reconstruction,
        "sequences": len(sequences),
        "errors": counts,
    }
