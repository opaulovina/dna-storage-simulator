import hashlib
import platform
import shutil
import subprocess
from dataclasses import dataclass
from importlib.metadata import version

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import SEQUENCE_LENGTH, encode
from dna_storage.correction.base import Corrector
from dna_storage.correction.mafft import MafftCorrector
from dna_storage.correction.reed_solomon import REDUNDANCY, RsHybridCorrector, add_parity
from dna_storage.metrics import mean_hamming, measure, recovery_rate

CONDITIONS = (MafftCorrector.name, RsHybridCorrector.name)


@dataclass(frozen=True)
# representa uma configuração de experimento, incluindo a taxa de erro, profundidade de leitura,
# semente aleatória, proporção I:D:S e redundância.
# imutável, como a ChannelConfig
class Trial:
    error_rate: float
    depth: int
    seed: int
    ids_ratio: tuple[float, float, float] = (1, 1, 1)
    redundancy: float = REDUNDANCY


def label(ids_ratio: tuple[float, float, float]) -> str:
    return ":".join(f"{part:g}" for part in ids_ratio)


# a fábrica do padrão Strategy. Recebe o nome da condição e devolve o objeto certo.
def corrector_for(condition: str, data_count: int, redundancy: float) -> Corrector:
    if condition == MafftCorrector.name:
        return MafftCorrector()
    if condition == RsHybridCorrector.name:
        return RsHybridCorrector(data_count, redundancy)
    raise ValueError(f"condição desconhecida: {condition}")


# 1. codifica o arquivo e acrescenta paridade (redundância)
# 2. gera as leituras a partir da semente
# 3. a condição A pega os 140 primeiros grupos e B os 161
# 4. mede (measure) e devolve uma linha com:
#    taxa de recuperação, Hamming, tempo, memória, moléculas e densidade
def run_condition(data: bytes, trial: Trial, condition: str) -> dict:
    sequences = encode(data)
    stored = add_parity(sequences, trial.redundancy)
    config = ChannelConfig(trial.error_rate, trial.depth, trial.seed, trial.ids_ratio)
    clusters = simulate_reads(stored, config)
    corrector = corrector_for(condition, len(sequences), trial.redundancy)
    molecules = len(stored) if isinstance(corrector, RsHybridCorrector) else len(sequences)
    result = measure(corrector, clusters[:molecules])
    return {
        "condition": corrector.name,
        "ids_ratio": label(trial.ids_ratio),
        "error_rate": trial.error_rate,
        "depth": trial.depth,
        "seed": trial.seed,
        "redundancy": trial.redundancy,
        "sequences": len(sequences),
        "molecules": molecules,
        "density": 8 * len(data) / (molecules * SEQUENCE_LENGTH),
        "recovery_rate": recovery_rate(result.recovered, sequences),
        "hamming": mean_hamming(result.recovered, sequences),
        "seconds": result.seconds,
        "memory_kib": result.memory_kib,
    }


# registra o SHA-256 do arquivo e as versões do MAFFT, Python, reedsolo e do próprio pacote
def environment(data: bytes) -> dict:
    mafft = subprocess.run(
        [shutil.which("mafft") or "mafft", "--version"], capture_output=True, text=True
    )
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "sequences": len(encode(data)),
        "dna_storage": version("dna-storage-simulator"),
        "reedsolo": version("reedsolo"),
        "mafft": mafft.stderr.strip(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
    }
