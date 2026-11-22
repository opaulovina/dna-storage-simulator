import hashlib
import platform
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import psutil

from dna_storage.channel import ChannelConfig, simulate_reads
from dna_storage.codec.payload import SEQUENCE_LENGTH, encode
from dna_storage.correction.base import Corrector
from dna_storage.correction.mafft import MafftCorrector
from dna_storage.correction.reed_solomon import REDUNDANCY, RsHybridCorrector, add_parity
from dna_storage.metrics import failure_site, mean_hamming, measure, recovery_rate

CONDITIONS = (MafftCorrector.name, RsHybridCorrector.name)
SITES = ("index", "payload", "both")
PACKAGES = (
    "dna-storage-simulator",
    "numpy",
    "pandas",
    "matplotlib",
    "psutil",
    "reedsolo",
    "pyspoa",
)


@dataclass(frozen=True)
# representa uma configuração de experimento, incluindo a taxa de erro, profundidade de leitura,
# semente aleatória, proporção I:D:S, redundância e reconstrução da Condição B.
# imutável, como a ChannelConfig
class Trial:
    error_rate: float
    depth: int
    seed: int
    ids_ratio: tuple[float, float, float] = (1, 1, 1)
    redundancy: float = REDUNDANCY
    reconstruction: str = "BMA"


def label(ids_ratio: tuple[float, float, float]) -> str:
    return ":".join(f"{part:g}" for part in ids_ratio)


def series_name(condition: str, reconstruction: str, redundancy: float) -> str:
    if condition == MafftCorrector.name:
        return condition
    if not redundancy:
        return f"{reconstruction} (sem RS)"
    return f"{reconstruction} + RS {redundancy:.0%}"


# a fábrica do padrão Strategy. Recebe o nome da condição e devolve o objeto certo.
def corrector_for(
    condition: str, data_count: int, redundancy: float, reconstruction: str = "BMA"
) -> Corrector:
    if condition == MafftCorrector.name:
        return MafftCorrector()
    if condition == RsHybridCorrector.name:
        return RsHybridCorrector(data_count, redundancy, reconstruction)
    raise ValueError(f"condição desconhecida: {condition}")


# 1. codifica o arquivo e acrescenta paridade (redundância)
# 2. gera as leituras a partir da semente
# 3. a condição A pega os 140 primeiros grupos e B os 161
# 4. mede (measure) e devolve uma linha com:
#    taxa de recuperação, Hamming, tempo, memória, moléculas e densidade
#    e também as sequências recuperadas (usadas pelo Streamlit)
def evaluate(data: bytes, trial: Trial, condition: str) -> tuple[dict, list[str]]:
    sequences = encode(data)
    stored = add_parity(sequences, trial.redundancy)
    config = ChannelConfig(trial.error_rate, trial.depth, trial.seed, trial.ids_ratio)
    clusters = simulate_reads(stored, config)
    corrector = corrector_for(condition, len(sequences), trial.redundancy, trial.reconstruction)
    hybrid = isinstance(corrector, RsHybridCorrector)
    molecules = len(stored) if hybrid else len(sequences)
    result = measure(corrector, clusters[:molecules])
    return {
        "condition": corrector.name,
        "ids_ratio": label(trial.ids_ratio),
        "error_rate": trial.error_rate,
        "depth": trial.depth,
        "seed": trial.seed,
        "redundancy": trial.redundancy,
        "reconstruction": trial.reconstruction if hybrid else "-",
        "sequences": len(sequences),
        "molecules": molecules,
        "density": 8 * len(data) / (molecules * SEQUENCE_LENGTH),
        "recovery_rate": recovery_rate(result.recovered, sequences),
        "hamming": mean_hamming(result.recovered, sequences),
        "seconds": result.seconds,
        "memory_kib": result.memory_kib,
    }, result.recovered


# só a linha, para o modo experimental
def run_condition(data: bytes, trial: Trial, condition: str) -> dict:
    return evaluate(data, trial, condition)[0]


def failure_counts(data: bytes, trial: Trial, condition: str) -> dict: # roda evaluate e para cada sequencia recuperada:
# sites conta quantas cairam em cada categoria
# devolve a linha com a configuracao e as tres contagens
    row, recovered = evaluate(data, trial, condition)
    sites = Counter(
        failure_site(new, old) for new, old in zip(recovered, encode(data), strict=True)
    )
    fields = ("condition", "ids_ratio", "error_rate", "depth", "seed", "redundancy")
    return {
        **{field: row[field] for field in fields},
        "reconstruction": row["reconstruction"],
        "sequences": row["sequences"],
        "recovery_rate": row["recovery_rate"],
        **{site: sites[site] for site in SITES},
    }


def processor() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor()


# registra o SHA-256 do arquivo, as versões do MAFFT e do Python, o modelo do processador, a
# quantidade de nucleos fisicos e logicos, a memoria total e a versão de cada biblioteca
def environment(data: bytes) -> dict:
    mafft = subprocess.run(
        [shutil.which("mafft") or "mafft", "--version"], capture_output=True, text=True
    )
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "sequences": len(encode(data)),
        "mafft": mafft.stderr.strip(),
        "python": platform.python_version(),
        "packages": {package: version(package) for package in PACKAGES},
        "platform": platform.platform(),
        "processor": processor(),
        "physical_cores": psutil.cpu_count(logical=False),
        "logical_cores": psutil.cpu_count(),
        "memory_gib": round(psutil.virtual_memory().total / 2**30, 1),
    }
