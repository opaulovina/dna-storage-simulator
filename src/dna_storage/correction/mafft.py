import shutil
import subprocess
from collections import Counter

from dna_storage.correction.base import Corrector

GAP = "-"
FFT_NS_2 = ("--retree", "2", "--maxiterate", "0")


# gera o formato FASTA a partir de uma lista de leituras
def to_fasta(reads: list[str]) -> str:
    return "".join(f">{number}\n{read}\n" for number, read in enumerate(reads))


# analisa o formato FASTA e retorna uma lista de leituras
def parse_fasta(text: str) -> list[str]:
    records = []
    for block in text.split(">")[1:]:
        lines = block.splitlines()
        records.append("".join(lines[1:]).upper())
    return records


# gera a sequência consenso a partir de um alinhamento múltiplo; a lacuna também vota, e a coluna
# em que ela vence sai do consenso
def consensus(alignment: list[str]) -> str:
    bases = []
    for column in zip(*alignment, strict=True):
        winner = Counter(column).most_common(1)[0][0]
        if winner != GAP:
            bases.append(winner)
    return "".join(bases)


# define a classe MafftCorrector, que herda de Corrector e implementa a correção de erros usando o
# MAFFT
class MafftCorrector(Corrector):
    name = "MAFFT"

    def __init__(self, executable: str = "mafft") -> None:
        path = shutil.which(executable)
        if path is None:
            raise RuntimeError(f"executável do MAFFT não encontrado: {executable}")
        self.executable = path

    def align(self, reads: list[str]) -> list[str]:
        result = subprocess.run(
            [self.executable, *FFT_NS_2, "--nuc", "--quiet", "-"],
            input=to_fasta(reads),
            capture_output=True,
            text=True,
            check=True,
        )
        return parse_fasta(result.stdout)

    def correct(self, clusters: list[list[str]]) -> list[str]:
        return [consensus(self.align(reads)) for reads in clusters]
