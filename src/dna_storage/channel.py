import random
from dataclasses import dataclass

from dna_storage.codec.mapping import BASES

ERROR_TYPES = ("I", "D", "S")


@dataclass(frozen=True) # campos obrigatórios criados e congelados
class ChannelConfig:
    error_rate: float # taxa por base
    depth: int # número de cópias
    seed: int #semente
    ids_ratio: tuple[float, float, float] = (1, 1, 1) # proporção ids 1:1:1

    #logo depois de criar a configuração, confere se ela faz sentido:
    def __post_init__(self) -> None:
        if not 0 <= self.error_rate <= 1: #taxa entre 0 e 1
            raise ValueError(f"taxa de erro fora de [0, 1]: {self.error_rate}")
        if self.depth < 1: #profundidade pelo menos 1
            raise ValueError(f"profundidade precisa ser pelo menos 1: {self.depth}")
        #proporção precisa ter três números, nenhum negativo e não todos zero.
        if len(self.ids_ratio) != 3 or min(self.ids_ratio) < 0 or sum(self.ids_ratio) == 0:
            raise ValueError(f"proporção I:D:S inválida: {self.ids_ratio}")


#aqui que geramos o ruido, com base na taxa de erro e na proporção de inserções, deleções e
#substituições
def corrupt(
    sequence: str,
    error_rate: float,
    ids_ratio: tuple[float, float, float],
    rng: random.Random,
) -> str:
    read = [] #a leitura que vai sendo montada
    for base in sequence:
        #sorteia um numero aleatório entre 0 e 1, se for maior ou igual à taxa de erro, a base é
        #mantida, não há erro nessa base
        if rng.random() >= error_rate:
            read.append(base) # então copie a base como está e passe para a próxima
            continue
        #se houve erro, 'error' recebe 'I', 'D' ou 'S' sorteado com os pesos da proporção
        error = rng.choices(ERROR_TYPES, weights=ids_ratio)[0]
        if error == "S":
            #se for substituição escreve uma base sorteada entre as que não é a base original
            read.append(rng.choice([other for other in BASES if other != base]))
        #se for inserção, sorteia uma base aleatória e escreve antes da base original
        elif error == "I":
            read.append(rng.choice(BASES))
            read.append(base)
        elif error == "D": #se for deleção, não escreve nada e passa para a próxima
            continue
    return "".join(read)


def simulate_reads(sequences: list[str], config: ChannelConfig) -> list[list[str]]:
    rng = random.Random(config.seed) # inicializa o gerador de números aleatórios com a semente
    #mesma semente, mesmos sorteios, mesmas leituras, sempre. é isso que torna a execução
    #reproduzível. Se não houvesse semente, cada execução geraria leituras diferentes, mesmo com a
    #mesma sequência e configuração.
    reads = []
    for sequence in sequences: #para cada sequence que a payload gerou.
        #gera depth leituras ruidosas e guarde no grupo
        copies = []
        for _ in range(config.depth):
            copies.append(corrupt(sequence, config.error_rate, config.ids_ratio, rng))
        reads.append(copies)
    return reads
