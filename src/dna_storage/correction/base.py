from abc import ABC, abstractmethod


class Corrector(ABC):
    name: str

    # quem herda escreve. toda condição recebe clusters, que é a lista de grupos de leituras, e
    # devolve uma sequência por grupo
    @abstractmethod
    def correct(self, clusters: list[list[str]]) -> list[str]: ...
