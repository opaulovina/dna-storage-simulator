from abc import ABC, abstractmethod


class Corrector(ABC):
    name: str

    @abstractmethod
    def correct(self, clusters: list[list[str]]) -> list[str]: ...
