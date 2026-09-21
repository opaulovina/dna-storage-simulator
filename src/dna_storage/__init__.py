"""Simulador de correcao de erros em armazenamento de dados em DNA.

Compara duas condicoes de correcao sob o mesmo conjunto de leituras ruidosas:

- Condicao A: MAFFT (alinhamento multiplo) + consenso, redundancia fisica apenas.
- Condicao B: Reed-Solomon hibrido (reconstrucao por consenso + RS externo),
  redundancia fisica e logica.

Ver README.md para a arquitetura e o desenho experimental.
"""

__version__ = "0.1.0"
