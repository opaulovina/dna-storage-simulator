"""Condicoes de correcao de erro, intercambiaveis atras de uma interface comum.

- ``base``: a interface ``Corrector``, que toda condicao implementa.
- ``mafft``: Condicao A, alinhamento multiplo pelo MAFFT seguido de consenso.
- ``bma``: reconstrucao de traco por Bitwise Majority Alignment (simples e dupla face).
- ``poa``: reconstrucao por alinhamento de ordem parcial com Needleman-Wunsch (NW-POA).
- ``reed_solomon``: Condicao B, BMA seguido de Reed-Solomon externo sobre a matriz de moleculas.
"""
