"""Condicoes de correcao de erro, intercambiaveis atras de uma interface comum.

- ``base``: a interface ``Corrector``, que toda condicao implementa.
- ``mafft``: Condicao A, alinhamento multiplo pelo MAFFT seguido de consenso.
- ``bma``: reconstrucao de traco por Bitwise Majority Alignment com janela de antecipacao.
- ``reed_solomon``: Condicao B, BMA seguido de Reed-Solomon externo sobre a matriz de moleculas.
"""
