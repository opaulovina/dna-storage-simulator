BASES = "ATGC" #a posição é o valor da base, ou seja, A=0, T=1, G=2, C=3
BASE_TO_VALUE = {base: value for value, base in enumerate(BASES)} #caminho inverso, como um dicionário, para mapear a base para o valor correspondente
RADIX = len(BASES)
BASES_PER_BYTE = 4 #um byte são oito bits, e cada base é representada por dois bits, então 8/2 = 4 bases por byte


def bytes_to_bases(data: bytes) -> str: #recebe bytes e retorna a sequência de bases correspondente
    bases = []
    for byte in data:
        for shift in (6, 4, 2, 0): #para cada shift, value recebe o par de bits daquela posição, o >> empurra para a direita e o &0b11 (& 3) pega apenas os dois últimos bits, que correspondem a uma base
            value = (byte >> shift) & 0b11
            bases.append(BASES[value]) #adiciona a base correspondente à lista
    return "".join(bases) #cola tudo num único texto


def bases_to_bytes(sequence: str) -> bytes:
    if len(sequence) % BASES_PER_BYTE != 0: #um comprimento que não seja múltiplo de 4 não é válido, pois cada byte representa 4 bases
        raise ValueError(
            f"comprimento {len(sequence)} não é múltiplo de {BASES_PER_BYTE}: " #por que recusar comprimento errado: indel muda o comprimento e causa uma perda de sincronização (Seção 2.2). Esta camada roda depois da correção. Se chega um comprimento errado aqui, a correção falhou naquela sequência, e devolver bytes errados em silêncio contaminaria a métrica
            "a sequência não está sincronizada"
        )

    data = bytearray()
    for start in range(0, len(sequence), BASES_PER_BYTE): 
        byte = 0
        for base in sequence[start : start + BASES_PER_BYTE]: #para cada base na sequência, obtém o valor correspondente e atualiza o byte
            try:
                value = BASE_TO_VALUE[base] #tenta obter um valor correspondente à base, ou seja, A=0, T=1, G=2, C=3
            except KeyError:
                raise ValueError(f"base inválida na sequência: {base!r}") from None #erro de base inválida, ou seja, não corresponde a nenhuma das letras do dicionário
            byte = (byte << 2) | value #o espelho da ida, desloca o byte dois bits para a esquerda e adiciona o valor da base. (|) é o operador OR bit a bit, que combina os bits do byte com os bits do valor da base
        data.append(byte)
    return bytes(data)


def int_to_bases(value: int, length: int) -> str: #aqui o valor é convertido em uma sequência de bases, com o comprimento especificado. Se o valor for negativo ou não couber no comprimento especificado, uma exceção é lançada
    if value < 0:
        raise ValueError(f"valor negativo: {value}")

    capacity = RADIX**length
    if value >= capacity:
        raise ValueError(
            f"valor {value} não cabe em {length} bases (capacidade {capacity})"
        )

    digits = []
    for _ in range(length):
        value, remainder = divmod(value, RADIX)
        digits.append(BASES[remainder])
    return "".join(reversed(digits))


def bases_to_int(sequence: str) -> int:
    value = 0
    for base in sequence:
        try:
            value = value * RADIX + BASE_TO_VALUE[base]
        except KeyError:
            raise ValueError(f"base inválida na sequência: {base!r}") from None
    return value
