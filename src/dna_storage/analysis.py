import json
from pathlib import Path

import pandas as pd

from dna_storage.experiment import series_name

POINT = ["series", "ids_ratio", "error_rate", "depth"]
TRIAL = ["ids_ratio", "error_rate", "depth", "seed"]
OLD_RECONSTRUCTION = {"MAFFT": "-", "RS híbrido": "BMA"}


# lê o SHA-256 do JSON do experimento
def checksum(path: Path) -> str:
    return json.loads(path.with_suffix(".json").read_text())["sha256"]


# para cada CSV:
# table preenche a reconstrucao que falta nos CSVs antigos
# runs junta tudo e cria series com o nome da configuracao
# read_cost profundidade dividida pela densidade: nucleotideos lidos por bit
# remove execuções repetidas mantendo a primeira
def load_runs(paths: list[Path]) -> pd.DataFrame:
    checksums = {checksum(path) for path in paths}
    if len(checksums) > 1:
        raise ValueError("os resultados vêm de arquivos de entrada diferentes")
    tables = []
    for path in paths:
        table = pd.read_csv(path)
        if "reconstruction" not in table:
            table["reconstruction"] = table["condition"].map(OLD_RECONSTRUCTION)
        tables.append(table)
    runs = pd.concat(tables, ignore_index=True)
    runs["series"] = [
        series_name(*values)
        for values in zip(
            runs["condition"], runs["reconstruction"], runs["redundancy"], strict=True
        )
    ]
    runs["read_cost"] = runs["depth"] / runs["density"]
    return runs.drop_duplicates(subset=[*POINT, "seed"], keep="first", ignore_index=True)


# para cada ponto - serie, proporcao, taxa, profundidade:
# calcula media, minimo e maximo entre as sementes
def by_point(runs: pd.DataFrame, metric: str = "recovery_rate") -> pd.DataFrame:
    return runs.groupby(POINT)[metric].agg(["mean", "min", "max", "count"]).reset_index()


# para cada serie, proporcao e profundidade:
def full_recovery_limit(runs: pd.DataFrame) -> pd.DataFrame:
    # a pior semente em cada taxa
    worst = runs.groupby(POINT)["recovery_rate"].min().reset_index()
    rows = []
    for (series, ids_ratio, depth), group in worst.groupby(["series", "ids_ratio", "depth"]):
        group = group.sort_values("error_rate")
        # verdadeiro ate a primeira taxa que falha
        complete = group["recovery_rate"].eq(1.0).cummin()
        # a maior taxa antes dessa falha
        limit = group.loc[complete, "error_rate"].max() if complete.any() else float("nan")
        rows.append({"series": series, "ids_ratio": ids_ratio, "depth": depth, "limit": limit})
    return pd.DataFrame(rows)


def paired_difference(runs: pd.DataFrame, series: str, baseline: str) -> pd.DataFrame:
    # uma coluna por serie para cada execucao
    # mesmas leituras
    wide = runs.pivot_table(index=TRIAL, columns="series", values="recovery_rate")
    # devolve a serie menos a de referencia
    return (wide[series] - wide[baseline]).rename("difference").dropna().reset_index()


def cheapest_depth(runs: pd.DataFrame) -> pd.DataFrame:
    # a pior semente em cada ponto
    worst = (
        runs.groupby(POINT)
        .agg(worst=("recovery_rate", "min"), read_cost=("read_cost", "first"))
        .reset_index()
    )
    # so os pontos com 100% em todas as sementes
    complete = worst[worst["worst"] == 1.0]
    # a menor profundidade desses pontos com o custo de leitura
    first = complete.groupby(["series", "ids_ratio", "error_rate"])["depth"].idxmin()
    columns = ["series", "ids_ratio", "error_rate", "depth", "read_cost"]
    return complete.loc[first, columns].reset_index(drop=True)
