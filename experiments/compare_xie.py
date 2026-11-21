import argparse
from pathlib import Path

import pandas as pd

MEASURES = {
    "sequências idênticas": lambda table: table["recovery_rate"],
    "bases certas (1 - Hamming)": lambda table: 1 - table["hamming"],
    "erros corrigidos (1 - Hamming / taxa)": lambda table: (
        1 - table["hamming"] / table["error_rate"]
    ),
}


def main():
    parser = argparse.ArgumentParser(description="Compara a validação com as Figuras 1 e 2a de Xie")
    parser.add_argument("--validation", type=Path, default=Path("results/validation_xie.csv"))
    parser.add_argument("--figures", type=Path, default=Path("data/xie_figuras.csv"))
    parser.add_argument("--fig2-ratio", default="1:1:1")
    args = parser.parse_args()

    ours = pd.read_csv(args.validation)
    ours = ours.groupby(["ids_ratio", "error_rate", "depth"], as_index=False)[
        ["recovery_rate", "hamming"]
    ].mean()
    xie = pd.read_csv(args.figures, dtype={"ids_ratio": str}).fillna({"ids_ratio": ""})
    xie.loc[xie["figure"] == "2a", "ids_ratio"] = args.fig2_ratio
    merged = xie.merge(ours, on=["ids_ratio", "error_rate", "depth"], how="inner")
    for figure, group in merged.groupby("figure"):
        print(f"\nFigura {figure} ({len(group)} pontos em comum)")
        for name, measure in MEASURES.items():
            difference = measure(group) - group["recovery_accuracy"]
            rmse = (difference**2).mean() ** 0.5
            print(f"  {name:<40} RMSE {rmse:.3f}  maior diferença {difference.abs().max():.3f}")


if __name__ == "__main__":
    main()
