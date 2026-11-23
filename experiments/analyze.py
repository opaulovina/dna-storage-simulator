import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import NullLocator, PercentFormatter, ScalarFormatter

from dna_storage.analysis import (
    by_point,
    cheapest_depth,
    checksum,
    full_recovery_limit,
    load_runs,
    paired_difference,
)
from dna_storage.experiment import series_name

TEXT = ["results/final_texto.csv", "results/ablacao.csv", "results/reconstrucoes.csv"]
MAIN = ["MAFFT", "BMA (sem RS)", "BMA + RS 15%"]
LIMITS = ["MAFFT", "BMA (sem RS)", "BMA + RS 15%", "DBMA (sem RS)", "NW-POA (sem RS)"]
RECONSTRUCTIONS = ["BMA (sem RS)", "DBMA (sem RS)", "NW-POA (sem RS)"]
REDUNDANCIES = ["BMA + RS 5%", "BMA + RS 10%", "BMA + RS 15%", "BMA + RS 25%", "BMA + RS 50%"]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#eda100"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
STYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]
INK, MUTED = "#0b0b0b", "#52514e"


ENTITIES = [
    "MAFFT",
    "BMA (sem RS)",
    "BMA + RS 15%",
    "BMA + RS 50%",
    "DBMA (sem RS)",
    "NW-POA (sem RS)",
]


def style(series):
    position = ENTITIES.index(series)
    return {
        "color": COLORS[position],
        "marker": MARKERS[position],
        "linestyle": STYLES[position],
        "linewidth": 2,
        "markersize": 5,
    }


def depth_axis(axis, depths):
    axis.set_xscale("log")
    axis.set_xticks(depths, [str(depth) for depth in depths])
    axis.xaxis.set_minor_locator(NullLocator())


def finish(axis, xlabel, ylabel):
    axis.set_xlabel(xlabel, color=MUTED)
    axis.set_ylabel(ylabel, color=MUTED)
    axis.grid(alpha=0.25)
    axis.tick_params(colors=MUTED)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)


def save(figure, folder, name):
    figure.tight_layout()
    for suffix in (".pdf", ".png"):
        figure.savefig(folder / f"{name}{suffix}", dpi=200, bbox_inches="tight")
    plt.close(figure)


def recovery_figure(points, ratio, folder):
    depths = sorted(points["depth"].unique())
    figure, axes = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
    for axis, depth in zip(axes.flat, depths, strict=True):
        for series in MAIN:
            line = points[
                (points["series"] == series)
                & (points["ids_ratio"] == ratio)
                & (points["depth"] == depth)
            ]
            axis.plot(line["error_rate"], line["mean"], label=series, **style(series))
        axis.set_title(f"{depth} leituras", color=INK, fontsize=10)
        axis.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        finish(axis, "taxa de erro por base", "taxa de recuperação")
    for axis in axes[0]:
        axis.set_xlabel("")
    for axis in axes[:, 1:].flat:
        axis.set_ylabel("")
    axes[0, 0].legend(frameon=False, fontsize=8)
    save(figure, folder, f"recuperacao_{ratio.replace(':', '')}")


def limit_figure(limits, ratio, folder):
    figure, axis = plt.subplots(figsize=(6.5, 4.2))
    for series in LIMITS:
        line = limits[(limits["series"] == series) & (limits["ids_ratio"] == ratio)]
        line = line.sort_values("depth")
        axis.plot(line["depth"], line["limit"], label=series, **style(series))
    depth_axis(axis, sorted(limits["depth"].unique()))
    axis.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    axis.legend(frameon=False, fontsize=8)
    finish(axis, "leituras por sequência", "maior taxa com recuperação completa")
    save(figure, folder, f"limites_{ratio.replace(':', '')}")


def gain_figure(runs, ratio, folder):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for axis, series in zip(axes, ["BMA + RS 15%", "BMA + RS 50%"], strict=True):
        gain = paired_difference(runs, series, "BMA (sem RS)")
        gain = gain[gain["ids_ratio"] == ratio]
        grid = gain.pivot_table(index="error_rate", columns="depth", values="difference")
        image = axis.imshow(
            grid.values, cmap="Blues", vmin=0, vmax=0.2, aspect="auto", origin="lower"
        )
        axis.set_xticks(range(len(grid.columns)), grid.columns)
        axis.set_yticks(range(len(grid.index)), [f"{rate:.0%}" for rate in grid.index])
        for row, rate in enumerate(grid.index):
            for column, depth in enumerate(grid.columns):
                value = grid.loc[rate, depth]
                if value >= 0.005:
                    color = "white" if value > 0.12 else INK
                    label = f"{value:.2f}"
                    axis.text(column, row, label, ha="center", va="center", fontsize=7, color=color)
        axis.set_title(f"{series} contra BMA (sem RS)", color=INK, fontsize=10)
        axis.set_xlabel("leituras por sequência", color=MUTED)
    axes[0].set_ylabel("taxa de erro por base", color=MUTED)
    figure.colorbar(image, ax=axes, label="ganho na taxa de recuperação", shrink=0.8)
    for suffix in (".pdf", ".png"):
        figure.savefig(
            folder / f"ganho_rs_{ratio.replace(':', '')}{suffix}", dpi=200, bbox_inches="tight"
        )
    plt.close(figure)


def reconstruction_figure(points, ratio, folder, rates=(0.08, 0.12, 0.16, 0.20)):
    figure, axes = plt.subplots(1, len(rates), figsize=(12, 3.4), sharey=True)
    for axis, rate in zip(axes, rates, strict=True):
        for series in RECONSTRUCTIONS:
            line = points[
                (points["series"] == series)
                & (points["ids_ratio"] == ratio)
                & (points["error_rate"] == rate)
            ]
            axis.plot(line["depth"], line["mean"], label=series, **style(series))
        depth_axis(axis, sorted(points["depth"].unique()))
        axis.set_title(f"erro de {rate:.0%}", color=INK, fontsize=10)
        finish(axis, "leituras por sequência", "")
    axes[0].set_ylabel("taxa de recuperação", color=MUTED)
    axes[0].legend(frameon=False, fontsize=8)
    save(figure, folder, f"reconstrucoes_{ratio.replace(':', '')}")


def cost_figure(cheapest, folder):
    ratios = sorted(cheapest["ids_ratio"].unique())
    figure, axes = plt.subplots(1, len(ratios), figsize=(10, 3.8), sharey=True)
    for axis, ratio in zip(axes, ratios, strict=True):
        for series in MAIN:
            line = cheapest[(cheapest["series"] == series) & (cheapest["ids_ratio"] == ratio)]
            axis.plot(line["error_rate"], line["read_cost"], label=series, **style(series))
        axis.set_yscale("log")
        axis.set_yticks([3, 5, 10, 20, 50])
        axis.yaxis.set_major_formatter(ScalarFormatter())
        axis.yaxis.set_minor_locator(NullLocator())
        axis.set_xticks(sorted(cheapest["error_rate"].unique())[::2])
        axis.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        axis.set_title(f"I:D:S {ratio}", color=INK, fontsize=10)
        finish(axis, "taxa de erro por base", "")
    axes[0].set_ylabel("nucleotídeos lidos por bit", color=MUTED)
    axes[0].legend(frameon=False, fontsize=8)
    save(figure, folder, "custo_leitura")


def time_figure(runs, folder):
    final = runs[runs["series"].isin(["MAFFT", "BMA + RS 15%"])].copy()
    final["per_sequence"] = final["seconds"] / final["sequences"]
    means = final.groupby(["series", "depth"])["per_sequence"].mean().reset_index()
    figure, axis = plt.subplots(figsize=(6, 3.8))
    names = ["MAFFT", "BMA + RS 15%"]
    for series in names:
        line = means[means["series"] == series]
        axis.plot(line["depth"], line["per_sequence"], label=series, **style(series))
    axis.set_yscale("log")
    axis.set_yticks([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2])
    axis.yaxis.set_major_formatter(ScalarFormatter())
    axis.yaxis.set_minor_locator(NullLocator())
    depth_axis(axis, sorted(means["depth"].unique()))
    axis.legend(frameon=False, fontsize=8)
    finish(axis, "leituras por sequência", "segundos por sequência")
    save(figure, folder, "tempo")


def against_mafft(runs):
    rows = []
    for series in sorted(set(runs["series"]) - {"MAFFT"}):
        difference = paired_difference(runs, series, "MAFFT")["difference"]
        rows.append(
            {
                "series": series,
                "acima": int((difference > 0).sum()),
                "igual": int((difference == 0).sum()),
                "abaixo": int((difference < 0).sum()),
                "maior_perda": round(difference.min(), 3),
            }
        )
    return pd.DataFrame(rows)


def rs_gain(runs):
    rows = []
    for series in REDUNDANCIES:
        gain = paired_difference(runs, series, "BMA (sem RS)")
        mean = gain.groupby(["ids_ratio", "error_rate", "depth"])["difference"].mean()
        rows.append(
            {
                "series": series,
                "densidade": round(runs.loc[runs["series"] == series, "density"].iloc[0], 3),
                "ganho_medio": round(mean.mean(), 4),
                "maior_ganho": round(mean.max(), 3),
                "pontos_com_ganho": int((mean > 0).sum()),
                "pontos_com_perda": int((mean < 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def resources(runs):
    final = runs[runs["series"].isin(["MAFFT", "BMA + RS 15%"])]
    final = final.assign(per_sequence=final["seconds"] / final["sequences"])
    return final.groupby(["series", "depth"])[["per_sequence", "memory_kib"]].mean().round(4)


def image_against_text(image, points):
    point = ["series", "ids_ratio", "error_rate", "depth"]
    common = by_point(image).merge(points, on=point, suffixes=("_imagem", "_texto"))
    return common[[*point, "mean_imagem", "mean_texto"]].round(3)


def failure_summary(failures):
    sites = ["index", "payload", "both"]
    failures = failures.assign(failed=failures[sites].sum(axis=1))
    failures["series"] = [
        series_name(*values)
        for values in zip(
            failures["condition"], failures["reconstruction"], failures["redundancy"], strict=True
        )
    ]
    failures["only_index"] = (failures["failed"] > 0) & (failures["index"] == failures["failed"])
    groups = failures.groupby(["series", "ids_ratio"])
    summary = groups[[*sites, "failed"]].sum()
    summary["so_indice"] = (summary["index"] / summary["failed"]).round(3)
    summary["execucoes_com_falha"] = groups["failed"].apply(lambda failed: (failed > 0).sum())
    summary["execucoes_so_indice"] = groups["only_index"].sum()
    payload = failures.assign(recovery_rate=(failures["payload"] + failures["both"] == 0) * 1.0)
    return summary, full_recovery_limit(payload)


def mafft_length_summary(diagnostic, points):
    point = ["series", "ids_ratio", "error_rate", "depth"]
    fixed = diagnostic.assign(recovery_rate=diagnostic["recovery_rate_120"])
    runs = pd.concat(
        [diagnostic.assign(series="MAFFT"), fixed.assign(series="MAFFT (120 colunas)")]
    )
    means = by_point(runs).pivot_table(index=point[1:], columns="series", values="mean")
    bma = points[points["series"] == "BMA (sem RS)"].set_index(point[1:])["mean"]
    means["BMA (sem RS)"] = bma.reindex(means.index)
    return means.round(3), full_recovery_limit(runs)


def by_depth(limits):
    return limits.pivot_table(index=["ids_ratio", "series"], columns="depth", values="limit")


def main():
    parser = argparse.ArgumentParser(description="Tabelas e figuras da análise (etapa 11)")
    parser.add_argument("--text", type=Path, nargs="+", default=[Path(path) for path in TEXT])
    parser.add_argument("--image", type=Path, default=Path("results/final_imagem.csv"))
    parser.add_argument("--failures", type=Path, default=Path("results/falhas.csv"))
    parser.add_argument("--mafft-length", type=Path, default=Path("results/mafft_120.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/analise"))
    args = parser.parse_args()

    folder = args.output
    folder.mkdir(parents=True, exist_ok=True)
    runs = load_runs(args.text)
    points = by_point(runs)
    limits = full_recovery_limit(runs)
    cheapest = cheapest_depth(runs)
    points.to_csv(folder / "pontos.csv", index=False)
    by_depth(limits).to_csv(folder / "limites.csv")
    cheapest.to_csv(folder / "custo_minimo.csv", index=False)
    against_mafft(runs).to_csv(folder / "contra_mafft.csv", index=False)
    rs_gain(runs).to_csv(folder / "ganho_rs.csv", index=False)
    resources(runs).to_csv(folder / "tempo_memoria.csv")
    if args.image.exists():
        image = load_runs([args.image])
        image_against_text(image, points).to_csv(folder / "imagem_texto.csv", index=False)
    if args.failures.exists():
        if checksum(args.failures) != checksum(args.text[0]):
            raise SystemExit(f"{args.failures} foi gerado com outro arquivo de entrada")
        summary, payload_limits = failure_summary(pd.read_csv(args.failures))
        summary.to_csv(folder / "falhas_resumo.csv")
        by_depth(payload_limits).to_csv(folder / "limites_so_payload.csv")
    if args.mafft_length.exists():
        if checksum(args.mafft_length) != checksum(args.text[0]):
            raise SystemExit(f"{args.mafft_length} foi gerado com outro arquivo de entrada")
        means, length_limits = mafft_length_summary(pd.read_csv(args.mafft_length), points)
        means.to_csv(folder / "mafft_120.csv")
        by_depth(length_limits).to_csv(folder / "limites_mafft_120.csv")

    for ratio in sorted(runs["ids_ratio"].unique()):
        recovery_figure(points, ratio, folder)
        limit_figure(limits, ratio, folder)
        gain_figure(runs, ratio, folder)
        reconstruction_figure(points, ratio, folder)
    cost_figure(cheapest, folder)
    time_figure(runs, folder)
    print(f"tabelas e figuras em {folder}")


if __name__ == "__main__":
    main()
