import io
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from PIL import Image

from dna_storage.codec.mapping import RADIX
from dna_storage.codec.payload import BYTES_PER_SEQUENCE, INDEX_LENGTH, encode, rebuild
from dna_storage.correction.reed_solomon import ADDRESSES, RECONSTRUCTIONS, molecule_count
from dna_storage.experiment import CONDITIONS, Trial, evaluate, series_name
from dna_storage.metrics import hamming

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
ASSETS = Path(__file__).resolve().parent / "assets"
MAX_BYTES = RADIX**INDEX_LENGTH * BYTES_PER_SEQUENCE
IMAGE_SUFFIXES = {
    ".png", ".jpg", ".jpeg", ".bmp", ".gif", ".tif", ".tiff", ".webp", ".pgm", ".ppm", ".pnm"
}
AUDIO_SUFFIXES = {".wav", ".mp3", ".ogg", ".flac", ".m4a"}
VIDEO_SUFFIXES = {".mp4", ".webm", ".mov"}
TEXT_SUFFIXES = {
    ".txt", ".csv", ".tsv", ".md", ".json", ".xml", ".html", ".htm", ".fasta", ".fa", ".fastq",
    ".fq", ".py", ".log", ".yaml", ".yml", ".toml", ".ini",
}
DIAGNOSTICS = {
    "pilot",
    "validation_xie",
    "falhas",
    "mafft_120",
    "organick_cobertura",
    "organick_fronteira",
    "perfil_posicional",
}
TEXT_PREVIEW = 3000
MAP_LIMIT = 2000
DEPTHS = (1, 5, 10, 20, 30, 50, 100)
METRICS = {
    "Taxa de recuperação": "recovery_rate",
    "Distância de Hamming normalizada": "hamming",
    "Tempo (s)": "seconds",
    "Memória (KiB)": "memory_kib",
}
SEED_HELP = (
    "A semente define o sorteio dos erros. Com a mesma semente e os mesmos ajustes, o simulador "
    "gera exatamente as mesmas leituras, e o resultado pode ser repetido. Troque o número para "
    "ver outro sorteio."
)
HIDE_STREAMLIT_PARTS = """
<style>
[data-testid="stStatusWidget"] {visibility: hidden;}
[data-testid="stFileUploaderDropzoneInstructions"] {display: none;}
</style>
"""


PGM_HEADER = re.compile(rb"P5\s+(\d+)\s+(\d+)\s+(\d+)\s")


@st.cache_data(show_spinner=False)
# monta o trial e chama evaluate para cada condicao
# devolve um dicionario, o cache guarda o resultado
# entao os mesmos parametros nao rodam de novo
def run(data, error_rate, depth, seed, ratio, redundancy, reconstruction):
    trial = Trial(error_rate, depth, seed, ratio, redundancy, reconstruction)
    # devolve a linha de métricas e as sequências recuperadas para cada condição
    return {
        series_name(condition, reconstruction, redundancy): evaluate(data, trial, condition)
        for condition in CONDITIONS
    }


def estimated_seconds(sequences, depth):
    return sequences * (0.072 + 0.0008 * depth)


def duration(seconds):
    if seconds < 60:
        return f"{seconds:.0f} s"
    if seconds < 3600:
        return f"{seconds / 60:.0f} min"
    return f"{seconds / 3600:.1f} h".replace(".", ",")


# le largura e altura do cabecalho do original com uma expressao regular
# depois monta a imagem com os bytes recuperados que vem depois do cabecalho
# so é usada quando o arquivo recuperado nao abre
def pgm_pixels(recovered, original):
    header = PGM_HEADER.match(original)
    if header is None:
        return None
    width, height = int(header.group(1)), int(header.group(2))
    pixels = recovered[header.end() :].ljust(width * height, b"\0")[: width * height]
    return Image.frombytes("L", (width, height), pixels)


# na barra lateral, envia um arquivo
def choose_input():
    upload = st.sidebar.file_uploader("Arquivo")
    st.sidebar.caption("Qualquer formato, até 1,75 MB.")
    if upload is None:
        return None, None
    data = upload.getvalue()
    if not data:
        st.sidebar.error("O arquivo está vazio.")
        return None, None
    if len(data) > MAX_BYTES:
        st.sidebar.error(
            "O arquivo passa de 1,75 MB. O índice de 8 bases numera no máximo 65.536 sequências "
            "de 28 bytes."
        )
        return None, None
    return upload.name, data


def kind(name, original):
    suffix = Path(name).suffix.lower()
    if suffix in IMAGE_SUFFIXES:
        return "image"
    if suffix in AUDIO_SUFFIXES:
        return "audio"
    if suffix in VIDEO_SUFFIXES:
        return "video"
    if suffix in TEXT_SUFFIXES:
        return "text"
    try:
        original.decode("utf-8")
    except UnicodeDecodeError:
        return "other"
    return "text"


def show_image(data, original, width):
    try:
        image = Image.open(io.BytesIO(data))
        image.load()
    except Exception:
        image = pgm_pixels(data, original)
        if image is None:
            st.warning("O arquivo recuperado não abre como imagem.")
            return
        st.caption(
            "O cabeçalho do PGM (sequência 0) veio corrompido. Os pixels aparecem com as "
            "dimensões do original; o arquivo baixado não abre num visualizador comum."
        )
    if image.width < width:
        height = round(width * image.height / image.width)
        image = image.resize((width, height), Image.NEAREST)
    st.image(image, width=width)


# se a extensao for de imagem, tentar abrir com o pillow
# se falhar usa pgm_pixels e avisa
def show_file(name, data, original, width=256):
    content = kind(name, original)
    if content == "image":
        show_image(data, original, width)
    elif content == "audio":
        st.audio(data)
    elif content == "video":
        st.video(data)
    elif content == "text":
        text = data.decode("utf-8", errors="replace")
        st.text(text[:TEXT_PREVIEW])
        if len(text) > TEXT_PREVIEW:
            st.caption(f"Mostrando os primeiros {TEXT_PREVIEW} caracteres.")
    else:
        st.caption("Este formato não tem pré-visualização. Baixe o arquivo para abrir.")


# para cada condição, calcula o hamming de cada sequencia (0 se identica)
# e desenha uma faixa por condicao, com uma coluna por sequencia
# branco é recuperada, vermelho é errada
def sequence_map(results, sequences):
    exact = len(sequences) <= MAP_LIMIT
    rows = [
        [
            0.0 if new == old else (hamming(new, old) if exact else 0.05)
            for new, old in zip(recovered, sequences, strict=True)
        ]
        for _, recovered in results.values()
    ]
    colors = plt.get_cmap("Reds").with_extremes(under="white")
    figure, axis = plt.subplots(figsize=(10, 1.6))
    axis.imshow(rows, cmap=colors, vmin=1e-9, vmax=0.05, aspect="auto", interpolation="nearest")
    axis.set_yticks(range(len(results)), list(results))
    axis.set_xlabel("sequência (índice)")
    figure.tight_layout()
    return figure, exact


def sidebar_controls():
    st.sidebar.header("Canal de erro")
    error_rate = st.sidebar.slider("Taxa de erro por base", 0.0, 0.30, 0.08, 0.01)
    depth = st.sidebar.select_slider("Profundidade (leituras por sequência)", DEPTHS, value=10)
    st.sidebar.markdown("Proporção I:D:S")
    st.sidebar.caption("Pesos de inserção (I), deleção (D) e substituição (S).")
    insertion, deletion, substitution = st.sidebar.columns(3)
    ratio = (
        insertion.number_input("I", min_value=0, value=1, step=1),
        deletion.number_input("D", min_value=0, value=1, step=1),
        substitution.number_input("S", min_value=0, value=1, step=1),
    )
    seed_column, info_column = st.sidebar.columns([3, 2], vertical_alignment="bottom")
    seed = seed_column.number_input("Semente", min_value=0, value=0, step=1)
    with info_column.popover(":material/info:"):
        st.write(SEED_HELP)
    st.sidebar.header("Condição B")
    reconstruction = st.sidebar.selectbox("Reconstrução", list(RECONSTRUCTIONS))
    redundancy = st.sidebar.slider("Redundância Reed-Solomon", 0.0, 0.50, 0.15, 0.05)
    return {
        "error_rate": error_rate,
        "depth": depth,
        "seed": int(seed),
        "ratio": tuple(int(part) for part in ratio),
        "redundancy": redundancy,
        "reconstruction": reconstruction,
    }


# le os controles na barra lateral;
# mostra quantos bytes e quantas sequencias tem o arquivo;
# ao clicar em simular chama run;
# monta a tabela de metricas;
# desenha o mapa de erros;
# mostra o original e os dois recuperados - rebuild. com botao de download para cada um.
def simulate_tab(name, data):
    settings = sidebar_controls()
    if data is None:
        st.info("Envie um arquivo na barra lateral para começar.")
        return
    if sum(settings["ratio"]) == 0:
        st.error("A proporção I:D:S precisa de pelo menos um valor maior que zero.")
        return
    sequences = encode(data)
    molecules = molecule_count(len(sequences), settings["redundancy"])
    if molecules > ADDRESSES:
        st.error(
            f"Com redundância de {settings['redundancy']:.0%}, a Condição B precisa de "
            f"{molecules:,} moléculas, e o índice de 8 bases numera no máximo {ADDRESSES:,}. "
            "Diminua a redundância ou use um arquivo menor.".replace(",", ".")
        )
        return
    seconds = estimated_seconds(len(sequences), settings["depth"])
    st.write(f"**{name}**: {len(data)} bytes, {len(sequences)} sequências de 120 bases.")
    st.caption(
        f"Tempo estimado no computador dos experimentos: cerca de {duration(seconds)}, quase todo "
        "gasto pelo MAFFT."
    )
    if not st.button("Simular", type="primary"):
        return
    with st.spinner("Simulando..."):
        results = run(data, **settings)

    table = pd.DataFrame([row for row, _ in results.values()], index=list(results))
    table.index.name = "condição"
    st.subheader("Métricas")
    st.dataframe(
        table[["recovery_rate", "hamming", "molecules", "density", "seconds", "memory_kib"]]
        .rename(
            columns={
                "recovery_rate": "taxa de recuperação",
                "hamming": "Hamming normalizada",
                "molecules": "moléculas",
                "density": "densidade (bits/nt)",
                "seconds": "tempo (s)",
                "memory_kib": "memória (KiB)",
            }
        ),
        width="stretch",
    )
    st.caption(
        "Tempo e memória mudam com o computador e com o que mais estiver aberto. Os números da "
        "monografia vêm do modo experimental."
    )

    st.subheader("Onde estão os erros")
    figure, exact = sequence_map(results, sequences)
    st.pyplot(figure)
    if exact:
        st.caption(
            "Cada coluna é uma sequência. Branco: voltou igual à original. Vermelho: voltou com "
            "erros; quanto mais escuro, maior a distância de Hamming (a escala vai até 0,05)."
        )
    else:
        st.caption(
            f"Cada coluna é uma sequência. Com mais de {MAP_LIMIT} sequências, o mapa só "
            "separa as que voltaram iguais (branco) das que voltaram com erros (vermelho)."
        )

    st.subheader("Arquivo recuperado")
    columns = st.columns(1 + len(results))
    with columns[0]:
        st.markdown("**Original**")
        show_file(name, data, data)
    for column, letter, (series, (_, recovered)) in zip(
        columns[1:], "AB", results.items(), strict=True
    ):
        recovered_file = rebuild(recovered, len(data))
        with column:
            st.markdown(f"**{series}**")
            show_file(name, recovered_file, data)
            st.download_button(
                "Baixar",
                recovered_file,
                file_name=f"recuperado_{letter}_{name}",
                key=letter,
                on_click="ignore",
            )


# le um ou mais CSVs. preenche a coluna reconstruction se o CSV for antigo
# cria a coluna series com series_name e remove linhas repetidas entre arquivos
def load_results(paths):
    tables = []
    for path in paths:
        table = pd.read_csv(path)
        if "reconstruction" not in table:
            table["reconstruction"] = "BMA"
        table["series"] = [
            series_name(*values)
            for values in zip(
                table["condition"], table["reconstruction"], table["redundancy"], strict=True
            )
        ]
        tables.append(table)
    table = pd.concat(tables, ignore_index=True)
    keys = ["series", "ids_ratio", "error_rate", "depth", "seed"]
    return table.drop_duplicates(subset=keys)


# le o SHA-256 do JSON correspondente
# serve para impedir que se combinem resultados de arquivos de entrada diferentes
def input_checksum(path):
    metadata = path.with_suffix(".json")
    return json.loads(metadata.read_text()).get("sha256") if metadata.exists() else None


# escolhe os experimentos, a metrica, a proporcao, as profundidades
# e as series, tira a media das sementes e desenha um grafico
# com a cor pela profundidade e o estilo da linha pela serie
# embaixo mostra a tabela
def results_tab():
    files = sorted(RESULTS.glob("*.csv"), key=lambda file: (file.stem != "final_texto", file.stem))
    files = [file for file in files if file.stem not in DIAGNOSTICS]
    if not files:
        st.info("Rode experiments/run_experiments.py para gerar resultados.")
        return
    paths = st.multiselect(
        "Experimentos (podem ser combinados se usarem o mesmo arquivo de entrada)",
        files,
        default=files[:1],
        format_func=lambda file: file.stem,
    )
    if not paths:
        return
    if len({input_checksum(path) for path in paths}) > 1:
        st.error("Os experimentos escolhidos usaram arquivos de entrada diferentes.")
        return
    table = load_results(paths)
    left, middle, right = st.columns(3)
    metric = left.selectbox("Métrica", list(METRICS))
    ratio = middle.selectbox("Proporção I:D:S", sorted(table["ids_ratio"].unique()))
    depths = right.multiselect(
        "Profundidades", sorted(table["depth"].unique()), default=sorted(table["depth"].unique())
    )
    available = list(dict.fromkeys(table["series"]))
    names = st.multiselect("Séries", available, default=available[:4])
    subset = table[
        (table["ids_ratio"] == ratio) & table["depth"].isin(depths) & table["series"].isin(names)
    ]
    means = subset.groupby(["series", "depth", "error_rate"])[METRICS[metric]].mean()

    figure, axis = plt.subplots(figsize=(9, 5))
    colors = plt.colormaps["viridis"].resampled(max(len(depths), 2))
    styles = ["-", "--", ":", "-."]
    for position, depth in enumerate(depths):
        for index, name in enumerate(names):
            if (name, depth) not in means.index.droplevel("error_rate"):
                continue
            series = means.loc[(name, depth)]
            axis.plot(
                series.index,
                series.values,
                styles[index % len(styles)],
                marker="o",
                color=colors(position),
                label=f"{name}, {depth} leituras",
            )
    axis.set_xlabel("taxa de erro por base")
    axis.set_ylabel(metric)
    axis.grid(alpha=0.3)
    axis.legend(fontsize="small", loc="upper left", bbox_to_anchor=(1.01, 1))
    st.pyplot(figure)
    st.caption("Cor: profundidade. Estilo da linha: condição (veja a legenda).")
    pivot = means.unstack(["series", "depth"]).round(3)
    pivot.index.name = "taxa de erro"
    pivot.columns.names = ["série", "leituras"]
    st.dataframe(pivot, width="stretch")


def main():
    st.set_page_config(
        page_title="Simulador de correção de erros em DNA",
        page_icon=str(ASSETS / "favicon.png"),
        layout="wide",
    )
    st.markdown(HIDE_STREAMLIT_PARTS, unsafe_allow_html=True)
    st.sidebar.image(str(ASSETS / "logo.png"), width=120)
    st.sidebar.title("Simulação")
    name, data = choose_input()

    st.title("Simulador de correção de erros em DNA")
    st.caption(
        "Compara o MAFFT (Condição A) e o Reed-Solomon híbrido (Condição B) nas mesmas leituras."
    )
    simulate, results = st.tabs(["Simular", "Resultados"])
    with simulate:
        simulate_tab(name, data)
    st.sidebar.caption("made by: Paulo Vinícius Silva - 2026")
    with results:
        results_tab()


if __name__ == "__main__":
    main()
