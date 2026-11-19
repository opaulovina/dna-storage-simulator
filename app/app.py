import io
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from PIL import Image

from dna_storage.codec.payload import encode, rebuild
from dna_storage.correction.reed_solomon import RECONSTRUCTIONS
from dna_storage.experiment import CONDITIONS, Trial, evaluate, series_name
from dna_storage.metrics import hamming

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SAMPLES = {
    "Abstract da monografia (texto, 3.920 bytes)": ROOT / "data" / "abstract_3920.txt",
    "Imagem 128 x 128 (PGM, 16 KB)": ROOT / "data" / "imagem.pgm",
}
MAX_BYTES = 32 * 1024
IMAGE_SUFFIXES = {".pgm", ".png", ".bmp", ".jpg", ".jpeg"}
DEPTHS = (1, 5, 10, 20, 30, 50, 100)
RATIOS = {"1:1:1": (1, 1, 1), "1:1:2": (1, 1, 2)}
METRICS = {
    "Taxa de recuperação": "recovery_rate",
    "Distância de Hamming normalizada": "hamming",
    "Tempo (s)": "seconds",
    "Memória (KiB)": "memory_kib",
}


PGM_HEADER = re.compile(rb"P5\s+(\d+)\s+(\d+)\s+(\d+)\s")


@st.cache_data(show_spinner=False)
# monta o trial e chama evaluate para cada condicao
# devolve um dicionario, o cache guarda o resultado
# entao os mesmos parametros nao rodam de novo
def run(data, error_rate, depth, seed, ratio, redundancy, reconstruction):
    trial = Trial(error_rate, depth, seed, RATIOS[ratio], redundancy, reconstruction)
    # devolve a linha de métricas e as sequências recuperadas para cada condição
    return {
        series_name(condition, reconstruction, redundancy): evaluate(data, trial, condition)
        for condition in CONDITIONS
    }


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


# na barra lateral, escolhe um arquivo de entrada entre os exemplos ou envia um arquivo
def choose_input():
    st.sidebar.title("Simulação")
    source = st.sidebar.radio("Arquivo de entrada", [*SAMPLES, "Enviar um arquivo"])
    if source in SAMPLES:
        path = SAMPLES[source]
        return path.name, path.read_bytes()
    upload = st.sidebar.file_uploader(f"Até {MAX_BYTES // 1024} KB")
    if upload is None:
        return None, None
    data = upload.getvalue()
    if not data or len(data) > MAX_BYTES:
        st.sidebar.error(f"O arquivo precisa ter entre 1 byte e {MAX_BYTES // 1024} KB.")
        return None, None
    return upload.name, data


# se a extensao for de imagem, tentar abrir com o pillow
# se falhar usa pgm_pixels e avisa
# senao, mostra como texto, bytes invalidos viram "�"
def show_file(name, data, original=None, width=256):
    if Path(name).suffix.lower() in IMAGE_SUFFIXES:
        try:
            image = Image.open(io.BytesIO(data))
            image.load()
        except Exception:
            image = pgm_pixels(data, original) if original else None
            if image is None:
                st.warning("O arquivo não abre como imagem.")
                return
            st.caption(
                "O cabeçalho do PGM (sequência 0) foi corrompido. Os pixels são mostrados com as "
                "dimensões do original; o arquivo baixado não abre num visualizador comum."
            )
        height = round(width * image.height / image.width)
        st.image(image.resize((width, height), Image.NEAREST))
    else:
        st.text(data.decode("utf-8", errors="replace"))


# para cada condição, calcula o hamming de cada sequencia (0 se identica)
# e desenha uma faixa por condicao, com uma coluna por sequencia
# branco é recuperada, vermelho é errada
def sequence_map(results, sequences):
    rows = [
        [
            0.0 if new == old else hamming(new, old)
            for new, old in zip(recovered, sequences, strict=True)
        ]
        for _, recovered in results.values()
    ]
    figure, axis = plt.subplots(figsize=(10, 1.6))
    axis.imshow(rows, cmap="Reds", vmin=0, vmax=0.05, aspect="auto", interpolation="nearest")
    axis.set_yticks(range(len(results)), list(results))
    axis.set_xlabel("sequência (índice)")
    figure.tight_layout()
    return figure


# le os controles na barra lateral;
# mostra quantos bytes e quantas sequencias tem o arquivo;
# ao clicar em simular chama run;
# monta a tabela de metricas;
# desenha o mapa de erros;
# mostra o original e os dois recuperados - rebuild. com botao de download para cada um.
def simulate_tab():
    name, data = choose_input()
    st.sidebar.header("Canal de erro")
    error_rate = st.sidebar.slider("Taxa de erro por base", 0.0, 0.30, 0.08, 0.01)
    depth = st.sidebar.select_slider("Profundidade (leituras por sequência)", DEPTHS, value=10)
    ratio = st.sidebar.radio("Proporção I:D:S", list(RATIOS), horizontal=True)
    seed = st.sidebar.number_input("Semente", min_value=0, value=0, step=1)
    st.sidebar.header("Condição B")
    reconstruction = st.sidebar.selectbox("Reconstrução", list(RECONSTRUCTIONS))
    redundancy = st.sidebar.slider("Redundância Reed-Solomon", 0.0, 0.50, 0.15, 0.05)

    if data is None:
        st.info("Escolha ou envie um arquivo na barra lateral.")
        return
    sequences = encode(data)
    st.write(f"**{name}**: {len(data)} bytes, {len(sequences)} sequências de 120 bases.")
    if not st.button("Simular", type="primary"):
        return

    with st.spinner("Gerando leituras e corrigindo nas duas condições..."):
        results = run(data, error_rate, depth, int(seed), ratio, redundancy, reconstruction)

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
        "Tempo e memória aqui são indicativos: a interface divide o computador com o navegador. "
        "Os valores da monografia vêm do modo experimental (experiments/run_experiments.py)."
    )

    st.subheader("Onde estão os erros")
    st.pyplot(sequence_map(results, sequences))
    st.caption(
        "Cada coluna é uma sequência. Branco: recuperada por inteiro; vermelho: com erros "
        "(quanto mais escuro, maior a distância de Hamming; escala até 0,05)."
    )

    st.subheader("Arquivo recuperado")
    columns = st.columns(1 + len(results))
    with columns[0]:
        st.markdown("**Original**")
        show_file(name, data)
    for column, letter, (series, (_, recovered)) in zip(
        columns[1:], "AB", results.items(), strict=True
    ):
        recovered_file = rebuild(recovered, len(data))
        with column:
            st.markdown(f"**{series}**")
            show_file(name, recovered_file, original=data)
            st.download_button(
                "Baixar", recovered_file, file_name=f"recuperado_{letter}_{name}", key=letter
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
    files = [file for file in files if file.stem not in {"pilot", "validation_xie"}]
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
    axis.legend(fontsize="small", ncol=2)
    st.pyplot(figure)
    st.caption("Cor: profundidade. Estilo da linha: condição (veja a legenda).")
    st.dataframe(means.unstack(["series", "depth"]).round(3), width="stretch")


# figura de validacao
def validation_tab():
    figure = RESULTS / "validation_xie.png"
    if figure.exists():
        st.image(str(figure))
        st.caption("Taxa de recuperação da Condição A nos parâmetros de Xie et al. (2023).")
    else:
        st.info("Rode experiments/validate_xie.py para gerar a validação.")


def main():
    st.set_page_config(page_title="Simulador de correção de erros em DNA", layout="wide")
    st.title("Correção de erros em armazenamento de dados em DNA")
    st.caption("MAFFT (Condição A) x Reed-Solomon híbrido (Condição B), sob o mesmo ruído.")
    simulate, results, validation = st.tabs(["Simular", "Resultados", "Validação (Xie)"])
    with simulate:
        simulate_tab()
    with results:
        results_tab()
    with validation:
        validation_tab()


if __name__ == "__main__":
    main()
