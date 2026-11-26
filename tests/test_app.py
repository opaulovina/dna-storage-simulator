import importlib.util
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = Path(__file__).parents[1] / "app" / "app.py"
SPEC = importlib.util.spec_from_file_location("app", APP)
app_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app_module)


def test_app_opens_asking_for_a_file():
    app = AppTest.from_file(str(APP)).run(timeout=60)
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Simular", "Resultados"]
    assert "Envie um arquivo" in app.info[0].value


def test_sidebar_has_the_channel_and_condition_b_controls():
    app = AppTest.from_file(str(APP)).run(timeout=60)
    assert app.sidebar.selectbox[0].options == ["BMA", "DBMA", "NW-POA"]
    assert [box.label for box in app.sidebar.number_input] == ["I", "D", "S", "Semente"]


def test_simulation_shows_both_conditions():
    app = AppTest.from_file(str(APP), default_timeout=120).run()
    app.sidebar.file_uploader[0].set_value(("teste.txt", bytes(range(56)), "text/plain"))
    app.run()
    app.sidebar.slider[0].set_value(0.0)
    app.sidebar.select_slider[0].set_value(1)
    app.run()
    app.button[0].click().run()
    assert not app.exception
    table = app.dataframe[0].value
    assert list(table.index) == ["MAFFT", "BMA + RS 15%"]
    assert list(table["taxa de recuperação"]) == [1.0, 1.0]
    assert [button.label for button in app.get("download_button")] == ["Baixar", "Baixar"]


def test_empty_and_oversized_files_are_refused():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    app.sidebar.file_uploader[0].set_value(("vazio.txt", b"", "text/plain")).run()
    assert "vazio" in app.sidebar.error[0].value
    big = bytes(app_module.MAX_BYTES + 1)
    app.sidebar.file_uploader[0].set_value(("grande.bin", big, "application/octet-stream")).run()
    assert "1,75 MiB" in app.sidebar.error[0].value


def test_footer_is_at_the_bottom_of_the_sidebar():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    assert app.sidebar.caption[-1].value == "made by: Paulo Vinícius Silva - 2026"


def test_file_kinds_cover_common_formats():
    assert app_module.kind("foto.JPG", b"") == "image"
    assert app_module.kind("som.wav", b"") == "audio"
    assert app_module.kind("leituras.fastq", b"") == "text"
    assert app_module.kind("sem_extensao", b"texto comum") == "text"
    assert app_module.kind("planilha.xlsx", bytes([0xFF, 0xFE, 0x00, 0x81])) == "other"


def test_upload_limit_is_the_capacity_of_the_index():
    assert app_module.MAX_BYTES == 65_536 * 28


def test_pgm_pixels_survive_a_corrupted_header():
    original = b"P5\n2 2\n255\n" + bytes([10, 20, 30, 40])
    corrupted = b"XX\n2 2\n255\n" + bytes([10, 20, 30, 40])
    image = app_module.pgm_pixels(corrupted, original)
    assert image.size == (2, 2)
    assert image.tobytes() == bytes([10, 20, 30, 40])


def test_a_valid_but_wrong_pgm_header_keeps_the_original_dimensions():
    original = b"P5\n4 2\n255\n" + bytes(range(8))
    shifted = b"P5\n4 1\n255\n" + bytes(range(8))
    image, damaged = app_module.open_image(shifted, original)
    assert damaged
    assert image.size == (4, 2)
    assert image.tobytes() == bytes(range(8))


def test_an_intact_image_opens_normally():
    original = b"P5\n4 2\n255\n" + bytes(range(8))
    image, damaged = app_module.open_image(original, original)
    assert not damaged
    assert image.size == (4, 2)


def test_condition_b_refuses_more_molecules_than_the_index_addresses():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    big = bytes(60_000 * 28)
    app.sidebar.file_uploader[0].set_value(("grande.bin", big, "application/octet-stream")).run()
    assert "65.536" in app.error[0].value
