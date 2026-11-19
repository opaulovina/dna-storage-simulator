import importlib.util
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = Path(__file__).parents[1] / "app" / "app.py"
SPEC = importlib.util.spec_from_file_location("app", APP)
app_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(app_module)


def test_app_opens_with_the_three_tabs():
    app = AppTest.from_file(str(APP)).run(timeout=60)
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Simular", "Resultados", "Validação (Xie)"]


def test_simulation_shows_both_conditions():
    app = AppTest.from_file(str(APP)).run(timeout=60)
    app.sidebar.slider[0].set_value(0.0)
    app.sidebar.select_slider[0].set_value(1)
    app.button[0].click().run(timeout=300)
    assert not app.exception
    table = app.dataframe[0].value
    assert list(table.index) == ["MAFFT", "BMA + RS 15%"]
    assert list(table["taxa de recuperação"]) == [1.0, 1.0]


def test_reconstruction_can_be_chosen():
    app = AppTest.from_file(str(APP)).run(timeout=60)
    assert app.sidebar.selectbox[0].options == ["BMA", "DBMA", "NW-POA"]


def test_pgm_pixels_survive_a_corrupted_header():
    original = b"P5\n2 2\n255\n" + bytes([10, 20, 30, 40])
    corrupted = b"XX\n2 2\n255\n" + bytes([10, 20, 30, 40])
    image = app_module.pgm_pixels(corrupted, original)
    assert image.size == (2, 2)
    assert image.tobytes() == bytes([10, 20, 30, 40])
