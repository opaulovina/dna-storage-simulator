import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "experiments" / "pilot.py"


def test_pilot_runs_isolated_and_parallel_and_reports(tmp_path):
    source = tmp_path / "input.bin"
    source.write_bytes(bytes(range(56)))
    output = tmp_path / "pilot"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(source),
            "--workers", "2",
            "--error-rates", "0.05",
            "--depths", "3",
            "--seeds", "0",
            "--output", str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "resultados idênticos isolado x paralelo: sim" in result.stdout
    metadata = json.loads(output.with_suffix(".json").read_text())
    assert metadata["identical_results"] is True
    assert set(metadata["projection"]) == {"MAFFT", "RS híbrido"}
    assert len(output.with_suffix(".csv").read_text().splitlines()) == 1 + 2 * 2
