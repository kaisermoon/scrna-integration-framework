import ast
import json
from pathlib import Path

NOTEBOOK_PATH = Path(__file__).parents[1] / "notebooks" / "06b_per_cluster.ipynb"


def _params_cell_source() -> str:
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    params_cells = [
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
        and "# === PARAMS ===" in "".join(cell["source"])
    ]
    assert len(params_cells) == 1
    return params_cells[0]


def test_06b_params_cell_runs_in_fresh_namespace() -> None:
    source = _params_cell_source()
    namespace: dict[str, object] = {}

    exec(compile(source, str(NOTEBOOK_PATH), "exec"), namespace)

    assert namespace["OUTPUT_VERSION"] == "v1"
    assert namespace["UPSTREAM_PATH"] == "results/06_annotated_v1.h5ad"
    assert namespace["LABEL_COL"] == "cell_type_final_v1"

    assignments = [
        node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Name)
        and isinstance(node.ctx, ast.Store)
        and node.id == "OUTPUT_VERSION"
    ]
    assert len(assignments) == 1
