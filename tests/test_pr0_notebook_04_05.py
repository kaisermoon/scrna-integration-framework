import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _metrics_assignment(notebook_path: Path) -> ast.Assign:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    matches = []
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        if isinstance(source, list):
            source = "".join(source)
        tree = ast.parse(source)
        matches.extend(
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "metrics_str"
                for target in node.targets
            )
        )

    assert len(matches) == 1
    return matches[0]


@pytest.mark.parametrize(
    "notebook_name",
    ["04_embedded.ipynb", "05_clustered.ipynb"],
)
def test_metrics_summary_uses_current_iteration_metrics(notebook_name):
    assignment = _metrics_assignment(PROJECT_ROOT / "notebooks" / notebook_name)
    module = ast.fix_missing_locations(ast.Module(body=[assignment], type_ignores=[]))
    namespace = {
        "_m": {"score": 1.25},
        "np": SimpleNamespace(isnan=math.isnan),
    }

    exec(compile(module, notebook_name, "exec"), namespace)

    assert namespace["metrics_str"] == "score=1.2500"
