#!/usr/bin/env python3
"""Compile all EMC article cells and execute the non-optional ones in order.

Run only on a trusted article: its Python cells are executed as ordinary code.
This checks the numerical examples, not the complete Quarto/HTML render.
"""
from __future__ import annotations

import argparse
import ast
from importlib.metadata import version
from pathlib import Path
import re
import sys


def check_prediction_contract(code: str) -> None:
    """Catch the original observed-value bug without pretending to run NumPyro."""
    tree = ast.parse(code)
    model = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name == "mineral_model")
    defaults = dict(zip([a.arg for a in model.args.args][-len(model.args.defaults):],
                        model.args.defaults))
    assert isinstance(defaults["y"], ast.Constant) and defaults["y"].value is None
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)]
    prediction = next(n for n in calls if isinstance(n.func, ast.Name)
                      and n.func.id == "predict")
    y = next(k.value for k in prediction.keywords if k.arg == "y")
    assert isinstance(y, ast.Constant) and y.value is None
    assert any(isinstance(n, ast.Import) and any(a.name == "arviz" for a in n.names)
               for n in tree.body), "Optional example must import ArviZ itself."


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[1] / "posts/element-to-mineral-conversion/index.qmd"
    parser.add_argument("article", nargs="?", type=Path, default=default)
    parser.add_argument("--figure-dir", type=Path)
    args = parser.parse_args()
    if not args.article.is_file():
        parser.error(f"Article not found: {args.article}")
    text = args.article.read_text(encoding="utf-8")
    assert re.search(r"^  error: false$", text, re.M), "Render must fail on code errors."
    cells = re.findall(r"^```\{python\}\s*\n(.*?)^```\s*$", text, re.M | re.S)
    if not cells:
        raise RuntimeError("No executable Python fences found.")
    if args.figure_dir:
        args.figure_dir.mkdir(parents=True, exist_ok=True)

    import IPython.display
    from plotnine import ggplot
    original_display = IPython.display.display
    current_label = ""
    figure_count = 0

    def capture_display(*objects: object, **kwargs: object) -> None:
        nonlocal figure_count
        for obj in objects:
            if isinstance(obj, ggplot):
                # Draw even when no image output was requested: catch plot errors.
                fig = obj.draw()
                figure_count += 1
                if args.figure_dir:
                    fig.savefig(args.figure_dir / f"{current_label}.png", dpi=160,
                                bbox_inches="tight")
                import matplotlib.pyplot as plt
                plt.close(fig)
                print(f"FIGURE OK: {current_label}")
            elif hasattr(obj, "to_string"):
                print(obj.to_string(index=False))
            else:
                print(obj)

    namespace: dict[str, object] = {"__name__": "__emc_article__"}
    executed, skipped, labels = [], [], set()
    IPython.display.display = capture_display
    try:
        for number, code in enumerate(cells, 1):
            label_match = re.search(r"^#\| label:\s*(.+)$", code, re.M)
            if not label_match:
                raise RuntimeError(f"Cell {number} has no label.")
            current_label = label_match.group(1).strip()
            assert current_label not in labels, "Duplicate code-cell label."
            labels.add(current_label)
            compiled = compile(code, f"{args.article}::{current_label}", "exec")
            if re.search(r"^#\| eval:\s*false\s*$", code, re.M):
                check_prediction_contract(code)
                skipped.append(current_label)
                print(f"SYNTAX / PREDICTION CONTRACT OK; NOT EXECUTED: {current_label}")
                continue
            exec(compiled, namespace)
            executed.append(current_label)
            print(f"PASS: {current_label}")
    finally:
        IPython.display.display = original_display

    # Independent checks of the scientific interpretation, not just code execution.
    import numpy as np
    from scipy.optimize import linprog
    from scipy.stats import beta
    C = namespace["C"]
    v = namespace["v"]
    h = np.array([0., 1., 0.])
    assert h @ v == 1 and np.linalg.matrix_rank(np.vstack([C, h])) == 3
    closure_example = np.array([[1., 0.]])
    assert np.linalg.matrix_rank(np.vstack([closure_example, np.ones(2)])) == 2
    # A rank-deficient equality system can still have one feasible boundary point.
    A_eq = np.array([[0., 1., 1.], [1., 1., 1.]])
    for j in range(3):
        for sign in (-1., 1.):
            fit = linprog(sign * np.eye(3)[j], A_eq=A_eq, b_eq=[0., 1.],
                          bounds=(0., None), method="highs")
            assert fit.success
            np.testing.assert_allclose(fit.x, [1., 0., 0.], atol=1e-9)
    assert not np.isclose(beta.var(2, 8), 0.4 ** 2 * beta.var(2, 2))
    mean = namespace["posterior_mean"]
    lo, hi = namespace["ci_low"], namespace["ci_high"]
    assert 0.26 < mean < 0.28 and 0. < lo < mean < hi < 0.40
    assert len(namespace["sensitivity_rows"]) == 5
    print("PASS: independent phase-row, closure, boundary-uniqueness and conditional-prior checks")
    print(f"RESULT: {len(executed)} main cells executed; {figure_count} figures drawn; "
          f"{len(skipped)} optional cell skipped.")
    print("LIMIT: This is not a full Quarto render or a NumPyro runtime test.")
    print("Environment: Python", sys.version.split()[0], "; ".join(
        f"{package}={version(package)}" for package in
        ("numpy", "scipy", "pandas", "plotnine", "IPython")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
