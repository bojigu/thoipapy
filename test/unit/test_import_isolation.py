"""Guard against the import bug that silently broke every standalone prediction.

`predict.py` used to do `import thoipapy.features as tf` and then reach for `tf.pssm`,
`tf.entropy` and eight more submodules. That worked only while something else had already
imported them as a side effect. Once every `__init__.py` was emptied, nothing did, and
`import thoipapy; thoipapy.run_THOIPA_prediction(...)` raised AttributeError on the first
feature call.

The full test suite did not catch it. pytest collects `test/functional/` before
`test/regression/`, and `test_ml_pipeline.py` imports `thoipapy.run`, which pulls in all ten
feature submodules. So the regression tests passed on the side effect of an unrelated test
having run first, and failed the moment they ran alone.

These tests use a subprocess so that nothing already imported into the pytest process can mask
the problem. That isolation is the entire point; do not rewrite them as plain imports.

Note which test does the real work. Reintroducing the bug was verified to leave the subprocess
import tests passing and to fail only `test_predict_does_not_rely_on_package_attribute_access`.
`import thoipapy` succeeds either way -- the AttributeError does not surface until a feature
function is actually called, minutes into a prediction. The source-level check is therefore not
redundant with the import checks; it is the only one of them that catches this.
"""

import subprocess
import sys
import textwrap

import pytest


def run_in_fresh_interpreter(code: str) -> subprocess.CompletedProcess:
    """Execute code in a new interpreter, so this process's imports cannot influence the result."""
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        timeout=180,
    )


def test_public_entry_point_is_importable_in_a_fresh_interpreter():
    """`import thoipapy` alone must expose a usable run_THOIPA_prediction."""
    result = run_in_fresh_interpreter("""
        import thoipapy
        assert callable(thoipapy.run_THOIPA_prediction)
        """)
    assert result.returncode == 0, f"`import thoipapy` failed in a clean interpreter:\n{result.stderr}"


@pytest.mark.parametrize(
    "submodule",
    [
        "combine_features",
        "entropy",
        "freecontact",
        "lipophilicity",
        "lips",
        "motifs",
        "physical_parameters",
        "pssm",
        "rate4site",
        "relative_position",
    ],
)
def test_every_feature_submodule_predict_uses_is_importable(submodule):
    """Every feature submodule predict.py calls must resolve without a side-effect import."""
    result = run_in_fresh_interpreter(f"""
        import importlib
        import thoipapy
        importlib.import_module("thoipapy.features.{submodule}")
        """)
    assert result.returncode == 0, f"thoipapy.features.{submodule} not importable:\n{result.stderr}"


def test_predict_does_not_rely_on_package_attribute_access():
    """predict.py must import names directly rather than reaching through a package object.

    `thoipapy.features` has an empty __init__.py by house style, so attribute access on it only
    works by accident. Catch a reintroduction at the source level.
    """
    from pathlib import Path

    import thoipapy

    source = (Path(thoipapy.__file__).parent / "predict.py").read_text()
    assert "import thoipapy.features as" not in source, (
        "predict.py imports the features package as a namespace and reaches through it with "
        "attribute access. thoipapy/features/__init__.py is empty, so this only resolves if "
        "another module happened to import the submodules first. Import the functions directly."
    )
