"""Shared helpers for the regression (golden-file) tests.

These tests exist to prove that refactoring does not change what the pipeline computes. They
compare freshly computed output against reference files under ``test/regression_data/golden/``,
which are DVC-tracked so CI can fetch them without carrying them in git.

Regenerating the references is deliberate, never automatic. When a change is *meant* to alter the
output, run::

    python test/regression/regenerate_golden.py

review the reported diff, then ``dvc add test/regression_data && dvc push`` and commit the updated
``.dvc`` file together with the code change that caused it.
"""

from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parents[2]
REGRESSION_DATA = REPO_ROOT / "test/regression_data"
GOLDEN_DIR = REGRESSION_DATA / "golden"

# Predicted probabilities are floating point sums over 50 trees. Bitwise equality is not a
# reasonable expectation across BLAS builds, but anything above this tolerance is a real change.
PREDICTION_TOLERANCE = 1e-9


def require_golden(relative_path: str) -> Path:
    """Return a golden file path, skipping the test with an actionable message if it is absent."""
    path = GOLDEN_DIR / relative_path
    if not path.is_file():
        pytest.skip(
            f"golden reference {relative_path} not present. Run 'dvc pull' to fetch "
            f"test/regression_data, or 'python test/regression/regenerate_golden.py' to create it."
        )
    return path


def assert_frame_matches_golden(
    actual: pd.DataFrame, relative_path: str, tolerance: float = PREDICTION_TOLERANCE
) -> None:
    """Compare a dataframe against its golden reference, with a diff that says what moved.

    pandas' own assertion message is not much help when a numeric column has drifted, so the
    largest offending column and row are reported explicitly.
    """
    expected = pd.read_csv(require_golden(relative_path), index_col=0)

    assert list(actual.columns) == list(expected.columns), (
        f"{relative_path}: column names or order changed.\n"
        f"  expected: {list(expected.columns)}\n"
        f"  actual:   {list(actual.columns)}"
    )
    assert len(actual) == len(expected), f"{relative_path}: row count changed, {len(expected)} -> {len(actual)}"

    numeric = [c for c in expected.columns if pd.api.types.is_numeric_dtype(expected[c])]
    worst_col, worst_delta, worst_row = None, 0.0, None
    for col in numeric:
        delta = actual[col].to_numpy() - expected[col].to_numpy()
        abs_delta = abs(delta)
        if abs_delta.max() > worst_delta:
            worst_delta = float(abs_delta.max())
            worst_col = col
            worst_row = str(actual.index[abs_delta.argmax()])

    assert worst_delta <= tolerance, (
        f"{relative_path}: numeric output changed beyond tolerance {tolerance:g}.\n"
        f"  largest change in column '{worst_col}' at row '{worst_row}': {worst_delta:g}\n"
        f"  If this change is intended, run test/regression/regenerate_golden.py, review the "
        f"diff, then dvc add + dvc push test/regression_data."
    )

    for col in expected.columns:
        if col not in numeric:
            mismatches = (actual[col].astype(str) != expected[col].astype(str)).sum()
            assert mismatches == 0, f"{relative_path}: column '{col}' differs in {mismatches} rows"
