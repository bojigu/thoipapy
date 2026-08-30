"""Regenerate the golden reference files used by the regression tests.

Run this ONLY when a change is meant to alter what the pipeline computes. It reports a diff
against the existing references so the change can be reviewed before it is accepted.

    python test/regression/regenerate_golden.py

Then review the printed diff, and if it is what you intended::

    dvc add test/regression_data
    dvc push
    git add test/regression_data.dvc

Commit the updated .dvc file in the same commit as the code change that caused it, so the reason
for the new reference values is recorded next to them.

Hardcoded paths, no argparse, per the house style for pipeline scripts.
"""

import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(REPO_ROOT))

from test.helpers.helpers import TestProtein  # noqa: E402
from thoipapy import run_THOIPA_prediction  # noqa: E402
from thoipapy.utils import make_sure_path_exists  # noqa: E402

GOLDEN_DIR = REPO_ROOT / "test/regression_data/golden"
STANDALONE_PROTEINS = ["1xioA4", "4hksA1"]


def regenerate_standalone(workdir: Path) -> list[tuple[str, pd.DataFrame]]:
    """Run the standalone predictor for each test protein and return (relative_path, df)."""
    produced = []
    for protein in STANDALONE_PROTEINS:
        tp = TestProtein()
        getattr(tp, f"with_{protein}")()
        out_dir = workdir / protein
        make_sure_path_exists(out_dir / "datafiles")
        tarball = REPO_ROOT / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
        shutil.copyfile(tarball, out_dir / "datafiles/BLAST_results.xml.tar.gz")
        print(f"  running standalone prediction for {protein} ...")
        run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=False)
        # THOIPA_full_out.csv is the machine-readable output: full float precision and every
        # feature column, so it regresses the whole single-protein feature path, not just the
        # final probability. THOIPA_out.csv is whitespace-padded to 3 decimals for humans.
        df = pd.read_csv(out_dir / "datafiles/THOIPA_full_out.csv", index_col=0)
        produced.append((f"standalone/{protein}_THOIPA_full_out.csv", df))
    return produced


def report_diff(relative_path: str, new: pd.DataFrame) -> bool:
    """Print how the new output differs from the existing golden. Returns True if it changed."""
    old_path = GOLDEN_DIR / relative_path
    if not old_path.is_file():
        print(f"  NEW      {relative_path}  ({len(new)} rows, {len(new.columns)} cols)")
        return True

    old = pd.read_csv(old_path, index_col=0)
    if list(old.columns) != list(new.columns):
        print(f"  COLUMNS  {relative_path}")
        print(f"    was: {list(old.columns)}")
        print(f"    now: {list(new.columns)}")
        return True
    if len(old) != len(new):
        print(f"  ROWS     {relative_path}: {len(old)} -> {len(new)}")
        return True

    changed = False
    for col in new.columns:
        if pd.api.types.is_numeric_dtype(new[col]):
            delta = abs(new[col].to_numpy() - old[col].to_numpy())
            if delta.max() > 0:
                print(
                    f"  CHANGED  {relative_path} [{col}]  max|delta| = {delta.max():.3e}  " f"mean = {delta.mean():.3e}"
                )
                changed = True
        else:
            n = (old[col].astype(str) != new[col].astype(str)).sum()
            if n:
                print(f"  CHANGED  {relative_path} [{col}]  {n} rows differ")
                changed = True
    if not changed:
        print(f"  same     {relative_path}")
    return changed


def main() -> int:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        produced = regenerate_standalone(Path(tmp))

        print("\nDiff against existing golden references:")
        any_changed = False
        for relative_path, df in produced:
            any_changed |= report_diff(relative_path, df)

        for relative_path, df in produced:
            target = GOLDEN_DIR / relative_path
            target.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(target)

    print(f"\nWrote {len(produced)} reference files to {GOLDEN_DIR}")
    if any_changed:
        print("\nOutput CHANGED. Review the diff above. If it is intended:")
        print("  dvc add test/regression_data && dvc push && git add test/regression_data.dvc")
    else:
        print("\nNo change against the existing references.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
