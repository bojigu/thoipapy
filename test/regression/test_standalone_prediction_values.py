"""Value-level regression tests for the standalone predictor.

``run_THOIPA_prediction`` is the function the THOIPA webserver calls and the one users call
directly. Until now the only assertions on it were that four output files exist -- the suite would
have passed unchanged if every prediction became 0.5. These tests pin the actual numbers.

Both proteins use homologue tarballs already committed under ``test/test_inputs``, so no network
access is needed. rate4site, freecontact and cd-hit must be installed.
"""

import sys
from pathlib import Path
from shutil import copyfile

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from test.helpers.helpers import TestProtein  # noqa: E402
from test.regression.conftest import assert_frame_matches_golden  # noqa: E402
from thoipapy import run_THOIPA_prediction  # noqa: E402
from thoipapy.utils import make_sure_path_exists  # noqa: E402

REPO_ROOT = Path(__file__).parents[2]


def _run_prediction_with_pre_downloaded_homologues(tp: TestProtein, out_dir: Path) -> pd.DataFrame:
    """Run the standalone predictor offline and return the parsed THOIPA output."""
    make_sure_path_exists(out_dir / "datafiles")
    tarball = REPO_ROOT / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
    assert tarball.is_file(), f"missing pre-downloaded homologues: {tarball}"
    copyfile(tarball, out_dir / "datafiles/BLAST_results.xml.tar.gz")

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=False)

    out_csv = out_dir / "datafiles/THOIPA_full_out.csv"
    assert out_csv.is_file(), "predictor did not produce THOIPA_full_out.csv"
    return pd.read_csv(out_csv, index_col=0)


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
@pytest.mark.parametrize("protein", ["1xioA4", "4hksA1"])
def test_standalone_prediction_values_match_golden(protein, tmp_path):
    tp = TestProtein()
    getattr(tp, f"with_{protein}")()
    df = _run_prediction_with_pre_downloaded_homologues(tp, tmp_path / protein)
    assert_frame_matches_golden(df, f"standalone/{protein}_THOIPA_full_out.csv")


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_standalone_prediction_is_reproducible_within_a_machine(tmp_path):
    """Two runs of the same input must give identical predictions.

    Catches unseeded randomness reintroduced anywhere in the single-protein feature path,
    independently of whether the golden files happen to be present.
    """
    tp = TestProtein()
    tp.with_1xioA4()
    first = _run_prediction_with_pre_downloaded_homologues(tp, tmp_path / "run1")
    second = _run_prediction_with_pre_downloaded_homologues(tp, tmp_path / "run2")
    pd.testing.assert_frame_equal(first, second)


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_residue_numbering_is_one_based(tmp_path):
    """res_num_full_seq must be the 1-based UniProt residue number.

    The training pipeline (common.py) has always used 1-based numbering, while predict.py used
    0-based, so every res_num_full_seq the predictor emitted was one lower than the true residue
    number. This pins the corrected convention.
    """
    tp = TestProtein()
    tp.with_1xioA4()
    df = _run_prediction_with_pre_downloaded_homologues(tp, tmp_path / "numbering")
    # Not a skip. The column is load-bearing: predict.py selects it when building the files the
    # user downloads, so its absence is a failure, not a reason to stop checking.
    assert "res_num_full_seq" in df.columns, "the predictor no longer emits res_num_full_seq"

    expected_start = tp.full_seq.index(tp.tmd_seq) + 1
    actual_start = int(df["res_num_full_seq"].iloc[0])
    assert actual_start == expected_start, (
        f"first TMD residue should be numbered {expected_start} (1-based UniProt), " f"got {actual_start}"
    )
    # and the residue letter at that position must actually match the sequence
    assert tp.full_seq[actual_start - 1] == df["residue_name"].iloc[0]
