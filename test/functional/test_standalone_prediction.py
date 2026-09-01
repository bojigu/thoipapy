from pathlib import Path
from shutil import copyfile, rmtree

import pandas as pd
import pytest
import thoipapy
from thoipapy import run_THOIPA_prediction
from thoipapy.utils import make_sure_path_exists

from test.helpers.helpers import TestProtein


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_standalone_prediction_with_pre_downloaded_homologues_1xioA4():
    # tmd with few homologues, not requiring cd-hit
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / f"test/test_outputs/pre_downloaded_homologues_{tp.acc}"
    make_sure_path_exists(out_dir / "datafiles")
    pre_downloaded_homologue_xml_tar_gz: Path = (
        thoipapy_module_path / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
    )
    assert pre_downloaded_homologue_xml_tar_gz.is_file()
    temp_homologue_xml_tar_gz: Path = out_dir / "datafiles/BLAST_results.xml.tar.gz"
    copyfile(pre_downloaded_homologue_xml_tar_gz, temp_homologue_xml_tar_gz)

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)
    assert_the_output_a_user_receives_is_usable(out_dir, tp)

    # cleanup
    rmtree(out_dir)


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_standalone_prediction_with_pre_downloaded_homologues_4hksA1():
    # tmd with medium number of homologues, requiring cd-hit
    tp: TestProtein = TestProtein()
    tp.with_4hksA1()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / f"test/test_outputs/pre_downloaded_homologues_{tp.acc}"
    make_sure_path_exists(out_dir / "datafiles")
    pre_downloaded_homologue_xml_tar_gz: Path = (
        thoipapy_module_path / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
    )
    assert pre_downloaded_homologue_xml_tar_gz.is_file()
    temp_homologue_xml_tar_gz: Path = out_dir / "datafiles/BLAST_results.xml.tar.gz"
    copyfile(pre_downloaded_homologue_xml_tar_gz, temp_homologue_xml_tar_gz)

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)
    assert_the_output_a_user_receives_is_usable(out_dir, tp)

    # cleanup
    rmtree(out_dir)


@pytest.mark.network
@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_standalone_prediction_for_tmd_with_few_homologues():
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)
    assert_the_output_a_user_receives_is_usable(out_dir, tp)

    # cleanup
    rmtree(out_dir)


@pytest.mark.network
@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_standalone_prediction_for_tmd_with_many_homologues():
    # tmds with many homologues will require a functional cd-hit before rate4site
    tp: TestProtein = TestProtein()
    tp.with_4ryiA2()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)
    assert_the_output_a_user_receives_is_usable(out_dir, tp)

    # cleanup
    rmtree(out_dir)


def assert_output_files_exist(out_dir):
    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()


def assert_the_output_a_user_receives_is_usable(out_dir, tp: TestProtein):
    """The two files a user actually downloads, checked as a reader would open them.

    THOIPA_out.csv was named .csv, written with tabs, and padded with spaces so it lined up in a
    terminal, so anything that guesses the separator from the extension -- Excel, mail-client
    previews -- showed the whole prediction as one column. And the residue number in the submitted
    protein was computed, written to the internal file, and then dropped from both of these.
    """
    expected_columns = ["residue number", "residue name", "residue number in full sequence", "THOIPA"]

    for path in [out_dir / "THOIPA_out.csv", out_dir / "THOIPA_out.xlsx"]:
        df = pd.read_csv(path) if path.suffix == ".csv" else pd.read_excel(path)
        assert list(df.columns) == expected_columns, f"{path.name} columns"
        assert len(df) == len(tp.tmd_seq), f"{path.name} should have one row per TMD residue"
        assert "".join(df["residue name"].tolist()) == tp.tmd_seq, f"{path.name} residues"
        assert list(df["residue number"]) == list(range(1, len(tp.tmd_seq) + 1)), f"{path.name} TMD numbering"

        # the residue number in the submitted protein, 1-based, and it must point at the right residue
        first = int(df["residue number in full sequence"].iloc[0])
        assert first == tp.full_seq.index(tp.tmd_seq) + 1, f"{path.name} full-sequence numbering"
        assert tp.full_seq[first - 1] == df["residue name"].iloc[0]

        assert df["THOIPA"].between(0, 1).all(), f"{path.name} scores must be probabilities"


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_a_tmd_that_occurs_twice_in_the_full_sequence_raises(tmp_path):
    """It used to log a warning and return, having written nothing.

    A caller could not tell that apart from success: the webserver waited for output files that
    were never going to appear.
    """
    tp = TestProtein()
    tp.with_1xioA4()
    doubled_full_seq = tp.full_seq + tp.tmd_seq

    with pytest.raises(ValueError, match="occurs 2 times"):
        run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, doubled_full_seq, tmp_path, create_heatmap=False)


@pytest.mark.requires_tool("rate4site", "freecontact", "cd-hit")
def test_a_tmd_that_is_not_in_the_full_sequence_raises(tmp_path):
    tp = TestProtein()
    tp.with_1xioA4()

    with pytest.raises(ValueError, match="was not found in the full protein sequence"):
        run_THOIPA_prediction(
            tp.protein_name, tp.md5, "WWWWWWWWWWWWWWWWWWWW", tp.full_seq, tmp_path, create_heatmap=False
        )
