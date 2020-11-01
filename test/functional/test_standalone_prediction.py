from pathlib import Path
from shutil import rmtree, copyfile

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy import run_THOIPA_prediction
from thoipapy.utils import make_sure_path_exists


def test_standalone_prediction_with_pre_downloaded_homologues_1xioA4():
    # tmd with few homologues, not requiring cd-hit
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / f"test/test_outputs/pre_downloaded_homologues_{tp.acc}"
    make_sure_path_exists(out_dir / "datafiles")
    pre_downloaded_homologue_xml_tar_gz: Path = thoipapy_module_path / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
    assert pre_downloaded_homologue_xml_tar_gz.is_file()
    temp_homologue_xml_tar_gz: Path = out_dir / "datafiles/BLAST_results.xml.tar.gz"
    copyfile(pre_downloaded_homologue_xml_tar_gz, temp_homologue_xml_tar_gz)

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)

    # cleanup
    rmtree(out_dir)


def test_standalone_prediction_with_pre_downloaded_homologues_4hksA1():
    # tmd with medium number of homologues, requiring cd-hit
    tp: TestProtein = TestProtein()
    tp.with_4hksA1()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / f"test/test_outputs/pre_downloaded_homologues_{tp.acc}"
    make_sure_path_exists(out_dir / "datafiles")
    pre_downloaded_homologue_xml_tar_gz: Path = thoipapy_module_path / f"test/test_inputs/blast_data_valid/{tp.acc}.surr20.BLAST.xml.tar.gz"
    assert pre_downloaded_homologue_xml_tar_gz.is_file()
    temp_homologue_xml_tar_gz: Path = out_dir / "datafiles/BLAST_results.xml.tar.gz"
    copyfile(pre_downloaded_homologue_xml_tar_gz, temp_homologue_xml_tar_gz)

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)

    # cleanup
    rmtree(out_dir)


def test_standalone_prediction_for_tmd_with_few_homologues():
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)

    # cleanup
    rmtree(out_dir)


def test_standalone_prediction_for_tmd_with_many_homologues():
    # tmds with many homologues will require a functional cd-hit before rate4site
    tp: TestProtein = TestProtein()
    tp.with_4ryiA2()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert_output_files_exist(out_dir)

    # cleanup
    rmtree(out_dir)


def assert_output_files_exist(out_dir):
    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()