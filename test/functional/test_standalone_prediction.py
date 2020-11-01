from pathlib import Path
from shutil import rmtree

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy import run_THOIPA_prediction


def test_standalone_prediction_for_tmd_with_few_homologues():
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()

    # cleanup
    rmtree(out_dir)


def test_standalone_prediction_for_tmd_with_many_homologues():
    # tmds with many homologues will require a functional cd-hit before rate4site
    tp: TestProtein = TestProtein()
    tp.with_4ryiA2()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict" / tp.acc

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()

    # cleanup
    rmtree(out_dir)
