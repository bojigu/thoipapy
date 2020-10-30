from pathlib import Path

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy import run_THOIPA_prediction


def test_standalone_prediction_for_tmd_with_few_homologues():
    tp: TestProtein = TestProtein()
    tp.with_BNIP3()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict/BNIP3"

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()


def test_standalone_prediction_for_tmd_with_many_homologues():
    tp: TestProtein = TestProtein()
    tp.with_4ryiA2()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict/4ryiA2"

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()
