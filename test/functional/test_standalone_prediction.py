from pathlib import Path

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy import run_THOIPA_prediction
from thoipapy.predict import get_md5_checksum


def test_standalone_prediction():
    tp: TestProtein = TestProtein()
    tp.with_BNIP3()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    out_dir = thoipapy_module_path / "test/test_outputs/test_predict"

    run_THOIPA_prediction(tp.protein_name, tp.md5, tp.tmd_seq, tp.full_seq, out_dir, create_heatmap=True)

    assert (out_dir / "datafiles").is_dir()
    assert (out_dir / "heatmap.pdf").is_file()
    assert (out_dir / "heatmap.png").is_file()
    assert (out_dir / "THOIPA_out.csv").is_file()
    assert (out_dir / "THOIPA_out.xlsx").is_file()
