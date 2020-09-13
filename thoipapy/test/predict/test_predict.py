import os
from pathlib import Path

import thoipapy
from thoipapy import run_THOIPA_prediction
from thoipapy.predict import get_md5_checksum


def test_predict_erbb3():
    protein_name = "ERBB3"
    TMD_seq = "MALTVIAGLVVIFMMLGGTFL"
    full_seq = "MVQNECRPCHENCTQGCKGPELQDCLGQTLVLIGKTHLTMALTVIAGLVVIFMMLGGTFLYWRGRRIQNKRAMRRYLERGESIEPLDPSEKANKVLA"
    md5 = get_md5_checksum(TMD_seq, full_seq)
    thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    out_dir = Path(thoipapy_module_path) / "test/test_outputs/test_predict"

    run_THOIPA_prediction(protein_name, md5, TMD_seq, full_seq, out_dir, create_heatmap=True, set_number=5)

