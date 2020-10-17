import re
from pathlib import Path

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy.homologues.NCBI_download import download_homologues_from_ncbi
from thoipapy.homologues.NCBI_parser import parse_NCBI_xml_to_csv
from thoipapy.utils import make_sure_path_exists, LogOnlyToConsole


def test_download_homologues_from_ncbi():
    tp: TestProtein = TestProtein()
    tp.with_BNIP3()
    thoipapy_module_path = Path(thoipapy.__file__).parents[1]
    datafiles_dir = thoipapy_module_path / "test/test_outputs/test_download/datafiles"
    blast_xml_file = datafiles_dir / "BLAST_results.xml"
    xml_txt = datafiles_dir / "BLAST_results_details.txt"
    xml_tar_gz = datafiles_dir / "BLAST_results.xml.tar.gz"
    make_sure_path_exists(blast_xml_file, isfile=True)
    expect_value = 0.01
    hit_list_size = 100
    logging = LogOnlyToConsole()
    download_homologues_from_ncbi(tp.acc, tp.full_seq, blast_xml_file, xml_txt, xml_tar_gz, expect_value, hit_list_size, logging)

    assert xml_tar_gz.is_file()

    BLAST_csv_tar = datafiles_dir / "BLAST.csv.tar.gz"
    e_value_cutoff = 0.01
    parse_NCBI_xml_to_csv(tp.acc, xml_tar_gz, BLAST_csv_tar, tp.tmd_start, tp.tmd_end, e_value_cutoff, logging)
    assert BLAST_csv_tar.is_file()


