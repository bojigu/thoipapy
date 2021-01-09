from pathlib import Path
from shutil import rmtree

import thoipapy
from test.helpers.helpers import TestProtein
from thoipapy.homologues.NCBI_download import download_homologues_from_ncbi
from thoipapy.homologues.NCBI_parser import parse_NCBI_xml_to_csv
from thoipapy.utils import make_sure_path_exists, LogOnlyToConsole


def test_download_homologues_from_ncbi():
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    test_download_path = Path(thoipapy.__file__).parents[1] / "test/test_outputs/test_download"
    blast_xml_file = test_download_path / f"datafiles/{tp.acc}.surr20.BLAST.xml"
    xml_txt = test_download_path / f"{tp.acc}.surr20.BLAST_details.txt"
    xml_tar_gz = test_download_path / f"{tp.acc}.surr20.BLAST.xml.tar.gz"
    make_sure_path_exists(blast_xml_file, isfile=True)
    expect_value = 0.01
    hit_list_size = 10
    logging = LogOnlyToConsole()
    db = "pdb"
    download_homologues_from_ncbi(tp.acc, tp.full_seq, blast_xml_file, xml_txt, xml_tar_gz, expect_value, hit_list_size, logging, db=db)

    assert xml_tar_gz.is_file()

    BLAST_csv_tar = test_download_path / f"{tp.acc}.surr20.BLAST.csv.tar.gz"
    e_value_cutoff = 0.01
    parse_NCBI_xml_to_csv(tp.acc, xml_tar_gz, BLAST_csv_tar, tp.tmd_start, tp.tmd_end, e_value_cutoff, logging)
    assert BLAST_csv_tar.is_file()

    # cleanup
    rmtree(test_download_path)
