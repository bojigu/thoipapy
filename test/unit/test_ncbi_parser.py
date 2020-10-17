from pathlib import Path
from shutil import rmtree

from test.helpers.helpers import TestProtein
from thoipapy.homologues.NCBI_parser import parse_NCBI_xml_to_csv
from thoipapy.utils import make_sure_path_exists, LogOnlyToConsole


def test_parse_NCBI_xml_to_csv():
    # example xml tarball includes "CREATE_VIEW" string indicating malformed xml that needs cleaning, as often retrieved from NCBI
    blast_xml_tar: Path = Path(__file__).parents[1] / "test_inputs/blast_data/BLAST_results.xml.tar.gz"
    blast_csv_tar = Path(__file__).parents[1] / "test_outputs/test_parse_NCBI_xml_to_csv/BLAST_results.csv.tar.gz"
    make_sure_path_exists(blast_csv_tar, isfile=True)
    tp: TestProtein = TestProtein()
    tp.with_BNIP3()
    logging = LogOnlyToConsole()
    assert not blast_csv_tar.is_file()
    parse_NCBI_xml_to_csv(tp.acc, blast_xml_tar, blast_csv_tar, tp.tmd_start, tp.tmd_end, 0.0001, logging)
    assert blast_csv_tar.is_file()
    if blast_csv_tar.parent.is_dir():
        rmtree(blast_csv_tar.parent)


def test_cleanup():
    # delete any unpacked xml or txt files in test directory
    test_inputs_dir: Path = Path(__file__).parents[1] / "test_inputs"
    files_to_delete = list(test_inputs_dir.glob("**/*.xml"))
    files_to_delete.extend(list(test_inputs_dir.glob("**/*.txt")))
    for file in files_to_delete:
        file.unlink()