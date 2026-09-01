from pathlib import Path
from shutil import rmtree

from thoipapy.homologues.NCBI_parser import parse_NCBI_xml_to_csv
from thoipapy.utils import LogOnlyToConsole, make_sure_path_exists

from test.helpers.helpers import TestProtein


def test_parse_NCBI_xml_to_csv():
    # example xml tarball includes "CREATE_VIEW" string indicating malformed xml that needs cleaning, as often retrieved from NCBI
    blast_xml_tar: Path = (
        Path(__file__).parents[1] / "test_inputs/blast_data_with_malformed_xml/BLAST_results.xml.tar.gz"
    )
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


# Only tarballs are committed to these two directories; anything else in them was unpacked beside
# its archive during a test run.
DIRECTORIES_HOLDING_ONLY_TARBALLS = ["blast_data_valid", "blast_data_with_malformed_xml"]


def test_no_unpacked_archives_are_left_in_the_test_inputs():
    """Remove the xml and csv files the BLAST parsers unpack next to their tarballs.

    This used to glob "**/*.xml" and "**/*.txt" across the whole of test_inputs and unlink
    everything it found, so any committed fixture with one of those extensions was deleted by
    running the test suite -- and, since it ran as a test, the deletion looked like a pass.
    """
    test_inputs_dir: Path = Path(__file__).parents[1] / "test_inputs"

    for directory in DIRECTORIES_HOLDING_ONLY_TARBALLS:
        for unpacked in (test_inputs_dir / directory).iterdir():
            if unpacked.is_file() and not unpacked.name.endswith(".tar.gz"):
                unpacked.unlink()

    assert not list(test_inputs_dir.glob("**/BLAST_results.xml"))
