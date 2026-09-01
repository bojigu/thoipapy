"""Fixtures shared by the unit tests."""

from pathlib import Path

import pytest

# Only tarballs are committed to these directories. Anything else in them was unpacked beside its
# archive by a parser under test.
DIRECTORIES_HOLDING_ONLY_TARBALLS = ["blast_data_valid", "blast_data_with_malformed_xml"]

TEST_INPUTS_DIR = Path(__file__).parents[1] / "test_inputs"


@pytest.fixture(autouse=True, scope="session")
def remove_unpacked_archives_from_the_test_inputs():
    """Delete the files the BLAST parsers unpack next to their tarballs, after the session.

    This was a function named test_cleanup that globbed "**/*.xml" and "**/*.txt" across the whole
    of test_inputs and unlinked everything it found, so running the suite deleted committed
    fixtures -- and, because it ran as a test, the deletion was reported as a pass. It also only
    cleaned up if it happened to run last, which was guaranteed by nothing but declaration order.
    """
    yield
    for directory in DIRECTORIES_HOLDING_ONLY_TARBALLS:
        unpacked_into = TEST_INPUTS_DIR / directory
        if not unpacked_into.is_dir():
            continue
        for unpacked in unpacked_into.iterdir():
            if unpacked.is_file() and not unpacked.name.endswith(".tar.gz"):
                unpacked.unlink()
