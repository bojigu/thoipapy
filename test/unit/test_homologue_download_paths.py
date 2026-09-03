"""The homologue search must create its output directory before writing into it.

``download_homologues_from_ncbi`` wrote its BLAST_details.txt two statements before calling
``make_sure_path_exists``, so it only worked where the directory already existed. Every directory
exists in a data directory restored by ``dvc pull``, which is why this was invisible; on a rebuild
from a clean checkout it stopped the pipeline at the first protein.
"""

import inspect
import tarfile

import pytest
from thoipapy.homologues import NCBI_download


def test_the_output_directory_is_created_before_either_file_is_written():
    """Both files land in the same directory, so both mkdir calls must precede both writes."""
    source = inspect.getsource(NCBI_download.download_homologues_from_ncbi)
    lines = source.splitlines()

    first_write = next(n for n, line in enumerate(lines) if "open(xml_txt" in line)
    mkdir_lines = [n for n, line in enumerate(lines) if "make_sure_path_exists" in line]

    assert mkdir_lines, "download_homologues_from_ncbi no longer creates its output directory"
    assert max(mkdir_lines) < first_write, (
        "make_sure_path_exists runs after the details file is written, so the search only "
        "works where the directory already exists"
    )


def test_a_search_into_a_missing_directory_reaches_the_blast_call(tmp_path, monkeypatch):
    """End to end for the failure that was reported: nothing on disk but the data root."""
    missing = tmp_path / "homologues" / "xml" / "crystal"
    assert not missing.exists()

    monkeypatch.setenv(NCBI_download.LOCAL_BLAST_DB_ENV_VAR, str(tmp_path / "somedb"))
    called = {}

    def fake_blastp(query_fasta_string, blast_xml_file, expect_value, hit_list_size, db, logging):
        called["db"] = db
        # The real blastp writes this; the parser downstream needs it to exist.
        open(blast_xml_file, "w").write("<BlastOutput></BlastOutput>")

    monkeypatch.setattr(NCBI_download, "run_local_blastp", fake_blastp)

    class Logger:
        def info(self, *a, **k):
            pass

        def warning(self, *a, **k):
            pass

        def exception(self, *a, **k):
            raise AssertionError("the search raised")

    NCBI_download.download_homologues_from_ncbi(
        acc="1xioA4",
        TMD_seq_pl_surr="GKVIGIAVMLTGISALTLLIGTVSNMFQ",
        blast_xml_file=missing / "1xioA4.surr20.BLAST.xml",
        xml_txt=missing / "1xioA4.surr20.BLAST_details.txt",
        xml_tar_gz=missing / "1xioA4.surr20.BLAST.xml.tar.gz",
        expect_value=10,
        hit_list_size=10,
        logging=Logger(),
    )

    assert called["db"] == str(tmp_path / "somedb")

    # The xml and its details file are archived into the tar.gz and then deleted, so the
    # provenance has to be read back out of the archive. That record is the only way to tell a
    # 2020 nr download from a current UniRef90 search after the fact.
    archive = missing / "1xioA4.surr20.BLAST.xml.tar.gz"
    assert archive.is_file()
    with tarfile.open(archive) as tar:
        names = tar.getnames()
        assert "1xioA4.surr20.BLAST_details.txt" in names
        details = tar.extractfile("1xioA4.surr20.BLAST_details.txt").read().decode()
    assert "1xioA4" in details
    assert "database" in details and "download_date" in details
    assert str(tmp_path / "somedb") in details


def test_the_details_file_records_the_database_that_was_searched(tmp_path, monkeypatch):
    monkeypatch.delenv(NCBI_download.LOCAL_BLAST_DB_ENV_VAR, raising=False)
    assert NCBI_download.get_local_blast_db() == ""

    monkeypatch.setenv(NCBI_download.LOCAL_BLAST_DB_ENV_VAR, "/blastdb/uniref90")
    assert NCBI_download.get_local_blast_db() == "/blastdb/uniref90"


@pytest.mark.parametrize("value,expected", [("  ", ""), ("/db/uniref90  ", "/db/uniref90")])
def test_a_blank_local_database_setting_means_query_ncbi(value, expected, monkeypatch):
    monkeypatch.setenv(NCBI_download.LOCAL_BLAST_DB_ENV_VAR, value)
    assert NCBI_download.get_local_blast_db() == expected
