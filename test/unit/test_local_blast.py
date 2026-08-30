"""The local BLAST path, which lets a prediction run without querying NCBI.

NCBI queues automated BLAST clients, and the queue is sometimes hours deep: a probe of a
three-residue query returned RTOE=11165, an estimate of 3.1 hours. A webserver job with a
one-hour ceiling cannot survive that, and a test suite certainly cannot depend on it.

Pointing THOIPA_LOCAL_BLAST_DB at a local BLAST database makes the predictor search that
instead. Swiss-Prot is ~340 MB built and answers in well under a second.
"""

import os

import pytest
from thoipapy.homologues.NCBI_download import (
    LOCAL_BLAST_DB_ENV_VAR,
    get_local_blast_db,
    run_local_blastp,
)


def test_no_local_db_configured_means_ncbi(monkeypatch):
    monkeypatch.delenv(LOCAL_BLAST_DB_ENV_VAR, raising=False)
    assert get_local_blast_db() == ""


def test_configured_local_db_is_returned(monkeypatch):
    monkeypatch.setenv(LOCAL_BLAST_DB_ENV_VAR, "/some/path/swissprot")
    assert get_local_blast_db() == "/some/path/swissprot"


def test_whitespace_only_setting_is_treated_as_unset(monkeypatch):
    monkeypatch.setenv(LOCAL_BLAST_DB_ENV_VAR, "   ")
    assert get_local_blast_db() == ""


@pytest.mark.requires_tool("blastp")
def test_local_blastp_writes_parseable_xml(tmp_path):
    """The output must be BLAST XML, because NCBI_parser reads exactly that.

    Skipped unless a local database is configured, since it needs one to search.
    """
    db = os.environ.get(LOCAL_BLAST_DB_ENV_VAR, "").strip()
    if not db:
        pytest.skip(f"{LOCAL_BLAST_DB_ENV_VAR} not set; no local BLAST database to search")

    import logging

    xml = tmp_path / "BLAST_results.xml"
    query = ">test\nPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSP"
    run_local_blastp(query, xml, expect_value=10, hit_list_size=100, db=db, logging=logging)

    assert xml.is_file(), "blastp produced no output file"
    text = xml.read_text()
    assert "<BlastOutput>" in text, "output is not BLAST XML; NCBI_parser cannot read it"
    assert "<Hit>" in text, "no hits returned for a glycophorin A query"

    # The query fasta is a working file and must not be left behind.
    assert not xml.with_suffix(".query.fasta").exists()


@pytest.mark.requires_tool("blastp")
def test_local_blastp_raises_on_a_missing_database(tmp_path):
    """A missing database must fail loudly, not leave an empty file that looks like no hits."""
    import logging

    with pytest.raises(RuntimeError, match="local blastp failed"):
        run_local_blastp(
            ">t\nMKV",
            tmp_path / "out.xml",
            expect_value=10,
            hit_list_size=10,
            db=str(tmp_path / "does_not_exist"),
            logging=logging,
        )
