"""NCBI requires automated BLAST clients to identify themselves.

Clients that do not are throttled or refused: repeated identical queries from an anonymous client
get pushed into a slow queue, and a run that took 390 seconds on first submission took roughly 15
minutes once NCBI had seen it several times. Biopython's qblast has no email argument, so this is
set through the module-level NCBIWWW.email that qblast forwards as a URL parameter.
"""

import os

import pytest
from Bio.Blast import NCBIWWW

from thoipapy.homologues.NCBI_download import NCBI_CONTACT_EMAIL_ENV_VAR, identify_to_ncbi


@pytest.fixture(autouse=True)
def restore_ncbiwww_state():
    """NCBIWWW.email is module-level global state; put it back after each test."""
    email, tool = NCBIWWW.email, NCBIWWW.tool
    yield
    NCBIWWW.email, NCBIWWW.tool = email, tool


def test_contact_email_is_forwarded_to_biopython(monkeypatch):
    monkeypatch.setenv(NCBI_CONTACT_EMAIL_ENV_VAR, "someone@example.com")
    identify_to_ncbi()
    assert NCBIWWW.email == "someone@example.com"
    assert NCBIWWW.tool == "thoipapy"


def test_missing_contact_email_raises_rather_than_querying_anonymously(monkeypatch):
    """Failing loudly is the point: an anonymous query is what gets the server throttled."""
    monkeypatch.delenv(NCBI_CONTACT_EMAIL_ENV_VAR, raising=False)
    with pytest.raises(ValueError, match=NCBI_CONTACT_EMAIL_ENV_VAR):
        identify_to_ncbi()


def test_whitespace_only_contact_email_is_treated_as_unset(monkeypatch):
    monkeypatch.setenv(NCBI_CONTACT_EMAIL_ENV_VAR, "   ")
    with pytest.raises(ValueError):
        identify_to_ncbi()


def test_qblast_actually_forwards_the_module_level_email():
    """Guard the assumption this design rests on, in case Biopython changes it."""
    import inspect
    source = inspect.getsource(NCBIWWW.qblast)
    assert '"email": email' in source, (
        "Biopython's qblast no longer forwards the module-level NCBIWWW.email. "
        "identify_to_ncbi needs another route to identify this client to NCBI."
    )
