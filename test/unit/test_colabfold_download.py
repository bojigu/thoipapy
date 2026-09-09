"""Unit tests for the ColabFold MSA server client.

The retry loop and the download-and-validate step had no coverage, and both carried a real defect:
an unbounded retry that could hang a whole protein set, and a failed download that left a corrupt
archive the resume logic then treated as complete. Everything here runs against a fake transport,
so no test touches the network.
"""

import tarfile
from pathlib import Path

import pytest
import requests
from thoipapy.homologues import colabfold_download
from thoipapy.homologues.colabfold_download import (
    MAX_HTTP_RETRIES,
    _post_json_with_retries,
    _request_with_retries,
    mode_searches_env_db,
    validate_mode,
)
from thoipapy.utils import LogOnlyToConsole


class _Response:
    """Enough of requests.Response for these tests."""

    def __init__(self, payload=None, content=b"", status_code=200, text=""):
        self._payload = payload
        self.content = content
        self.status_code = status_code
        self.text = text

    def json(self):
        if self._payload is None:
            raise ValueError("no json")
        return self._payload


@pytest.fixture(autouse=True)
def _no_sleeping(monkeypatch):
    """Retries back off; the tests should not actually wait."""
    monkeypatch.setattr(colabfold_download.time, "sleep", lambda _seconds: None)


def test_a_connect_timeout_does_not_retry_forever(monkeypatch):
    """ConnectTimeout subclasses both ConnectionError and Timeout.

    It used to reach a Timeout branch that neither counted the attempt nor slept, so a host that
    accepted connections and never answered would spin at one attempt per timeout indefinitely.
    """
    attempts = []

    def always_times_out(*_args, **_kwargs):
        attempts.append(1)
        raise requests.exceptions.ConnectTimeout("connect timed out")

    monkeypatch.setattr(colabfold_download.requests, "request", always_times_out)

    with pytest.raises(RuntimeError, match="unreachable"):
        _request_with_retries("GET", "https://example.invalid/x", {}, LogOnlyToConsole())

    assert len(attempts) == MAX_HTTP_RETRIES


def test_a_read_timeout_is_also_bounded(monkeypatch):
    attempts = []

    def always_times_out(*_args, **_kwargs):
        attempts.append(1)
        raise requests.exceptions.ReadTimeout("read timed out")

    monkeypatch.setattr(colabfold_download.requests, "request", always_times_out)

    with pytest.raises(RuntimeError):
        _request_with_retries("GET", "https://example.invalid/x", {}, LogOnlyToConsole())

    assert len(attempts) == MAX_HTTP_RETRIES


def test_a_transient_failure_is_retried_and_then_succeeds(monkeypatch):
    calls = []

    def fail_once_then_work(*_args, **_kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise requests.exceptions.ConnectionError("boom")
        return _Response(payload={"status": "COMPLETE"})

    monkeypatch.setattr(colabfold_download.requests, "request", fail_once_then_work)

    result = _request_with_retries("GET", "https://example.invalid/x", {}, LogOnlyToConsole())
    assert result.json() == {"status": "COMPLETE"}
    assert len(calls) == 2


def test_an_html_error_page_is_retried_rather_than_aborting_the_set(monkeypatch):
    """A 502 from a load balancer arrives as a successful response with the wrong body."""
    calls = []

    def gateway_error_then_json(*_args, **_kwargs):
        calls.append(1)
        if len(calls) == 1:
            return _Response(payload=None, status_code=502, text="<html>502 Bad Gateway</html>")
        return _Response(payload={"id": "abc", "status": "PENDING"})

    monkeypatch.setattr(colabfold_download.requests, "request", gateway_error_then_json)

    out = _post_json_with_retries("https://example.invalid/ticket/msa", {}, {}, LogOnlyToConsole())
    assert out["id"] == "abc"
    assert len(calls) == 2


def test_persistent_non_json_eventually_raises(monkeypatch):
    monkeypatch.setattr(
        colabfold_download.requests,
        "request",
        lambda *_a, **_k: _Response(payload=None, status_code=502, text="<html>nope</html>"),
    )
    with pytest.raises(RuntimeError, match="did not reply with json"):
        _post_json_with_retries("https://example.invalid/ticket/msa", {}, {}, LogOnlyToConsole())


def _fake_server(monkeypatch, archive_bytes):
    """Drive submit -> poll -> download without a network."""
    monkeypatch.setattr(colabfold_download, "get_user_agent", lambda: "thoipapy/test test@example.com")

    def fake_request(method, url, *_args, **_kwargs):
        if method == "POST":
            return _Response(payload={"id": "TICKET", "status": "COMPLETE"})
        if "/ticket/" in url:
            return _Response(payload={"status": "COMPLETE"})
        return _Response(content=archive_bytes)

    monkeypatch.setattr(colabfold_download.requests, "request", fake_request)


def _valid_archive(tmp_path: Path) -> bytes:
    a3m = tmp_path / "uniref.a3m"
    a3m.write_text(">101\nAAAA\n")
    archive = tmp_path / "ok.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(a3m, arcname="uniref.a3m")
    return archive.read_bytes()


def test_a_successful_download_writes_the_archive_and_details(tmp_path, monkeypatch):
    _fake_server(monkeypatch, _valid_archive(tmp_path))
    out = tmp_path / "out" / "acc.surr20.env.a3m.tar.gz"
    details = tmp_path / "out" / "acc.surr20.env_details.txt"

    colabfold_download.download_homologues_from_colabfold("acc", "AAAA", out, details, "env", LogOnlyToConsole())

    assert out.is_file()
    assert "colabfold_mode\tenv" in details.read_text()


def test_a_corrupt_download_leaves_nothing_the_resume_logic_would_accept(tmp_path, monkeypatch):
    """The skip-if-exists check must never see a file that is not a usable archive.

    Writing straight to the final path meant an HTML error body or a truncated transfer left a
    file that later runs skipped, so the protein could never repair itself.
    """
    _fake_server(monkeypatch, b"<html>not a tarball</html>")
    out = tmp_path / "out" / "acc.surr20.env.a3m.tar.gz"
    details = tmp_path / "out" / "acc.surr20.env_details.txt"

    with pytest.raises(RuntimeError, match="not a tar archive"):
        colabfold_download.download_homologues_from_colabfold("acc", "AAAA", out, details, "env", LogOnlyToConsole())

    assert not out.exists()
    assert not out.with_suffix(out.suffix + ".partial").exists()
    assert not details.exists()


def test_an_archive_without_uniref_is_rejected_and_not_left_behind(tmp_path, monkeypatch):
    other = tmp_path / "other.a3m"
    other.write_text(">101\nAAAA\n")
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(other, arcname="something_else.a3m")

    _fake_server(monkeypatch, archive.read_bytes())
    out = tmp_path / "out" / "acc.surr20.env.a3m.tar.gz"

    with pytest.raises(RuntimeError, match="uniref.a3m"):
        colabfold_download.download_homologues_from_colabfold(
            "acc", "AAAA", out, tmp_path / "out" / "d.txt", "env", LogOnlyToConsole()
        )

    assert not out.exists()


def test_an_unknown_mode_never_reaches_the_server(tmp_path, monkeypatch):
    def explode(*_args, **_kwargs):
        raise AssertionError("the server must not be contacted for an invalid mode")

    monkeypatch.setattr(colabfold_download.requests, "request", explode)

    with pytest.raises(ValueError, match="not one of"):
        colabfold_download.download_homologues_from_colabfold(
            "acc", "AAAA", tmp_path / "a.a3m.tar.gz", tmp_path / "d.txt", "environmental", LogOnlyToConsole()
        )


def test_mode_helpers():
    assert mode_searches_env_db("env") is True
    assert mode_searches_env_db("env-nofilter") is True
    assert mode_searches_env_db("all") is False
    assert mode_searches_env_db("nofilter") is False
    with pytest.raises(ValueError):
        validate_mode("nonsense")
