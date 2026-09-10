"""The homologue source a prediction uses can be set per deployment.

A server chooses its homologue search with an environment variable, next to the other deployment
settings, rather than by editing a file inside the installed package. The shipped default stays
"ncbi", because that is what the golden-file tests pin and what reproduces the published
alignments, so an upgrade must not silently change an existing installation's output.
"""

import pytest
from thoipapy.predict import HOMOLOGUE_SOURCE_ENV_VAR, get_homologue_source


def test_the_settings_value_is_used_when_the_environment_is_silent(monkeypatch):
    monkeypatch.delenv(HOMOLOGUE_SOURCE_ENV_VAR, raising=False)
    assert get_homologue_source({"homologue_source": "ncbi"}) == "ncbi"


def test_the_shipped_default_is_ncbi_when_the_setting_is_absent(monkeypatch):
    monkeypatch.delenv(HOMOLOGUE_SOURCE_ENV_VAR, raising=False)
    assert get_homologue_source({}) == "ncbi"


def test_the_environment_overrides_the_settings_file(monkeypatch):
    monkeypatch.setenv(HOMOLOGUE_SOURCE_ENV_VAR, "colabfold")
    assert get_homologue_source({"homologue_source": "ncbi"}) == "colabfold"


def test_whitespace_and_an_empty_value_fall_back_rather_than_failing(monkeypatch):
    """An unset variable and one exported as empty should behave the same way."""
    for value in ("", "   "):
        monkeypatch.setenv(HOMOLOGUE_SOURCE_ENV_VAR, value)
        assert get_homologue_source({"homologue_source": "ncbi"}) == "ncbi"


def test_an_unrecognised_source_is_refused_rather_than_silently_ignored(monkeypatch):
    monkeypatch.setenv(HOMOLOGUE_SOURCE_ENV_VAR, "uniref90")
    with pytest.raises(ValueError, match="not one of"):
        get_homologue_source({})


def test_an_unrecognised_settings_value_is_refused_too(monkeypatch):
    monkeypatch.delenv(HOMOLOGUE_SOURCE_ENV_VAR, raising=False)
    with pytest.raises(ValueError, match="not one of"):
        get_homologue_source({"homologue_source": "blast"})
