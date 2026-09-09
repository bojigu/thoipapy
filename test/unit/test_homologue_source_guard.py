"""Tests for the homologue-source / stage coupling guard."""

import pytest
from thoipapy.run import VALID_HOMOLOGUE_SOURCES, check_homologue_source_matches_stages
from thoipapy.run_settings import RunSettings


def test_shipped_defaults_are_accepted():
    check_homologue_source_matches_stages({"homologue_source": "ncbi"}, RunSettings())


def test_a_missing_homologue_source_defaults_to_ncbi():
    check_homologue_source_matches_stages({}, RunSettings())


def test_colabfold_stages_with_the_ncbi_source_are_refused():
    """This is the case that would overwrite the DVC-tracked 2020 nr artefacts."""
    for stages in (
        RunSettings(run_retrieve_homologues_from_colabfold=True),
        RunSettings(run_parse_colabfold_a3m_into_csv=True),
    ):
        with pytest.raises(ValueError, match="ColabFold homologue stage is enabled"):
            check_homologue_source_matches_stages({"homologue_source": "ncbi"}, stages)


def test_ncbi_stages_with_the_colabfold_source_are_refused():
    for stages in (
        RunSettings(run_retrieve_NCBI_homologues_with_blastp=True),
        RunSettings(run_parse_homologues_xml_into_csv=True),
    ):
        with pytest.raises(ValueError, match="NCBI homologue stage is enabled"):
            check_homologue_source_matches_stages({"homologue_source": "colabfold"}, stages)


def test_matching_source_and_stages_are_accepted():
    check_homologue_source_matches_stages(
        {"homologue_source": "colabfold"}, RunSettings(run_parse_colabfold_a3m_into_csv=True)
    )
    check_homologue_source_matches_stages(
        {"homologue_source": "ncbi"}, RunSettings(run_parse_homologues_xml_into_csv=True)
    )


def test_an_unrecognised_source_is_refused_rather_than_creating_a_directory():
    with pytest.raises(ValueError, match="not one of"):
        check_homologue_source_matches_stages({"homologue_source": "uniref90"}, RunSettings())


def test_the_valid_sources_are_the_two_directories_that_exist():
    assert VALID_HOMOLOGUE_SOURCES == ("ncbi", "colabfold")
