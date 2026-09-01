"""Tests for the reduction of a user-supplied protein name to something safe in a filename."""

import os

from thoipapy.utils import MAX_SAFE_NAME_LEN, safe_filename_component

# The longest suffix thoipapy appends to an accession when building an intermediate filename.
LONGEST_THOIPAPY_SUFFIX = ".surr20.gaps5.uniq.for_PSSM_FREECONTACT.txt"


def test_a_uniprot_fasta_header_survives_as_something_readable():
    """The single most likely thing a user pastes."""
    header = "sp|Q8WWF3|SSMM1_HUMAN Serine-rich single-pass membrane protein 1 OS=Homo sapiens OX=9606"
    assert safe_filename_component(header) == "sp_Q8WWF3_SSMM1_HUMAN_Serine-rich_single-pass_membrane_prote"


def test_shell_metacharacters_are_removed():
    assert safe_filename_component("x; touch /tmp/EXECUTED; echo") == "x_touch_tmp_EXECUTED_echo"
    assert safe_filename_component("`whoami`") == "whoami"
    assert safe_filename_component("$(id)") == "id"


def test_path_separators_and_traversal_are_removed():
    """A name is a filename component, never a path."""
    result = safe_filename_component("../../etc/passwd")
    assert "/" not in result
    assert not result.startswith(".")
    assert result == "etc_passwd"


def test_a_long_name_leaves_room_for_the_suffixes_thoipapy_appends():
    """NAME_MAX is 255; a full FASTA header plus a suffix used to risk ENAMETOOLONG."""
    result = safe_filename_component("A" * 500)
    assert len(result) == MAX_SAFE_NAME_LEN
    assert len(result + LONGEST_THOIPAPY_SUFFIX) < os.pathconf("/", "PC_NAME_MAX")


def test_a_name_with_nothing_usable_in_it_still_produces_a_filename():
    """An empty component would put the suffix at the start of the name, and collide."""
    assert safe_filename_component("...") == "protein"
    assert safe_filename_component("") == "protein"
    assert safe_filename_component("###") == "protein"


def test_case_and_the_separators_users_rely_on_are_kept():
    """Unlike the slugify it replaces, which lowercased and dropped dots and underscores."""
    assert safe_filename_component("P21860_ERBB3") == "P21860_ERBB3"
    assert safe_filename_component("1xioA4_rhodopsin") == "1xioA4_rhodopsin"
    assert safe_filename_component("Q9Y5-2.1") == "Q9Y5-2.1"


def test_a_name_that_has_already_been_sanitised_is_unchanged():
    """The THOIPA webserver applies the same rule at its own boundary."""
    already_safe = "sp_Q8WWF3_SSMM1_HUMAN_Serine-rich_single-pass_membrane_prote"
    assert safe_filename_component(already_safe) == already_safe
