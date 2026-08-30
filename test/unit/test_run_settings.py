"""Tests for the typed pipeline stage flags."""
from dataclasses import fields

import pytest

from thoipapy.common import create_settingdict
from thoipapy.paths import RUN_SETTINGS_CSV, STANDALONE_SETTINGS_CSV
from thoipapy.run_settings import RunSettings


def test_missing_flags_default_to_off():
    """A flag absent from the spreadsheet must mean "do not run", not KeyError mid-pipeline."""
    stages = RunSettings.from_settings({})
    assert stages.enabled_stages() == []
    assert stages.train_machine_learning_model is False


@pytest.mark.parametrize("raw,expected", [
    (True, True), (False, False),
    (1, True), (0, False),
    ("TRUE", True), ("FALSE", False),
    ("true", True), ("false", False),
    ("T", True), ("F", False),
    ("Yes", True), ("no", False),
])
def test_boolean_coercion(raw, expected):
    """The spreadsheet's type column is unreliable, so values are coerced here instead.

    The FALSE case is the one that matters: a non-empty string is truthy in python, so before
    this an explicitly disabled stage could still run.
    """
    stages = RunSettings.from_settings({"run_validation": raw})
    assert stages.run_validation is expected


def test_string_false_does_not_enable_a_stage():
    """Regression guard for the specific bug: "FALSE" is truthy as a bare string."""
    assert RunSettings.from_settings({"run_validation": "FALSE"}).run_validation is False


@pytest.mark.parametrize("bad", ["maybe", "", "2", None, 3.7, []])
def test_unrecognised_values_are_rejected_at_load_time(bad):
    with pytest.raises(ValueError, match="not recognisably true or false"):
        RunSettings.from_settings({"pssm_calculation": bad})


def test_is_frozen():
    stages = RunSettings.from_settings({})
    with pytest.raises(Exception):
        stages.pssm_calculation = True


def test_unknown_keys_in_settings_are_ignored():
    """The settings dict carries far more than stage flags; the rest must not break construction."""
    stages = RunSettings.from_settings({"pssm_calculation": True, "lipophilicity_scale": "Engelman(GES)"})
    assert stages.pssm_calculation is True


@pytest.mark.parametrize("settings_csv", [STANDALONE_SETTINGS_CSV, RUN_SETTINGS_CSV])
def test_builds_from_the_shipped_settings_files(settings_csv):
    """Both shipped settings files must produce a valid RunSettings, with no unreadable flag."""
    stages = RunSettings.from_settings(create_settingdict(settings_csv))
    assert isinstance(stages.enabled_stages(), list)
    # every field is a real bool, not a string that merely looks like one
    for field in fields(stages):
        assert isinstance(getattr(stages, field.name), bool), f"{field.name} is not a bool"


def test_standalone_settings_enable_the_feature_pipeline():
    """Pins what the shipped standalone settings actually switch on."""
    stages = RunSettings.from_settings(create_settingdict(STANDALONE_SETTINGS_CSV))
    for stage in ["pssm_calculation", "entropy_calculation", "rate4site_calculation",
                  "coevolution_calculation", "lips_score_calculation", "motifs_from_seq"]:
        assert getattr(stages, stage) is True, f"{stage} should be enabled"
    # validation is off, and known to be broken
    assert stages.run_validation is False
    assert stages.run_testset_trainset_validation is False
