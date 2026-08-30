"""Canonical paths for the thoipapy data directory.

The pipeline used to expect ``data_dir``, ``sets_dir`` and ``base_dir`` to be injected into the
settings dict by whatever called it. Neither shipped settings file defined them, so the
only thing that ever set them was the functional test, and the package could not be run
reproducibly from a clean checkout.

These constants replace that. ``data_dir`` now always resolves to the DVC-tracked ``data/``
directory of this repository, so a checkout plus ``dvc pull`` is enough to run the pipeline.

Deliberately NOT validated at import time. The webserver pip-installs thoipapy into
site-packages, where ``data/`` does not exist and is not meant to: the live prediction path
(``thoipapy.predict.run_THOIPA_prediction``) builds every path from its ``out_dir`` argument and
never reads any of these. Validating here would break that deployment. Call
:func:`assert_data_dir_available` instead at the top of anything that genuinely needs the
training data, so an installed copy fails loudly rather than silently writing into site-packages.
"""
from pathlib import Path

REPO_ROOT: Path = Path(__file__).parents[1]

# DVC-tracked. Contents come from `dvc pull`, not from the wheel.
DATA_DIR: Path = REPO_ROOT / "data"

# Protein set definitions (setNN.csv) and protein_names.csv. Small, human-authored inputs, so
# these live in git rather than DVC. Named "sets" because paper_figures/create_heatmap_from_
# merge_file.py looks for protein sets at base_dir/"sets" while run.py uses sets_dir; with
# base_dir == DATA_DIR, this single directory satisfies both conventions.
SETS_DIR: Path = DATA_DIR / "sets"

BASE_DIR: Path = DATA_DIR

SETTING_DIR: Path = Path(__file__).parent / "setting"

# The whitelist of features eligible for the ML model, and their groupings. Previously the
# "features" sheet of the settings workbook, which meant the settings file was simultaneously a
# key-value store and a data table.
MODEL_FEATURES_CSV: Path = SETTING_DIR / "model_features.csv"

# Which predictors appear in the comparison figures. Was the "selected_predictors" sheet of
# the settings workbook; it had no equivalent when the settings moved to CSV.
SELECTED_PREDICTORS_CSV: Path = SETTING_DIR / "selected_predictors.csv"

# Settings for a full training run, and for a single standalone prediction.
RUN_SETTINGS_CSV: Path = SETTING_DIR / "run_settings_example.csv"
STANDALONE_SETTINGS_CSV: Path = SETTING_DIR / "standalone_run_settings.csv"

ML_MODEL_DIR: Path = Path(__file__).parent / "ML_model"
TRAINED_MODEL_LPKL: Path = ML_MODEL_DIR / "THOIPA_trained_ML_model.lpkl"


def assert_data_dir_available() -> None:
    """Raise if the DVC-tracked training data is not present.

    Raises
    ------
    FileNotFoundError
        If ``data/`` is missing or empty, which means either `dvc pull` has not been run or
        thoipapy is being used from an installed copy that ships no training data.
    """
    if not DATA_DIR.is_dir() or not any(DATA_DIR.iterdir()):
        raise FileNotFoundError(
            f"The thoipapy data directory is not available at {DATA_DIR}.\n"
            "Run 'dvc pull' from the repository root to fetch it. Note that the training data is "
            "not shipped in the installed package, so pipeline code cannot be run from a pip "
            "install of thoipapy -- only from a checkout of the repository."
        )
