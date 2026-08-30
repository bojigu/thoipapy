import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier

from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import dropna_with_report

# Seed used for every ExtraTreesClassifier built by this module. The published 2020 model was
# trained with no seed at all, so it was not reproducible: retraining the same code on the same
# data moved the mean set07 ROC AUC over a range of roughly 0.61-0.66. Fixing the seed makes the
# shipped predictor reproducible from the DVC-tracked training data.
RANDOM_STATE = 0


def train_machine_learning_model(
    paths: ArtefactPaths, bind_column: str, min_n_homol_training: int, bootstrap: bool, logging
):
    """Train the machine learning model for a particular set.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Sawed Files
    -----------
    model_pkl : pickle
        Pickle containing the trained machine learning model.

    """
    logging.info("starting train_machine_learning_model")

    # inputs
    train_data_after_first_feature_seln_csv = paths.train_data_after_first_feature_seln_csv()
    tuned_ensemble_parameters_csv = paths.tuned_ensemble_parameters_csv()
    # outputs
    model_pkl = paths.ml_model_lpkl()

    df_data = pd.read_csv(train_data_after_first_feature_seln_csv, index_col=0)

    if min_n_homol_training != 0:
        df_data = df_data.loc[df_data.n_homologues >= min_n_homol_training]

    df_data = dropna_with_report(df_data, "train_machine_learning_model", logging)

    cols_excluding_y = [c for c in df_data.columns if c != bind_column]
    X = df_data[cols_excluding_y]
    y = df_data[bind_column]

    if 1 not in y.tolist():
        raise ValueError("None of the residues are marked 1 for an interface residue!")

    X.shape[1]

    cls: ExtraTreesClassifier = return_classifier_with_loaded_ensemble_parameters(
        tuned_ensemble_parameters_csv, bootstrap
    )
    fit = cls.fit(X, y)
    joblib.dump(fit, model_pkl)

    tree_depths = np.array([estimator.tree_.max_depth for estimator in cls.estimators_])
    logging.info(f"tree depth mean = {tree_depths.mean()} ({tree_depths})")

    logging.info(f"finished training machine learning algorithm ({model_pkl})")


def return_classifier_with_loaded_ensemble_parameters(
    tuned_ensemble_parameters_csv, bootstrap: bool, totally_randomized_trees=False, n_jobs=-1, random_state=RANDOM_STATE
) -> ExtraTreesClassifier:
    df_tuned_ensemble_parameters: pd.DataFrame = pd.read_csv(tuned_ensemble_parameters_csv, index_col=0)
    ensemble_parameters_ser: pd.Series = df_tuned_ensemble_parameters["GridSearchSlowMethod"]

    n_estimators = int(ensemble_parameters_ser["n_estimators"])
    criterion = ensemble_parameters_ser["criterion"]
    oob_score = False

    if totally_randomized_trees:
        max_features = 1
        min_samples_leaf = 1
        max_depth = None
        bootstrap = False
    else:
        max_features = ensemble_parameters_ser["max_features"]
        if max_features == "auto":
            # scikit-learn dropped max_features="auto" in 1.3. For classifiers it had always been
            # defined as sqrt(n_features), so this mapping preserves the tuned 2020 behaviour
            # exactly. The tuned parameter CSVs are DVC-tracked inputs and still carry "auto".
            max_features = "sqrt"
        min_samples_leaf = int(ensemble_parameters_ser["min_samples_leaf"])
        if pd.isnull(ensemble_parameters_ser["max_depth"]):
            max_depth = None
        else:
            max_depth = int(ensemble_parameters_ser["max_depth"])
        bootstrap = bool(bootstrap)

    cls = ExtraTreesClassifier(
        n_estimators=n_estimators,
        n_jobs=n_jobs,
        criterion=criterion,
        min_samples_leaf=min_samples_leaf,
        max_depth=max_depth,
        oob_score=oob_score,
        bootstrap=bootstrap,
        max_features=max_features,
        random_state=random_state,
    )
    return cls
