import warnings
from pathlib import Path

from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor, RandomForestClassifier
import os
import joblib
import pickle


def THOIPA_classifier_with_settings_deprecated(s, n_features, totally_randomized_trees=False):
    """ Classifier or Regressor settings are specified in the settings file.
    """
    # convert max_features to python None if "None"
    max_features = None if s["max_features"] == "None" else s["max_features"]
    # default for min_samples_leaf is 1 (no filter)
    min_samples_leaf = 1 if s["min_samples_leaf"] == "None" else s["min_samples_leaf"]
    # default for min_samples_leaf is None (no restriction on tree size)
    max_depth = None if s["max_depth"] == "None" else s["max_depth"]

    # forest = RandomForestClassifier(n_estimators=s["RF_number_of_estimators"], n_jobs=s["n_CPU_cores"], criterion=s["criterion"],
    #                                 min_samples_leaf=min_samples_leaf,
    #                                 max_depth=max_depth,
    #                                 oob_score=True, max_features=max_features, bootstrap=bool(s["bootstrap"]),
    #                                 #random_state=s["random_state"]
    #                                 )

    if totally_randomized_trees:
        max_features = 1
        min_samples_leaf = 1
        max_depth = None

    if s["bind_column"] == "interface":
        cls = ExtraTreesClassifier(n_estimators=s["RF_number_of_estimators"], n_jobs=s["n_CPU_cores"], criterion=s["criterion"],
                                        min_samples_leaf=min_samples_leaf,
                                        max_depth=max_depth,
                                        oob_score=bool(s["oob_estimation"]),
                                        bootstrap=bool(s["bootstrap"]),
                                        max_features=max_features,
                                        random_state=s["random_state"]
                                        )

    # in case you want to try to run the ML as a regressor
    elif s["bind_column"] == "interface_score_norm":
        cls = ExtraTreesRegressor(n_estimators=s["RF_number_of_estimators"], n_jobs=s["n_CPU_cores"],
                                        #criterion=s["criterion"],
                                        min_samples_leaf=min_samples_leaf,
                                        max_depth=max_depth,
                                        oob_score=bool(s["oob_estimation"]),
                                        bootstrap=bool(s["bootstrap"]),
                                        max_features=max_features,
                                        random_state=s["random_state"]
                                        )
    else:
        raise ValueError("bind column in excel settings file is not recognised ({})".format(s["bind_column"]))

    return cls


def train_machine_learning_model(s, logging):
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
    logging.info('starting train_machine_learning_model')

    # inputs
    train_data_after_first_feature_seln_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
    tuned_ensemble_parameters_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"
    # outputs
    model_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_ML_model.lpkl".format(s["setname"]))

    df_data = pd.read_csv(train_data_after_first_feature_seln_csv, index_col=0)

    if s["min_n_homol_training"] != 0:
        df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    df_data = df_data.dropna()

    cols_excluding_y = [c for c in df_data.columns if c != s['bind_column']]
    X = df_data[cols_excluding_y]
    y = df_data["interface"]

    if 1 not in y.tolist():
        raise ValueError("None of the residues are marked 1 for an interface residue!")

    n_features = X.shape[1]

    cls: ExtraTreesClassifier = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv)
    fit = cls.fit(X, y)
    joblib.dump(fit, model_pkl)

    tree_depths = np.array([estimator.tree_.max_depth for estimator in cls.estimators_])
    logging.info("tree depth mean = {} ({})".format(tree_depths.mean(), tree_depths))

    logging.info('finished training machine learning algorithm ({})'.format(model_pkl))


def return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv, totally_randomized_trees=False) -> ExtraTreesClassifier:
    df_tuned_ensemble_parameters: pd.DataFrame = pd.read_csv(tuned_ensemble_parameters_csv, index_col=0)
    ensemble_parameters_ser: pd.Series = df_tuned_ensemble_parameters["GridSearchSlowMethod"]

    n_estimators = int(ensemble_parameters_ser["n_estimators"])
    n_jobs = int(s["n_CPU_cores"])
    criterion = ensemble_parameters_ser["criterion"]
    oob_score = False

    if totally_randomized_trees:
        max_features = 1
        min_samples_leaf = 1
        max_depth = None
        bootstrap = False
    else:
        max_features = ensemble_parameters_ser["max_features"]
        min_samples_leaf = int(ensemble_parameters_ser["min_samples_leaf"])
        if pd.isnull(ensemble_parameters_ser["max_depth"]):
            max_depth = None
        else:
            max_depth = int(ensemble_parameters_ser["max_depth"])
        bootstrap = bool(s["bootstrap"])

    cls = ExtraTreesClassifier(
        n_estimators=n_estimators,
        n_jobs=n_jobs,
        criterion=criterion,
        min_samples_leaf=min_samples_leaf,
        max_depth=max_depth,
        oob_score=oob_score,
        bootstrap=bootstrap,
        max_features=max_features
    )
    return cls


