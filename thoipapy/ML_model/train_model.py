from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
import os
import joblib
import pickle


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
    train_data_after_first_feature_seln_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
    tuned_ensemble_parameters_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"
    # outputs
    model_pkl = os.path.join(s["data_dir"], "results", s["setname"], "{}_ML_model.lpkl".format(s["setname"]))

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


def return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv, totally_randomized_trees=False, n_jobs=-1) -> ExtraTreesClassifier:
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
