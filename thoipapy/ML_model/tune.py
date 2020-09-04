from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV

from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML


def tune_ensemble_parameters_after_feature_seln(s, logging):
    """Tune the machine learning algorithm to obtain the best ensemble parameters.

    A broader quick search for parameters is first made with RandomizedSearchCV.
        The results are printed and saved, but not automatically used in the pipeline.
        These are designed to guide the manual setting of the range for the GridSearch parameters.
    A more narrow, slow search is made using GridSearch
        The best parameters are automatically used to create the ML model used for
        other predictions (e.g. testset/trainset predictions)

    """

    logging.info('starting tune_ensemble_parameters')
    # inputs
    train_data_after_first_feature_seln_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
    # output
    tuned_ensemble_parameters_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"

    tune_ensemble_parameters(s, train_data_after_first_feature_seln_csv, tuned_ensemble_parameters_csv, logging)


def tune_ensemble_parameters(s, train_data_csv, tuned_ensemble_parameters_csv, logging):
    df_data = pd.read_csv(train_data_csv, index_col=0)
    if True in df_data.columns.str.contains("Unnamed").tolist():
        raise ValueError(f"unnamed column found when reading {train_data_csv}")
    y = df_data[s["bind_column"]]
    X = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])
    #del X[s["bind_column"]]
    assert "interface" not in X.columns
    cls = ExtraTreesClassifier(oob_score=False)

    best_params_rs = get_optimised_ensemble_parameters_using_quick_randomizedsearchcv_method(X, cls, logging, y)

    best_params_gs = get_optimised_ensemble_parameters_using_slow_gridsearchcv_method(X, cls, logging, y)

    df_tuned_ensemble_parameters = pd.DataFrame()
    df_tuned_ensemble_parameters["RandomizedSearchQuickMethod"] = best_params_rs
    df_tuned_ensemble_parameters["GridSearchSlowMethod"] = best_params_gs
    df_tuned_ensemble_parameters.to_csv(tuned_ensemble_parameters_csv)
    logging.info(f"parameters tuned and saved ({tuned_ensemble_parameters_csv})")


def get_optimised_ensemble_parameters_using_slow_gridsearchcv_method(X, cls, logging, y):
    logging.info("Calculating best parameters using slow GridSearch method")
    ############# Parameter ranges for slow GridSearchCV method ###########
    criterion = ["gini", "entropy"]
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(50, 200, 6)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = list(range(5, 10)) + [None]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = range(3, 7)
    # Method of selecting samples for training each tree
    bootstrap = [False]
    # Create the random grid
    random_grid_gs = {'criterion': criterion,
                      'n_estimators': n_estimators,
                      'max_features': max_features,
                      'max_depth': max_depth,
                      'min_samples_leaf': min_samples_leaf,
                      'bootstrap': bootstrap}
    cls_gridsearch = GridSearchCV(cls, random_grid_gs, cv=5, verbose=1, n_jobs=-1, scoring="average_precision")
    cls_gridsearch.fit(X, y)
    logging.info(f"best parameters: {cls_gridsearch.best_params_}")
    logging.info(f"best score: {cls_gridsearch.best_score_}")
    logging.info("Use parameters from RandomizedSearch to tweak range for GridSearch if necessary.")
    best_params_gs = pd.Series(cls_gridsearch.best_params_)
    best_params_gs["best_score"] = cls_gridsearch.best_score_
    return best_params_gs


def get_optimised_ensemble_parameters_using_quick_randomizedsearchcv_method(X, cls, logging, y):
    logging.info("Calculating best parameters using quick RandomizedSearchCV method")
    ############# Parameter ranges for quick RandomizedSearchCV method ###########
    criterion = ["gini", "entropy"]
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(20, 200, 19)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = list(range(2, 10)) + [None]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = range(3, 10)
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid_rs_cv = {'criterion': criterion,
                         'n_estimators': n_estimators,
                         'max_features': max_features,
                         'max_depth': max_depth,
                         'min_samples_leaf': min_samples_leaf,
                         'bootstrap': bootstrap}
    cls_random = RandomizedSearchCV(estimator=cls, param_distributions=random_grid_rs_cv, n_iter=100, cv=5, verbose=1, n_jobs=-1, scoring="average_precision")
    cls_random.fit(X, y)
    logging.info(f"best parameters: {cls_random.best_params_}")
    logging.info(f"best score: {cls_random.best_score_}")
    best_params_rs = pd.Series(cls_random.best_params_)
    best_params_rs["best_score"] = cls_random.best_score_
    logging.info('finished tune_ensemble_parameters')
    return best_params_rs

