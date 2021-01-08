import sys
from pathlib import Path

import numpy as np
import pandas as pd

from thoipapy.ML_model.tune import tune_ensemble_parameters
from thoipapy.utils import make_sure_path_exists
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters


def get_initial_ensemble_parameters_before_feature_selection(s, logging):
    # logging.info('RF_variable_importance_calculate is running\n')
    train_data_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_train_data_orig.csv"
    # output
    tuned_ensemble_parameters_before_feature_seln_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_tuned_ensemble_parameters_before_feature_seln.csv"

    tune_ensemble_parameters(s, train_data_csv, tuned_ensemble_parameters_before_feature_seln_csv, logging)


def calc_feat_import_using_MDI_before_feature_seln(s, logging):
    """Calculate the variable importance (mean decrease gini) for all variables in THOIPA.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    mean_decrease_impurity_all_features_csv : csv
        List of variables, sorted by their importance to the algorithm.
        Also includes the standard deviation supplied by the machine learning algorithm
    """
    # logging.info('RF_variable_importance_calculate is running\n')
    train_data_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_train_data_orig.csv"
    tuned_ensemble_parameters_before_feature_seln_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_tuned_ensemble_parameters_before_feature_seln.csv"

    # mean_decrease_impurity_all_features_csv = Path(s["data_dir"]) / "results" / s["setname"] / "feat_imp/mean_decrease_impurity_all_features.csv"
    mean_decrease_impurity_all_features_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_feat_imp_MDI_before_feature_seln.csv"

    make_sure_path_exists(mean_decrease_impurity_all_features_csv, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    X = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])
    y = df_data["interface"]
    n_features = X.shape[1]

    # regular trees, or Totally Randomized Trees
    model_types = ["", "_TRT"]
    output_dfs = []

    for model_type in model_types:
        if model_type == "_TRT":
            forest = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_before_feature_seln_csv, totally_randomized_trees=True)
            logging.info("IMPORTANCES FOR TOTALLY RANDOMIZED TREES (max_features=1, max_depth=None, min_samples_leaf=1)")
        elif model_type == "":
            forest = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_before_feature_seln_csv)
            logging.info("Feature ranking:")
        else:
            raise ValueError("model type unknown")

        df_imp = calculate_mean_decrease_impurity_for_dataset(X, y, forest, model_type, logging)

        output_dfs.append(df_imp)

    df_imp2 = pd.concat(output_dfs, axis=1)
    df_imp2.sort_values("order_importance", inplace=True)
    df_imp2["original_order"] = df_imp2.index
    df_imp2.set_index("feature", inplace=True)

    df_imp2.to_csv(mean_decrease_impurity_all_features_csv)


def calculate_mean_decrease_impurity_for_dataset(X, y, forest, model_type, logging):
    forest.fit(X, y)
    importances_arr = forest.feature_importances_
    std_arr = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    indices_arr = np.argsort(importances_arr)[::-1]
    importances_text_list = X.columns.tolist()
    # order_list = [importances_arr[indices_arr[f]] for f in range(X.shape[1])]
    nested_dict = {}
    for f in range(X.shape[1]):
        if f < 10:
            logging.info("%d. feature %d (%f) %s" % (f + 1, indices_arr[f], importances_arr[indices_arr[f]], importances_text_list[indices_arr[f]]))
        single_feature_dict = {"original_order": indices_arr[f], "mean_decrease_impurity{}".format(model_type): importances_arr[indices_arr[f]], "feature{}".format(model_type): importances_text_list[indices_arr[f]],
                               "std{}".format(model_type): std_arr[f]}
        nested_dict[f + 1] = single_feature_dict
    sys.stdout.write("\n\n"), sys.stdout.flush()
    df_imp = pd.DataFrame(nested_dict).T
    df_imp["order_importance{}".format(model_type)] = df_imp.index
    # df_imp.set_index("feature", inplace=True)
    df_imp.set_index("original_order", inplace=True)
    return df_imp
