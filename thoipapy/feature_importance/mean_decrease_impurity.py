import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from thoipapy.utils import create_colour_lists, make_sure_path_exists
from thoipapy.feature_importance.plots import create_var_imp_plot
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML
from thoipapy.validation.train_model import THOIPA_classifier_with_settings


def calc_feat_import_from_mean_decrease_impurity(s, logging):
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
    #logging.info('RF_variable_importance_calculate is running\n')
    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))

    mean_decrease_impurity_all_features_csv = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / "feat_imp/mean_decrease_impurity_all_features.csv"
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
            forest = THOIPA_classifier_with_settings(s, n_features, totally_randomized_trees=True)
            logging.info("IMPORTANCES FOR TOTALLY RANDOMIZED TREES (max_features=1, max_depth=None, min_samples_leaf=1)")
        elif model_type == "":
            forest = THOIPA_classifier_with_settings(s, n_features)
            logging.info("Feature ranking:")
        else:
            raise ValueError("model type unknown")
        forest.fit(X, y)
        importances_arr = forest.feature_importances_
        std_arr = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
        indices_arr = np.argsort(importances_arr)[::-1]

        importances_text_list = X.columns.tolist()

        #order_list = [importances_arr[indices_arr[f]] for f in range(X.shape[1])]

        nested_dict = {}

        for f in range(X.shape[1]):
            if f < 10 :
                logging.info("%d. feature %d (%f) %s" % (f + 1, indices_arr[f], importances_arr[indices_arr[f]], importances_text_list[indices_arr[f]]))
            single_feature_dict = {"original_order" : indices_arr[f], "mean_decrease_impurity{}".format(model_type) : importances_arr[indices_arr[f]], "feature{}".format(model_type) : importances_text_list[indices_arr[f]],  "std{}".format(model_type) : std_arr[f]}
            nested_dict[f + 1] = single_feature_dict

        sys.stdout.write("\n\n"), sys.stdout.flush()

        df_imp = pd.DataFrame(nested_dict).T
        df_imp["order_importance{}".format(model_type)] = df_imp.index
        #df_imp.set_index("feature", inplace=True)
        df_imp.set_index("original_order", inplace=True)
        output_dfs.append(df_imp)

    df_imp2 = pd.concat(output_dfs, axis=1)
    df_imp2.sort_values("order_importance", inplace=True)
    df_imp2["original_order"] = df_imp2.index
    df_imp2.set_index("feature", inplace=True)

    df_imp2.to_csv(mean_decrease_impurity_all_features_csv)


def fig_feat_import_from_mean_decrease_impurity(s, logging):
    """Create figures showing ML feature importance.

    Fig1 : Barchart all features
    Fig2 : Barchart top features (currently set at 30)

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    plt.style.use('seaborn-whitegrid')
    plt.rcParams['errorbar.capsize'] = 1
    plt.rcParams.update({'font.size':4})
    colour_dict = create_colour_lists()

    mean_decrease_impurity_all_features_csv = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / "feat_imp/mean_decrease_impurity_all_features.csv"
    variable_importance_all_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "all_var_import.png".format(s["setname"]))
    variable_importance_top_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "top_var_import.png".format(s["setname"]))
    make_sure_path_exists(mean_decrease_impurity_all_features_csv, isfile=True)

    df_imp = pd.read_csv(mean_decrease_impurity_all_features_csv, index_col = 0)

    create_var_imp_plot(df_imp, colour_dict, variable_importance_all_png, df_imp.shape[0])
    #create_var_imp_plot(df_imp, colour_dict, variable_importance_top_png, 30)