import os
import sys
import time

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import thoipapy.utils
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML
from thoipapy.validation.train_model import THOIPA_classifier_with_settings
from thoipapy.validation.auc import calc_mean_AUC


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
    variable_importance_csv : csv
        List of variables, sorted by their importance to the algorithm.
        Also includes the standard deviation supplied by the machine learning algorithm
    """
    #logging.info('RF_variable_importance_calculate is running\n')
    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))

    variable_importance_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "variable_importance.csv".format(s["setname"]))
    thoipapy.utils.make_sure_path_exists(variable_importance_csv, isfile=True)

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

    df_imp2.to_csv(variable_importance_csv)


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
    #from thoipapy.utils import create_colour_lists
    from thoipapy.utils import create_colour_lists
    colour_dict = create_colour_lists()

    variable_importance_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "variable_importance.csv".format(s["setname"]))
    variable_importance_all_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "all_var_import.png".format(s["setname"]))
    variable_importance_top_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "top_var_import.png".format(s["setname"]))
    thoipapy.utils.make_sure_path_exists(variable_importance_csv, isfile=True)

    df_imp = pd.read_csv(variable_importance_csv, index_col = 0)

    create_var_imp_plot(df_imp, colour_dict, variable_importance_all_png, df_imp.shape[0])
    #create_var_imp_plot(df_imp, colour_dict, variable_importance_top_png, 30)


def create_var_imp_plot(df_imp, colour_dict, variable_importance_png, n_features_in_plot):
    """Plot function for fig_feat_import_from_mean_decrease_impurity, allowing a variable number of features.
    """
    df_sel = df_imp.iloc[:n_features_in_plot, :].copy()

    # regular trees, or Totally Randomized Trees
    model_types = ["", "_TRT"]

    for model_type in model_types:

        # add suffix for the totally randomised trees
        variable_importance_png = variable_importance_png[:-4] + model_type + ".png"

        df_sel.sort_values("mean_decrease_impurity{}".format(model_type), ascending=True, inplace=True)
        #min_ = df_sel.mean_decrease_impurity.min()

        # determine the plot height by the number of features
        # currently set for 30
        plot_height = 4 * n_features_in_plot / 30
        figsize = np.array([4.42, plot_height])
        fig, ax = plt.subplots(figsize=figsize)

        TUMblue = colour_dict["TUM_colours"]['TUMBlue']
        df_sel["mean_decrease_impurity{}".format(model_type)].plot(kind="barh", color="#17a8a5", ax=ax)# xerr=df_sel["std"]
        ax.errorbar(df_sel["mean_decrease_impurity{}".format(model_type)], range(len(df_sel.index)), xerr=df_sel["std{}".format(model_type)], fmt="none", ecolor="k", ls="none", capthick=0.5, elinewidth=0.5, capsize=1, label=None)

        ax.set_xlim(0)

        ax.set_ylabel("")
        ax.set_xlabel("variable importance\n(mean decrease impurity)")
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(variable_importance_png, dpi=240)
        #fig.savefig(variable_importance_png[:-4] + ".pdf")
        fig.savefig(thoipapy.utils.pdf_subpath(variable_importance_png))


def calc_feat_import_from_mean_decrease_accuracy(s, logging):
    """Calculate feature importances using mean decrease in accuracy.

    This method differs from calc_feat_import_from_mean_decrease_impurity.
    It's much slower, and involves the use of 10-fold cross-validation for each variable separately.

     - a feature (or group of features) is selected for randomisation
     - in theory, randomising important features will cause a drop in prediction accuracy
     - The feature (or group of features) is shuffled
     - precision-recall AUC and ROC-AUC is measured
     - the difference between the original AUC and the AUC with shuffled variable is measured
     - higher values suggest more important features

    feature groups:
    polarity_and_pssm_features : ['polarity', 'relative_polarity', 'polarity4mean', 'polarity3Nmean', 'polarity3Cmean', 'polarity1mean', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'CS', 'DE', 'KR', 'QN', 'LIV']
    coev_features : ['DImax', 'MImax', 'DItop4mean', 'MItop4mean', 'DItop8mean', 'MItop8mean', 'DI4max', 'MI4max', 'DI1mean', 'MI1mean', 'DI3mean', 'MI3mean', 'DI4mean', 'MI4mean', 'DI4cum', 'MI4cum']
    cons_features : ['conservation', 'cons4mean']
    motif_features : ['GxxxG', 'SmxxxSm']
    physical_features : ['branched', 'mass']
    TMD_features : ['residue_depth', 'n_TMDs', 'n_homologues']

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    feat_imp_MDA_csv : csv
        Comma separated values, showing decrease in AUC for each feature or group of features.
    """
    logging.info('calc_feat_import_from_mean_decrease_accuracy is running')
    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))
    # crossvalidation_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_10F_data.csv".format(s["setname"]))
    feat_imp_MDA_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "data", "feat_imp_MDA.csv".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(feat_imp_MDA_csv, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    X = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])
    y = df_data["interface"]

    n_features = X.shape[1]
    feat_list = X.columns.tolist()
    # coev_features = feat_list[0:16]
    # cons_features = feat_list[16:18]
    # polarity_features = feat_list[18:24]
    # motif_features = feat_list[24:26]
    # pssm_features = feat_list[26:51]
    # physical_features = feat_list[51:53]
    # TMD_features = feat_list[53:]

    # DEPRECATED in favour of combined polarity_and_pssm_features
    #features_nested_list = [coev_features, cons_features, polarity_features, motif_features, pssm_features, physical_features, TMD_features]
    #features_nested_namelist = ["coev_features", "cons_features", "polarity_features", "motif_features", "pssm_features", "physical_features", "TMD_features"]

    polarity_features = ["test_dropping_of_features_not_included", "polarity", "relative_polarity", "polarity4mean", "polarity3Nmean", "polarity3Cmean", "polarity1mean"]
    pssm_features = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "CS", "DE", "KR", "QN", "LIV"]
    coev_features = ["DImax", "MImax", "DItop4mean", "MItop4mean", "DItop8mean", "MItop8mean", "DI4max", "MI4max", "DI1mean", "MI1mean", "DI3mean", "MI3mean", "DI4mean", "MI4mean", "DI4cum", "MI4cum"]
    DI_features = ["DImax", "DItop4mean", "DItop8mean", "DI4max", "DI1mean", "DI3mean", "DI4mean", "DI4cum"]
    MI_features = ["MImax", "MItop4mean", "MItop8mean", "MI4max", "MI1mean", "MI3mean", "MI4mean", "MI4cum"]
    cons_features = ["conservation", "cons4mean", "rate4site"]
    motif_features =  ["GxxxG", "SmxxxSm"]
    physical_features = ["branched", "mass"]
    TMD_features = ["residue_depth", "n_TMDs", "n_homologues"]
    polarity_and_pssm_features = polarity_features + pssm_features
    features_nested_list = [polarity_and_pssm_features, coev_features, DI_features, MI_features, cons_features, motif_features, physical_features, TMD_features]
    features_nested_namelist = ["polarity_and_pssm_features", "coev_features",  "DI_features", "MI_features", "cons_features", "motif_features", "physical_features", "TMD_features"]

    for i in range(len(features_nested_list)):
        sys.stdout.write("\n{} : {}".format(features_nested_namelist[i], features_nested_list[i]))

    forest = THOIPA_classifier_with_settings(s, n_features)

    pr_auc_orig, roc_auc_orig = calc_mean_AUC(X, y, forest)

    start = time.clock()

    sys.stdout.write("\nmean : {:.03f}\n".format(pr_auc_orig)), sys.stdout.flush()

    decrease_PR_AUC_dict = {}
    decrease_ROC_AUC_dict = {}

    for feature_type, feature_list in zip(features_nested_namelist, features_nested_list):
        feature_list = list(set(feature_list).intersection(set(X.columns.tolist())))
        logging.info("{} : {}".format(feature_type, feature_list))
        X_t = X.copy()
        for feature in feature_list:
            # shuffle the data for that feature
            row_to_shuffle = X_t[feature].to_numpy()
            np.random.shuffle(row_to_shuffle)
            X_t[feature] = row_to_shuffle
        # calculate prediction performance after shuffling
        PR_AUC, ROC_AUC = calc_mean_AUC(X_t, y, forest)

        decrease_PR_AUC = pr_auc_orig - PR_AUC
        decrease_PR_AUC_dict[feature_type] = decrease_PR_AUC

        decrease_ROC_AUC = roc_auc_orig - ROC_AUC
        decrease_ROC_AUC_dict[feature_type] = decrease_ROC_AUC

        logging.info("  {} {:.03f} {:.03f}".format(feature_type, decrease_PR_AUC, decrease_ROC_AUC))

    for feature in X.columns:
        X_t = X.copy()
        # shuffle the data for that feature
        row_to_shuffle = X_t[feature].to_numpy()
        np.random.shuffle(row_to_shuffle)
        X_t[feature] = row_to_shuffle
        # calculate prediction performance after shuffling
        PR_AUC, ROC_AUC = calc_mean_AUC(X_t, y, forest)

        decrease_PR_AUC = pr_auc_orig - PR_AUC
        decrease_PR_AUC_dict[feature] = decrease_PR_AUC

        decrease_ROC_AUC = roc_auc_orig - ROC_AUC
        decrease_ROC_AUC_dict[feature] = decrease_ROC_AUC

        logging.info("  {} {:.03f} {:.03f}".format(feature, decrease_PR_AUC, decrease_ROC_AUC))

    df_fi = pd.DataFrame()
    df_fi["PR_AUC"] = pd.Series(decrease_PR_AUC_dict)
    df_fi["ROC_AUC"] = pd.Series(decrease_ROC_AUC_dict)

    df_fi.to_csv(feat_imp_MDA_csv)

    duration = time.clock() - start

    logging.info('{} calc_feat_import_from_mean_decrease_accuracy. PR_AUC({:.3f}). Time taken = {:.2f}.\nFeatures: {}'.format(s["setname"], pr_auc_orig, duration, X.columns.tolist()))


def fig_feat_import_from_mean_decrease_accuracy(s, logging):

    feat_imp_MDA_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "data", "feat_imp_MDA.csv".format(s["setname"]))
    feat_imp_MDA_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "feat_imp_MDA.png".format(s["setname"]))

    df_fi = pd.read_csv(feat_imp_MDA_csv, index_col=0)


    pass