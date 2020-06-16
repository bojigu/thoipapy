import os
import pickle
import sys
import time
from ast import literal_eval
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp
from sklearn.metrics import auc, roc_curve, precision_recall_curve

import thoipapy.figs
import thoipapy.utils
import thoipapy.validation
from thoipapy.utils import Log_Only_To_Console


class LooValidationData:

    def __init__(self):
        self.testdata_combined_file = None
        self.THOIPA_prediction_csv = None
        self.df_data = None
        self.excel_file_with_settings = None
        self.forest = None
        self.pred_colname = None
        self.i = None
        self.acc_db = None
        self.database = None
        self.bind_column = None
        self.n = None
        self.logger = None
        self.excel_file_with_settings = None


def run_LOO_validation(s, df_set, logging):
    """Run Leave-One-Out cross-validation for a particular set of TMDs (e.g. set05).

    The SAME SET is used for both training and cross-validation.
    ONLY USE IF YOUR DATASET IS NON-REDUNDANT! OR USE AUTO CD-HIT REDUNDANCY CHECKS!
    Each protein is uniquely identified by the acc and the database (acc_db).

    The training set (df_train) consists of all proteins except the one tested.
    The test dataset (df_test) contains only the interface and features of the one test protein

    The model is trained on the train data csv (e.g. set38_train_data.csv)
         - for crystal subset, hetero contacts (folding residues) are removed from this training set!

    The model is validated against each combined CSV with features (e.g. "D:\data_thoipapy\Features\combined\ETRA\Q12983.surr20.gaps5.combined_features.csv")
     - for crystal subset, folding residues (hetero contacts) are INCLUDED here.
     - the model created without the folding contacts is therefore validated against the full seq, including folding residues

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    LOO_crossvalidation_pkl, pickle
        Pickled dictionary (xv_dict) containing the results for each fold of validation.
        Also contains the mean ROC curve, and the mean AUC.
    BO_all_data_csv, csv
        CSV with the BO curve underlying data
    BO_data_excel, csv
        excel file with the processed BO-curve data
    """
    logging.info('Leave-One-Out cross validation is running')
    setname = s["setname"]
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))
    LOO_crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))
    BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation")
    BO_data_excel = os.path.join(BO_curve_folder, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))
    sim_matrix_xlsx = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / f"crossvalidation/clusters/{setname}_sim_matrix.xlsx"

    if not sim_matrix_xlsx.is_file():
        raise FileNotFoundError(f"The similarity matrix with clusters of putative homologues could not be found ({sim_matrix_xlsx})")

    thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    assert "Unnamed" not in ", ".join(df_data.columns.tolist())

    #df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    acc_db_list = df_data.acc_db.unique()
    logging.info(f"Dataset has {len(acc_db_list)} unique proteins for training.")
    start = time.clock()
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])

    n_features = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"]).shape[1]
    forest = thoipapy.validation.train_model.THOIPA_classifier_with_settings(s, n_features)

    #list of all dictionaries for each protein, for multiprocessing
    loo_validation_data_list = []

    if s["use_multiprocessing"]:
        # TURN LOGGING OFF BEFORE MULTIPROCESSING
        logger = Log_Only_To_Console()
    else:
        logger = logging

    df_clusters = pd.read_excel(sim_matrix_xlsx, sheet_name="reduced_clusters", index_col=0)
    df_clusters["reduced_clusters"] = df_clusters["reduced_clusters"].apply(lambda x: literal_eval(x))
    df_clusters["acc_database_list"] = df_clusters["reduced_clusters"].apply(lambda x: [y[2:] for y in list(x)])
    cluster = df_clusters.at[1, "acc_database_list"]

    val_list = []
    for n, i in enumerate(df_set.index):
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        loo_validation_data = LooValidationData()
        loo_validation_data.acc = acc
        loo_validation_data.acc_db = df_set.loc[i, "acc_db"]
        loo_validation_data.bind_column = s["bind_column"]
        loo_validation_data.database = database
        loo_validation_data.df_data = df_data
        loo_validation_data.excel_file_with_settings = s["excel_file_with_settings"]
        loo_validation_data.forest = forest
        loo_validation_data.i = i
        loo_validation_data.logger = logger
        loo_validation_data.n = n
        loo_validation_data.pred_colname = pred_colname
        loo_validation_data.testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
        loo_validation_data.THOIPA_prediction_csv = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / f"crossvalidation/leave_one_out/protein_data/{acc}.{database}.LOO.prediction.csv"

        thoipapy.utils.make_sure_path_exists(loo_validation_data.THOIPA_prediction_csv, isfile=True)
        #######################################################################################################
        #                                                                                                     #
        #      Train data is based on large training csv, after dropping the protein of interest              #
        #                   (positions without closedist/disruption data will be excluded)                    #
        #                                                                                                     #
        #######################################################################################################

        if not loo_validation_data.acc_db in acc_db_list:
            logging.warning("{} is in protein set, but not found in training data".format(loo_validation_data.acc_db))
            # skip protein
            continue
        # df_train = df_data.loc[df_data.acc_db != acc_db]
        # X_train = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, df_train, s["excel_file_with_settings"])
        # y_train = df_train[s["bind_column"]]

        #######################################################################################################
        #                                                                                                     #
        #           Pack all variables into an object compatible with multiprocessing Pool                    #
        #                                                                                                     #
        #######################################################################################################

        if s["use_multiprocessing"]:
            loo_validation_data_list.append(loo_validation_data)
        else:
            auc_dict, BO_df = LOO_single_prot(loo_validation_data)
            val_tuple = (auc_dict, BO_df)
            val_list.append(val_tuple)

    if s["use_multiprocessing"]:
        with Pool(processes=s["multiple_tmp_simultaneous"]) as pool:
            val_list = pool.map(LOO_single_prot, loo_validation_data_list)

    #######################################################################################################
    #                                                                                                     #
    #                            Get mean AUC etc for all proteins in list                                #
    #                                                                                                     #
    #######################################################################################################
    duration = time.clock() - start
    sys.stdout.write("\n")

    # copied from original mean_tpr code
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)

    BO_all_df = pd.DataFrame()
    all_roc_auc = []
    all_pr_auc = []
    xv_dict = {}
    acc_db_list = df_set.acc_db.tolist()

    #iterate through the output tuple (auc_dict, BO_df)
    for nn, val_tuple in enumerate(val_list):
        acc_db = acc_db_list[nn]
        auc_dict = val_tuple[0]
        all_roc_auc.append(auc_dict["roc_auc"])
        all_pr_auc.append(auc_dict["pr_auc"])
        BO_df = val_tuple[1]
        # join the data for all BO curves together
        if nn == 0:
            BO_all_df = BO_df
        else:
            BO_all_df = pd.concat([BO_all_df, BO_df], axis=1, join="outer")
        mean_tpr += interp(mean_fpr, auc_dict["fpr"], auc_dict["tpr"])
        mean_tpr[0] = 0.0

        xv_dict[acc_db] = {"fpr": auc_dict["fpr"], "tpr": auc_dict["tpr"], "roc_auc": auc_dict["roc_auc"]}

    # copied from original mean_tpr code
    mean_tpr /= df_set.shape[0]
    mean_tpr[-1] = 1.0
    mean_roc_auc_from_joined_data = auc(mean_fpr, mean_tpr)

    # calculate mean of each protein AUC separately
    mean_roc_auc_all_prot = np.array(all_roc_auc).mean()
    xv_dict["mean_roc_auc_all_prot"] = mean_roc_auc_all_prot
    mean_pr_auc_all_prot = np.array(all_pr_auc).mean()
    xv_dict["mean_pr_auc_all_prot"] = mean_pr_auc_all_prot

    # add to dict that can be used for figure creation later
    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["mean_roc_auc_from_joined_data"] = mean_roc_auc_from_joined_data
    xv_dict["mean_roc_auc_all_prot"] = mean_roc_auc_all_prot

    # save dict as pickle
    thoipapy.utils.make_sure_path_exists(LOO_crossvalidation_pkl, isfile=True)
    with open(LOO_crossvalidation_pkl, "wb") as f:
        pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    #######################################################################################################
    #                                                                                                     #
    #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
    #                                                                                                     #
    #######################################################################################################

    BO_all_df.to_csv(BO_all_data_csv)
    #names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

    #linechart_mean_obs_and_rand = thoipapy.figs.Create_Bo_Curve_files.analyse_bo_curve_underlying_data(BO_all_data_csv, crossvalidation_folder, names_excel_path)
    thoipapy.figs.create_BOcurve_files.parse_BO_data_csv_to_excel(BO_all_data_csv, BO_data_excel, logging)

    logging.info('{} LOO crossvalidation. Time taken = {:.2f}.'.format(s["setname"], duration))
    logging.info('---ROC_AUC(mean each protein : {:.2f})(from joined data {:.2f})---'.format(mean_roc_auc_all_prot, mean_roc_auc_from_joined_data))
    logging.info('---PR_AUC(mean each protein : {:.2f})---'.format(mean_pr_auc_all_prot))


def LOO_single_prot(d: LooValidationData):
    """Create Leave-One-Out cross-validation for a single protein in a dataset

    see docstring of run_LOO_validation

    Parameters
    ----------
    d : dict
        dictionary with all necessary values for the function
        having a single variable input allows python multiprocessing with Pool
    """
    logger = d.logger

    #######################################################################################################
    #                                                                                                     #
    #                 Train data is based on the large file of all residues in the dataset                #
    #                       - positions without closedist/disruption data will be EXCLUDED)               #
    #                       - crystal positions with "folding/hetero contacts" will be EXCLUDED           #
    #                                                                                                     #
    #######################################################################################################

    df_train = d.df_data.loc[d.df_data.acc_db != d.acc_db]
    X_train = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(d.logger, df_train, d.excel_file_with_settings, d.i)
    y_train = df_train[d.bind_column]

    #n, logger, excel_file_with_settings = d["n"], d["logger"], d["excel_file_with_settings"]

    #######################################################################################################
    #                                                                                                     #
    #                  Test data is based on the combined features file for that protein and TMD          #
    #                              - positions without closedist/disruption data will be INCLUDED)        #
    #                              - positions with "folding/hetero contacts" will be included            #
    #                                                                                                     #
    #######################################################################################################
    # df_test = df_data.loc[df_data.acc_db == acc_db]
    df_test = pd.read_csv(d.testdata_combined_file, index_col=0)
    assert "Unnamed" not in ",".join(df_test.columns.tolist())
    assert df_test.index.tolist() == list(range(df_test.shape[0]))

    X_test = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logger, df_test, d.excel_file_with_settings, d.i)
    y_test = df_test["interface"].fillna(0).astype(int)

    #######################################################################################################
    #                                                                                                     #
    #                  Run prediction and save output individually for each protein                       #
    #                                                                                                     #
    #######################################################################################################

    fitted = d.forest.fit(X_train, y_train)
    if d.bind_column == "interface":
        prediction = fitted.predict_proba(X_test)[:, 1]
    elif d.bind_column == "interface_score_norm":
        prediction = fitted.predict(X_test)#[:, 1]
    else:
        raise ValueError("bind_column in excel settings file is not recognised ({})".format(d.bind_column))
    # add the prediction to the combined file
    df_test[d.pred_colname] = prediction
    # save just the prediction alone to csv
    prediction_df = df_test[["residue_num", "residue_name", d.pred_colname]]
    prediction_df.to_csv(d.THOIPA_prediction_csv, index=False)

    fpr, tpr, thresholds = roc_curve(y_test, prediction, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    precision, recall, thresholds_PRC = precision_recall_curve(y_test, prediction)
    pr_auc = auc(recall, precision)

    auc_dict = {"fpr": fpr, "tpr": tpr, "roc_auc": roc_auc, "precision" : precision, "recall" : recall, "pr_auc" : pr_auc}

    if d.database == "crystal" or d.database == "NMR":
        # low closest distance means high importance at interface
        df_test["interface_score"] = -1 * df_test["interface_score"]

    BO_df = thoipapy.figs.fig_utils.calc_best_overlap(d.acc_db, df_test, experiment_col="interface_score", pred_col=d.pred_colname)

    if d.n == 0:
        tree_depths = np.array([estimator.tree_.max_depth for estimator in d.forest.estimators_])
        logger.info("tree depth mean = {} ({})".format(tree_depths.mean(), tree_depths))
    logger.info("{} AUC : {:.2f}".format(d.acc_db, roc_auc))

    return auc_dict, BO_df


def create_LOO_validation_fig(s, df_set, logging):
    """Create Leave-One-Out cross-validation for each TMD in a dataset.

    Training dataset = all residues in full dataset, except that being trained.

    ASSUMES YOUR DATASET IS NON-REDUNDANT! AUTO CD-HIT REDUNDANCY CHECKS ARE MOSTLY UNUSED!

    The model is trained on the train data csv (e.g. set38_train_data.csv)
         - for crystal subset, hetero contacts (folding residues) are removed from this training set!

    The model is validated against each combined CSV with features
     - for crystal subset, folding residues (hetero contacts) are INCLUDED here.
     - the model created without the folding contacts is therefore validated against the full seq, including folding residues

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    #plt.rcParams.update({'font.size': 7})
    LOO_crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    LOO_crossvalidation_ROC_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "{}_LOO_crossvalidation_ROC.png".format(s["setname"]))
    LOO_crossvalidation_AUC_bar_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "{}_LOO_crossvalidation_AUC_bar.png".format(s["setname"]))
    AUC_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_AUC.csv".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))
    BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation")
    BO_data_excel = os.path.join(BO_curve_folder, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))
    BO_linechart_png = os.path.join(BO_curve_folder, "{}_BO_linechart.png".format(s["setname"]))
    BO_barchart_png = os.path.join(BO_curve_folder, "{}_LOO_AUBOC10_barchart.png".format(s["setname"]))

    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    # open pickle file
    with open(LOO_crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    # due to problems on Bo's computer, set the figsize to double what we should be using for the publication?
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    auc_dict = {}

    for acc_db in df_set.acc_db:
        if acc_db in xv_dict:
            roc_auc = xv_dict[acc_db]["roc_auc"]
            auc_dict[acc_db] = roc_auc
            ax.plot(xv_dict[acc_db]["fpr"], xv_dict[acc_db]["tpr"], lw=1, label='{} ({:.2f})'.format(acc_db, roc_auc), alpha=0.8)
        else:
            logging.warning("{} not in xv_dict after LOO validation".format(acc_db))

    mean_roc_auc_all_prot = xv_dict["mean_roc_auc_all_prot"]

    ax.plot(xv_dict["false_positive_rate_mean"], xv_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_roc_auc_all_prot, lw=1.5)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(LOO_crossvalidation_ROC_png, dpi=240)
    #fig.savefig(LOO_crossvalidation_ROC_png[:-4] + ".pdf")
    fig.savefig(thoipapy.utils.pdf_subpath(LOO_crossvalidation_ROC_png))

    AUC_ser = pd.Series(auc_dict)
    AUC_ser.to_csv(AUC_csv)
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    AUC_ser.plot(ax=ax, kind="bar")
    ax.set_ylabel("performance (AUBOC10)")
    fig.tight_layout()
    fig.savefig(LOO_crossvalidation_AUC_bar_png, dpi=240)
    fig.savefig(thoipapy.utils.pdf_subpath(LOO_crossvalidation_AUC_bar_png))

    AUBOC10 = thoipapy.figs.create_BOcurve_files.save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser)

    if "TRUE" == "TRUE":
        other_figs_path = os.path.join(BO_curve_folder, "other_figs")
        thoipapy.utils.make_sure_path_exists(other_figs_path)
        thoipapy.figs.create_BOcurve_files.save_extra_BO_figs(BO_data_excel, other_figs_path)

    logging.info('{} LOO crossvalidation. AUBOC10({:.2f}).'.format(s["setname"], AUBOC10))
    logging.info("create_LOO_validation_fig finished ({})".format(LOO_crossvalidation_AUC_bar_png))