import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from random import shuffle
import pandas as pd
import numpy as np
import thoipapy
import sys
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
import os
import joblib
import pickle
import time
from thoipapy.utils import convert_truelike_to_bool, convert_falselike_to_bool, Log_Only_To_Console
from multiprocessing import Pool

# intersect function
def intersect(a, b):
     return list(set(a) & set(b))

def drop_cols_not_used_in_ML(logging, df_data, excel_file_with_settings, i=0):
    """Remove columns not used in machine learning training or testing.

    This includes
        1) text and name columns (e.g. residue_name)
        2) bind columns (interface, interface_score).
        3) columns of non-normalised features, or features that we don't want to include anymore.

    Parameters
    ----------
    logging : logging.Logger
        Python object with settings for logging to console and file.
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns = 'acc_db', "residue_num", "residue_name", etc
        rows = range of number of residues
    excel_file_with_settings : str
        Path to excel file with all settings.
        Necessary for getting the list of features included in the THOIPA algorithm.


    Returns
    -------
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns :
    """
    # read the features tab of the excel settings file
    features_df = pd.read_excel(excel_file_with_settings, sheet_name="features", index_col=0)
    features_df.drop("Notes", axis=0, inplace=True)
    # convert "WAHR" etc to true and false
    features_df["include"] = features_df["include"].apply(convert_truelike_to_bool, convert_nontrue=False)
    features_df["include"] = features_df["include"].apply(convert_falselike_to_bool)

    # print any features that are not shared between the two lists
    unused_features_in_excel = set(features_df.index) - set(df_data.columns)
    features_missing_from_excel = set(df_data.columns) - set(features_df.index)

    # only print this stuff for the first TMD analysed, if in a list
    if i == 0:
        if len(features_missing_from_excel) > 1:
            logging.info("\nfeatures_missing_from_excel, {}".format(features_missing_from_excel))
            if "Unnamed: 0" in features_missing_from_excel:
                raise IndexError("Unnamed column is in dataframe. Try opening csv with index_col=0.")
        if len(unused_features_in_excel) > 1:
            logging.info("\nunused_features_in_excel, {}".format(unused_features_in_excel))
            if len(unused_features_in_excel) > 2:
                logging.warning("There are more than two unused features in the excel settings file. "
                                "This is acceptable for standalone predictor, but otherwise the settings file "
                                "might need to be checked.")
    # drop any features that are not labeled TRUE for inclusion
    features_df = features_df.loc[features_df.include == True]
    # filter df_data to only keep the desired feature columns
    feature_list = features_df.index.tolist()
    df_data = df_data.loc[:, feature_list]
    return df_data

def THOIPA_classifier_with_settings(s, n_features, totally_randomized_trees=False):
    """ For tuning the RF parameters, they are always in one place, and determined by the settings file.

    Parameters
    ----------
    s

    Returns
    -------

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

    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))
    train_data_used_for_model_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data_used_for_model.csv".format(s["setname"]))
    model_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_ML_model.lpkl".format(s["setname"]))

    df_data = pd.read_csv(train_data_csv, index_col=0)

    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]
    df_data = df_data.dropna()

    y = df_data["interface"]
    X = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])

    X.to_csv(train_data_used_for_model_csv)

    if 1 not in y.tolist():
        raise ValueError("None of the residues are marked 1 for an interface residue!")

    n_features = X.shape[1]
    forest = THOIPA_classifier_with_settings(s, n_features)
    # save machine learning model into local driver
    # pkl_file = r'D:\thoipapy\RandomForest\ML_model.lpkl'
    fit = forest.fit(X, y)
    joblib.dump(fit, model_pkl)

    tree_depths = np.array([estimator.tree_.max_depth for estimator in forest.estimators_])
    logging.info("tree depth mean = {} ({})".format(tree_depths.mean(), tree_depths))

    logging.info('finished training machine learning algorithm ({})'.format(model_pkl))


def run_10fold_cross_validation(s, logging):
    """Run 10-fold cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.

    The number of folds is determined by "cross_validation_number_of_splits" in the settings file.

    The number of estimators in the machine learning algorithm is determined by "RF_number_of_estimators" in the settings file,
    therefore it should match the full training result of the whole dataset.

    IMPORTANT. CURRENTLY THERE IS NO AUTOMATIC REDUNDANCY CHECK.
     - homologues of the tested protein could be in the training dataset

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    crossvalidation_pkl : pickle
        Pickled dictionary (xv_dict) containing the results for each fold of validation.
        Also contains the mean ROC curve, and the mean AUC.
    """
    logging.info('10-fold cross validation is running')
    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))
    #crossvalidation_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_10F_data.csv".format(s["setname"]))
    crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_10F_data.pkl".format(s["setname"]))
    features_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_test_features.csv".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(crossvalidation_pkl, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    X = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])
    y = df_data["interface"]

    skf = StratifiedKFold(n_splits=s["cross_validation_number_of_splits"])
    cv = list(skf.split(X, y))

    n_features = X.shape[1]
    forest = THOIPA_classifier_with_settings(s, n_features)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    # save all outputs to a cross-validation dictionary, to be saved as a pickle file
    xv_dict = {}

    start = time.clock()

    for i, (train, test) in enumerate(cv):
        sys.stdout.write("f{}.".format(i+1)), sys.stdout.flush()
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1], drop_intermediate=False)
        xv_dict["fpr{}".format(i)] = fpr
        xv_dict["tpr{}".format(i)] = tpr
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    sys.stdout.write("\n"), sys.stdout.flush()

    logging.info("tree depths : {}".format([estimator.tree_.max_depth for estimator in forest.estimators_]))

    duration = time.clock() - start

    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0

    ROC_AUC = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["ROC_AUC"] = ROC_AUC

    # save dict as pickle
    with open(crossvalidation_pkl, "wb") as f:
        pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    features_ser = pd.Series(X.columns)
    features_ser.to_csv(features_csv)
    logging.info('{} 10-fold validation. AUC({:.3f}). Time taken = {:.2f}.\nFeatures: {}'.format(s["setname"], ROC_AUC, duration, X.columns.tolist()))

def create_10fold_cross_validation_fig(s, logging):
    """Create figure showing ROC curve for each fold in a 10-fold validation.

    The underlying data is created by run_10fold_cross_validation. If this has not been run,
    it will return a file-not-found error.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    #plt.rcParams.update({'font.size': 7})
    crossvalidation_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "{}_10F_ROC.png".format(s["setname"]))
    crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_10F_data.pkl".format(s["setname"]))

    # open pickle file
    with open(crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for i in range(s["cross_validation_number_of_splits"]):
        roc_auc = auc(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)])
        ax.plot(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)], lw=1, label='fold %d (area = %0.2f)' % (i, roc_auc), alpha=0.8)

    ROC_AUC = xv_dict["ROC_AUC"]

    ax.plot(xv_dict["false_positive_rate_mean"], xv_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % ROC_AUC, lw=1.5)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(crossvalidation_png, dpi=240)
    #fig.savefig(crossvalidation_png[:-4] + ".pdf")
    fig.savefig(thoipapy.utils.pdf_subpath(crossvalidation_png))


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
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))
    LOO_crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))
    BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation")
    BO_data_excel = os.path.join(BO_curve_folder, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    acc_db_list = df_data.acc_db.unique()
    start = time.clock()
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])

    n_features = thoipapy.validation.validation.drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"]).shape[1]
    forest = thoipapy.validation.validation.THOIPA_classifier_with_settings(s, n_features)

    #list of all dictionaries for each protein, for multiprocessing
    d_list = []

    if s["use_multiprocessing"]:
        # TURN LOGGING OFF BEFORE MULTIPROCESSING
        logger = Log_Only_To_Console()
    else:
        logger = logging

    val_list = []
    for n, i in enumerate(df_set.index):
        acc, acc_db, database  = df_set.loc[i, "acc"], df_set.loc[i, "acc_db"], df_set.loc[i, "database"]
        THOIPA_prediction_csv = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.{}.LOO.prediction.csv".format(acc, database, s["setname"]))
        testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
        thoipapy.utils.make_sure_path_exists(THOIPA_prediction_csv, isfile=True)
        #######################################################################################################
        #                                                                                                     #
        #      Train data is based on large training csv, after dropping the protein of interest              #
        #                   (positions without closedist/disruption data will be excluded)                    #
        #                                                                                                     #
        #######################################################################################################

        if not acc_db in acc_db_list:
            logging.warning("{} is in protein set, but not found in training data".format(acc_db))
            # skip protein
            continue
        # df_train = df_data.loc[df_data.acc_db != acc_db]
        # X_train = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, df_train, s["excel_file_with_settings"])
        # y_train = df_train[s["bind_column"]]

        #######################################################################################################
        #                                                                                                     #
        #         Pack all variables into a dictionary compatible with multiprocessing Pool                   #
        #                                                                                                     #
        #######################################################################################################

        d = {}
        d["testdata_combined_file"], d["THOIPA_prediction_csv"] = testdata_combined_file, THOIPA_prediction_csv
        #d["X_train"], d["y_train"] = X_train, y_train
        d["df_data"] = df_data
        d["excel_file_with_settings"] = s["excel_file_with_settings"]
        d["forest"], d["pred_colname"] = forest, pred_colname
        d["i"], d["acc_db"], d["database"], d["bind_column"] = i, acc_db, database, s["bind_column"]
        d["n"], d["logger"], d["excel_file_with_settings"] = n, logger, s["excel_file_with_settings"]

        if s["use_multiprocessing"]:
            d_list.append(d)
        else:
            auc_dict, BO_df = LOO_single_prot(d)
            val_tuple = (auc_dict, BO_df)
            val_list.append(val_tuple)

    if s["use_multiprocessing"]:
        with Pool(processes=s["multiple_tmp_simultaneous"]) as pool:
            val_list = pool.map(LOO_single_prot, d_list)

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

def LOO_single_prot(d):
    """Create Leave-One-Out cross-validation for a single protein in a dataset

    see docstring of run_LOO_validation

    Parameters
    ----------
    d : dict
        dictionary with all necessary values for the function
        having a single variable input allows python multiprocessing with Pool
    """
    testdata_combined_file, THOIPA_prediction_csv = d["testdata_combined_file"], d["THOIPA_prediction_csv"]
    #X_train, y_train = d["X_train"], d["y_train"]
    df_data, logger = d["df_data"], d["logger"]
    excel_file_with_settings = d["excel_file_with_settings"]
    forest, pred_colname = d["forest"], d["pred_colname"]
    i, acc_db, database, bind_column = d["i"], d["acc_db"], d["database"], d["bind_column"]

    #######################################################################################################
    #                                                                                                     #
    #                 Train data is based on the large file of all residues in the dataset                #
    #                       - positions without closedist/disruption data will be EXCLUDED)               #
    #                       - crystal positions with "folding/hetero contacts" will be EXCLUDED           #
    #                                                                                                     #
    #######################################################################################################

    df_train = df_data.loc[df_data.acc_db != acc_db]
    X_train = thoipapy.validation.validation.drop_cols_not_used_in_ML(logger, df_train, excel_file_with_settings, i)
    y_train = df_train[bind_column]

    n, logger, excel_file_with_settings = d["n"], d["logger"], d["excel_file_with_settings"]

    #######################################################################################################
    #                                                                                                     #
    #                  Test data is based on the combined features file for that protein and TMD          #
    #                              - positions without closedist/disruption data will be INCLUDED)        #
    #                              - positions with "folding/hetero contacts" will be included            #
    #                                                                                                     #
    #######################################################################################################
    # df_test = df_data.loc[df_data.acc_db == acc_db]
    df_test = pd.read_csv(testdata_combined_file)

    X_test = thoipapy.validation.validation.drop_cols_not_used_in_ML(logger, df_test, excel_file_with_settings, i)
    y_test = df_test["interface"].fillna(0).astype(int)

    #######################################################################################################
    #                                                                                                     #
    #                  Run prediction and save output individually for each protein                       #
    #                                                                                                     #
    #######################################################################################################

    fitted = forest.fit(X_train, y_train)
    if bind_column == "interface":
        prediction = fitted.predict_proba(X_test)[:, 1]
    elif bind_column == "interface_score_norm":
        prediction = fitted.predict(X_test)#[:, 1]
    else:
        raise ValueError("bind_column in excel settings file is not recognised ({})".format(bind_column))
    # add the prediction to the combined file
    df_test[pred_colname] = prediction
    # save just the prediction alone to csv
    prediction_df = df_test[["residue_num", "residue_name", pred_colname]]
    prediction_df.to_csv(THOIPA_prediction_csv, index=False)

    fpr, tpr, thresholds = roc_curve(y_test, prediction, drop_intermediate=False)
    roc_auc = auc(fpr, tpr)

    precision, recall, thresholds_PRC = precision_recall_curve(y_test, prediction)
    pr_auc = auc(recall, precision)

    auc_dict = {"fpr": fpr, "tpr": tpr, "roc_auc": roc_auc, "precision" : precision, "recall" : recall, "pr_auc" : pr_auc}

    if database == "crystal" or database == "NMR":
        # low closest distance means high importance at interface
        df_test["interface_score"] = -1 * df_test["interface_score"]

    BO_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, df_test, experiment_col="interface_score", pred_col=pred_colname)

    if n == 0:
        tree_depths = np.array([estimator.tree_.max_depth for estimator in forest.estimators_])
        logger.info("tree depth mean = {} ({})".format(tree_depths.mean(), tree_depths))
    logger.info("{} AUC : {:.2f}".format(acc_db, roc_auc))

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
    coev_features = feat_list[0:16]
    cons_features = feat_list[16:18]
    polarity_features = feat_list[18:24]
    motif_features = feat_list[24:26]
    pssm_features = feat_list[26:51]
    physical_features = feat_list[51:53]
    TMD_features = feat_list[53:]

    # DEPRECATED in favour of combined polarity_and_pssm_features
    #features_nested_list = [coev_features, cons_features, polarity_features, motif_features, pssm_features, physical_features, TMD_features]
    #features_nested_namelist = ["coev_features", "cons_features", "polarity_features", "motif_features", "pssm_features", "physical_features", "TMD_features"]

    polarity_and_pssm_features = polarity_features + pssm_features
    features_nested_list = [polarity_and_pssm_features, coev_features, cons_features, motif_features, physical_features, TMD_features]
    features_nested_namelist = ["polarity_and_pssm_features", "coev_features", "cons_features", "motif_features", "physical_features", "TMD_features"]

    for i in range(len(features_nested_list)):
        sys.stdout.write("\n{} : {}".format(features_nested_namelist[i], features_nested_list[i]))

    forest = THOIPA_classifier_with_settings(s, n_features)

    pr_auc_orig, roc_auc_orig = calc_mean_AUC(X, y, forest)

    start = time.clock()

    sys.stdout.write("\nmean : {:.03f}\n".format(pr_auc_orig)), sys.stdout.flush()

    decrease_PR_AUC_dict = {}
    decrease_ROC_AUC_dict = {}

    for feature_type, feature_list in zip(features_nested_namelist, features_nested_list):
        logging.info("{} : {}".format(feature_type, feature_list))
        X_t = X.copy()
        for feature in feature_list:
            # shuffle the data for that feature
            row_to_shuffle = X_t[feature].as_matrix()
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
        row_to_shuffle = X_t[feature].as_matrix()
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

def calc_mean_AUC(X, y, forest):
    """Calculate mean precision-recall and ROC AUC using 10-fold cross-validation.

    Parameters
    ----------
    X : pd.DataFrame
    y : pd.Series
    forest : sklearn.ensemble.ExtraTreesClassifier

    Returns
    -------
    (PR_AUC, ROC_AUC) : tuple
        PR_AUC - float of precision-recall AUC
        ROC_AUC - float of ROC-AUC
    """

    skf = StratifiedKFold(n_splits=10)
    cv = list(skf.split(X, y))
    # save precision-recall AUC of each fold
    pr_auc_list = []
    roc_auc_list = []
    for i, (train, test) in enumerate(cv):
        sys.stdout.write("f{}.".format(i + 1)), sys.stdout.flush()
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        
        precision, recall, thresholds_PRC = precision_recall_curve(y.iloc[test], probas_[:, 1])
        pred_auc = auc(recall, precision)
        pr_auc_list.append(pred_auc)

        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1], drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        roc_auc_list.append(roc_auc)
        
    PR_AUC = np.mean(pr_auc_list)
    ROC_AUC = np.mean(roc_auc_list)
    return PR_AUC, ROC_AUC

def fig_feat_import_from_mean_decrease_accuracy(s, logging):

    feat_imp_MDA_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "data", "feat_imp_MDA.csv".format(s["setname"]))
    feat_imp_MDA_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "feat_imp_MDA.png".format(s["setname"]))

    df_fi = pd.read_csv(feat_imp_MDA_csv, index_col=0)


    pass


def create_ROC_all_residues(s, df_set, logging):
    """Combine all residue predictions, so AUC can be calculated from a single array.

    Effectively stacks the CSVs on top of each other.

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
    predictions_csv : csv
        csv file with stacked predictions data for multiple proteins
        index = range(0, ..)
        columns =
    """
    logging.info('Starting combine_all_residue_predictions.')

    # output file with all predictions
    pred_all_res_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "ROC", "{}_pred_all_res.csv".format(s["setname"]))
    #all_res_ROC_data_dict_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "ROC", "{}_all_res_ROC_data_dict.pickle".format(s["setname"]))
    all_res_ROC_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "ROC", "{}_all_res_ROC_data.csv".format(s["setname"]))
    all_res_ROC_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "ROC", "{}_all_res_ROC.png".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(pred_all_res_csv, isfile=True)

    df_set_nonred = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    # set up a dataframe to hold the features for all proteins
    df_all = pd.DataFrame()
    for i in df_set_nonred.index:
        acc = df_set_nonred.loc[i, "acc"]
        database = df_set_nonred.loc[i, "database"]
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))

        df_merged_new_protein = pd.read_csv(merged_data_csv_path, index_col=0)
        df_merged_new_protein["acc_db"] = "{}-{}".format(acc, database)
#
        # for the first protein, replace the empty dataframe
        if df_all.empty:
            df_all = df_merged_new_protein
        else:
            # concatenate the growing dataframe of combined proteins and new dataframe
            df_all = pd.concat([df_all, df_merged_new_protein])

    # drop any positions where there is no interface_score (e.g. no mutations, or hetero contacts?)
    if "interface_score" in df_all.columns:
        df_all.dropna(subset=["interface_score"], inplace=True)
    else:
        logging.warning("No experimental data has been added to this dataset. Hope you're not trying to train with it!!!")

    # reset the index to be a range (0,...).
    df_all.index = range(df_all.shape[0])

    # reorder the columns
    column_list = ['acc_db', 'interface', 'interface_score', 'residue_num', 'residue_name']
    df_all = thoipapy.utils.reorder_dataframe_columns(df_all, column_list)

    df_all.to_csv(pred_all_res_csv)

    save_fig_ROC_all_residues(s, df_all, all_res_ROC_png, all_res_ROC_data_csv, logging)

    df_all["subset"] = df_all.acc_db.str.split("-").str[1]

    subsets = ["ETRA", "NMR", "crystal"]
    for subset in subsets:
        df_subset = df_all.loc[df_all.subset == subset]
        if not df_subset.empty:
            ROC_png = all_res_ROC_png[:-4] + "_{}_subset.png".format(subset)
            ROC_data_csv = all_res_ROC_data_csv[:-4] + "_{}_subset.csv".format(subset)
            save_fig_ROC_all_residues(s, df_subset, ROC_png, ROC_data_csv, logging)

    # with open(all_res_ROC_data_pkl, "wb") as f:
    #     pickle.dump(output_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Finished combine_all_residue_predictions.')


def save_fig_ROC_all_residues(s, df, all_res_ROC_png, all_res_ROC_data_csv, logging):
    fontsize=8
    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot([0, 1], [0, 1], color="0.5", linestyle="--", label="random", linewidth=1)
    THOIPA_predictor = "THOIPA_{}_LOO".format(s["set_number"])
    predictors = [THOIPA_predictor, "TMDOCK", "LIPS_surface_ranked", "PREDDIMER"]
    output_dict = {}
    for predictor in predictors:
        df_sel = df[["interface", predictor]].dropna()
        if predictor in ["TMDOCK", "PREDDIMER"]:
            pred = - df_sel[predictor]
            # pred = normalise_between_2_values(df_sel[predictor], 2.5, 8, invert=True)
        else:
            pred = df_sel[predictor]
        fpr, tpr, thresholds = roc_curve(df_sel.interface, pred, drop_intermediate=False)
        pred_auc = auc(fpr, tpr)
        #sys.stdout.write("{} AUC : {:.03f}\n".format(predictor, pred_auc))
        label = "{}. AUC : {:.03f}".format(predictor, pred_auc)
        ax.plot(fpr, tpr, label=label, linewidth=1)
        # output_dict["fpr_{}".format(predictor)] = fpr
        # output_dict["tpr_{}".format(predictor)] = tpr
        # output_dict["auc_{}".format(predictor)] = auc

        output_dict[predictor] = {"fpr" : list(fpr), "tpr" : list(tpr), "pred_auc" : pred_auc}
    ax.grid(False)
    ax.set_xlabel("false positive rate", fontsize=fontsize)
    ax.set_ylabel("true positive rate", fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(all_res_ROC_png, dpi=240)
    fig.savefig(all_res_ROC_png[:-4] + ".pdf")

    df_ROC_data = pd.DataFrame(output_dict).T
    df_ROC_data.to_csv(all_res_ROC_data_csv)

    logging.info("save_fig_ROC_all_residues finished ({})".format(all_res_ROC_data_csv))


def create_precision_recall_all_residues(s, df_set, logging):
    """Combine all residue predictions, so precision recall can be calculated from a single array.

    Effectively stacks the CSVs on top of each other.

    Code is directly copied and modified from create_ROC_all_residues

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
    predictions_csv : csv
        csv file with stacked predictions data for multiple proteins
        index = range(0, ..)
        columns =
    """
    logging.info('Starting combine_all_residue_predictions.')

    # output file with all predictions
    pred_all_res_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "precision_recall", "{}_pred_all_res.csv".format(s["setname"]))
    #all_res_precision_recall_data_dict_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "precision_recall", "{}_all_res_precision_recall_data_dict.pickle".format(s["setname"]))
    all_res_precision_recall_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "precision_recall", "{}_all_res_precision_recall_data.csv".format(s["setname"]))
    all_res_precision_recall_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "precision_recall", "{}_all_res_precision_recall.png".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(pred_all_res_csv, isfile=True)

    df_set_nonred = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    # set up a dataframe to hold the features for all proteins
    df_all = pd.DataFrame()
    for i in df_set_nonred.index:
        acc = df_set_nonred.loc[i, "acc"]
        database = df_set_nonred.loc[i, "database"]
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))

        df_merged_new_protein = pd.read_csv(merged_data_csv_path, index_col=0)
        df_merged_new_protein["acc_db"] = "{}-{}".format(acc, database)
#
        # for the first protein, replace the empty dataframe
        if df_all.empty:
            df_all = df_merged_new_protein
        else:
            # concatenate the growing dataframe of combined proteins and new dataframe
            df_all = pd.concat([df_all, df_merged_new_protein])

    # drop any positions where there is no interface_score (e.g. no mutations, or hetero contacts?)
    if "interface_score" in df_all.columns:
        df_all.dropna(subset=["interface_score"], inplace=True)
    else:
        logging.warning("No experimental data has been added to this dataset. Hope you're not trying to train with it!!!")

    # reset the index to be a range (0,...).
    df_all.index = range(df_all.shape[0])

    # reorder the columns
    column_list = ['acc_db', 'interface', 'interface_score', 'residue_num', 'residue_name']
    df_all = thoipapy.utils.reorder_dataframe_columns(df_all, column_list)

    df_all.to_csv(pred_all_res_csv)

    save_fig_precision_recall_all_residues(s, df_all, all_res_precision_recall_png, all_res_precision_recall_data_csv, logging)

    df_all["subset"] = df_all.acc_db.str.split("-").str[1]

    subsets = ["ETRA", "NMR", "crystal"]
    for subset in subsets:
        df_subset = df_all.loc[df_all.subset == subset]
        if not df_subset.empty:
            precision_recall_png = all_res_precision_recall_png[:-4] + "_{}_subset.png".format(subset)
            precision_recall_data_csv = all_res_precision_recall_data_csv[:-4] + "_{}_subset.csv".format(subset)
            save_fig_precision_recall_all_residues(s, df_subset, precision_recall_png, precision_recall_data_csv, logging)

    # with open(all_res_precision_recall_data_pkl, "wb") as f:
    #     pickle.dump(output_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Finished combine_all_residue_predictions.')


def save_fig_precision_recall_all_residues(s, df, all_res_precision_recall_png, all_res_precision_recall_data_csv, logging):
    """Save figure for precision recall plot of all residues joined together.

    Code is directly copied and modified from save_fig_ROC_all_residues

    """
    fontsize=8
    fig, ax = plt.subplots(figsize=(5,5))
    THOIPA_predictor = "THOIPA_{}_LOO".format(s["set_number"])
    predictors = [THOIPA_predictor, "TMDOCK", "LIPS_surface_ranked", "PREDDIMER", "random"]
    output_dict = {}
    interface_random = df.interface_score.tolist()
    shuffle(interface_random)
    df["random"] = interface_random

    for predictor in predictors:
        df_sel = df[["interface", predictor]].dropna()
        if predictor in ["TMDOCK", "PREDDIMER"]:
            pred = - df_sel[predictor]
            # pred = normalise_between_2_values(df_sel[predictor], 2.5, 8, invert=True)
        else:
            pred = df_sel[predictor]
        precision, recall, thresholds_PRC = precision_recall_curve(df_sel.interface, pred)

        pred_auc = auc(recall, precision)
        #sys.stdout.write("{} AUC : {:.03f}\n".format(predictor, pred_auc))
        label = "{}. AUC : {:.03f}".format(predictor, pred_auc)
        ax.plot(recall, precision, label=label, linewidth=1)

        output_dict[predictor] = {"precision" : list(precision), "recall" : list(recall), "pred_auc" : pred_auc}
    ax.grid(False)

    ax.set_xlabel("recall", fontsize=fontsize)
    ax.set_ylabel("precision", fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(all_res_precision_recall_png, dpi=240)
    fig.savefig(all_res_precision_recall_png[:-4] + ".pdf")

    df_precision_recall_data = pd.DataFrame(output_dict).T
    df_precision_recall_data.to_csv(all_res_precision_recall_data_csv)

    logging.info("save_fig_precision_recall_all_residues finished ({})".format(all_res_precision_recall_data_csv))