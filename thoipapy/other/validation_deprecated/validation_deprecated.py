import os
import pickle
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp
import joblib
from sklearn.metrics import roc_curve, auc

import thoipapy
import thoipapy.validation.feature_selection
from thoipapy.utils import get_testsetname_trainsetname_from_run_settings
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML


def run_LOO_validation_od_non_multiprocessing(s, df_set, logging):
    """Run Leave-One-Out cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.
    Each protein is uniquely identified by the acc and the database (acc_db).
    The test and training datasets are based on the following:
        1) the df_set derived from the set of protein sequences, e.g. set03
            - created manually
        2) the train_data csv, e.g."D:\data_thoipapy\Results\set03\set03_train_data.csv".
            - filtered according to redundancy etc by combine_all_train_data_for_machine_learning()
            - this requires a CD-HIT file for this dataset to remove redundant proteins
    If the acc_db is not in BOTH of these locations, it will not be used for training and validation.

    The training set (df_train) consists of all proteins except the one tested.
    The test dataset (df_test) contains only the interface and features of the one test protein

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
    crossvalidation_pkl : pickle
        Pickled dictionary (xv_dict) containing the results for each fold of validation.
        Also contains the mean ROC curve, and the mean AUC.
    """
    logging.info('Leave-One-Out cross validation is running')
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    train_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/01_train_data_orig.csv"
    crossvalidation_folder = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation")
    LOO_crossvalidation_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))
    BO_data_excel: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/data/{s['setname']}_thoipa_loo_bo_curve_data.xlsx"

    thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    if s["min_n_homol_training"] != 0:
        df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    acc_db_list = df_data.acc_db.unique()
    xv_dict = {}
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    start = time.clock()
    BO_all_df = pd.DataFrame()
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])

    n_features = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"]).shape[1]
    forest = thoipapy.validation.train_model.THOIPA_classifier_with_settings(s, n_features)

    for i in df_set.index:

        acc, acc_db, database  = df_set.loc[i, "acc"], df_set.loc[i, "acc_db"], df_set.loc[i, "database"]
        THOIPA_LOO_prediction_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/predictions/THOIPA_LOO/{database}.{acc}.LOO.prediction.csv"
        thoipapy.utils.make_sure_path_exists(THOIPA_LOO_prediction_csv, isfile=True)

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
        df_train = df_data.loc[df_data.acc_db != acc_db]
        X_train = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, df_train, s["excel_file_with_settings"], i)
        y_train = df_train[s["bind_column"]]

        #######################################################################################################
        #                                                                                                     #
        #                  Test data is based on the combined features file for that protein and TMD          #
        #                   (positions without closedist/disruption data will be INCLUDED)                    #
        #                                                                                                     #
        #######################################################################################################
        #df_test = df_data.loc[df_data.acc_db == acc_db]
        testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
        df_test = pd.read_csv(testdata_combined_file)

        X_test = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, df_test, s["excel_file_with_settings"])
        y_test = df_test["interface"].fillna(0).astype(int)

        #######################################################################################################
        #                                                                                                     #
        #                  Run prediction and save output individually for each protein                       #
        #                                                                                                     #
        #######################################################################################################

        prediction = forest.fit(X_train, y_train).predict_proba(X_test)[:, 1]
        # add the prediction to the combined file
        df_test[pred_colname] = prediction
        # save just the prediction alone to csv
        prediction_df = df_test[["residue_num", "residue_name", pred_colname]]
        prediction_df.to_csv(THOIPA_LOO_prediction_csv, index=False)

        fpr, tpr, thresholds = roc_curve(y_test, prediction, drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        xv_dict[acc_db] = {"fpr": fpr, "tpr": tpr, "roc_auc": roc_auc}
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0

        if database == "crystal" or database == "NMR":
            # (it is closest distance and low value means high propencity of interfacial)
            df_test["interface_score"] = -1 * df_test["interface_score"]

        BO_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, df_test, experiment_col="interface_score", pred_col=pred_colname)

        if BO_all_df.empty:
            BO_all_df = BO_df
        else:
            BO_all_df = pd.concat([BO_all_df, BO_df], axis=1, join="outer")
        logging.info("{} AUC : {:.2f}".format(acc_db, roc_auc))

        #######################################################################################################
        #                                                                                                     #
        #                          Get tree info, mean AUC for all proteins, etc                              #
        #                                                                                                     #
        #######################################################################################################

    tree_depths = np.array([estimator.tree_.max_depth for estimator in forest.estimators_])
    logging.info("tree depth mean = {} ({})".format(tree_depths.mean(), tree_depths))

    duration = time.clock() - start

    mean_tpr /= df_set.shape[0]
    mean_tpr[-1] = 1.0

    mean_roc_auc = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["mean_roc_auc"] = mean_roc_auc

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
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

    #linechart_mean_obs_and_rand = thoipapy.figs.Create_Bo_Curve_files.analyse_bo_curve_underlying_data(BO_all_data_csv, crossvalidation_folder, names_excel_path)
    thoipapy.figs.create_BOcurve_files.parse_BO_data_csv_to_excel(BO_all_data_csv, BO_data_excel, logging)

    logging.info('{} LOO crossvalidation. Time taken = {:.2f}.'.format(s["setname"], duration))
    logging.info('---AUC({:.2f})---'.format(mean_roc_auc))


def run_Rscipt_random_forest(s, output_file_loc, logging, Random_Forest_R_code_file):
    logging.info('begining to run machine learning R code')
    Rscript_loc = s["Rscript_dir"]
    #Random_Forest_R_code_file=s["Rcode"]
    train_data_file=os.path.join(s["thoipapy_data_folder"], "RandomForest","NoRedundPro/TRAINDATA68.csv")
    acc = s["tm_protein_name"]
    tmp_protein_test_data = os.path.join(s["thoipapy_data_folder"], "RandomForest", "TestData/%s/%s.mem.2gap.physipara.testdata.csv") % (s["Datatype"],acc)
    #out_put_file_loc_handle=open(output_file_loc,"w")
    if os.path.isfile(tmp_protein_test_data):
        prediction_output_file = os.path.join(s["thoipapy_data_folder"], "RandomForest","%s.pred.out") % acc
        prediction_output_file = os.path.join("/home/students/zeng/workspace/test2/out", "%s.pred.out") % acc
        prediction_output_file_handle=open(prediction_output_file,"w")
        exect_str = "{Rscript} {Random_Forest_R_code} {train_data} {test_data} {output}".format(Rscript=Rscript_loc, Random_Forest_R_code=Random_Forest_R_code_file,train_data=train_data_file,test_data=tmp_protein_test_data,output=output_file_loc)
        sys.stdout.write(exect_str)
        class Command(object):
            def __init__(self, cmd):
                self.cmd = cmd
                self.process = None

            def run(self, timeout):
                def target():
                    sys.stdout.write('Thread started')
                    self.process = subprocess.Popen(self.cmd,shell=True,stdout=subprocess.PIPE)
                    #subprocess.call(self.cmd,shell=True,stdout=prediction_output_file_handle)
                    subprocess.call(self.cmd, shell=True)
                    #self.process.communicate()
                    # sys.stdout.write(self.process.communicate())
                    sys.stdout.write('Thread finished')

                thread = threading.Thread(target=target)
                thread.start()

                thread.join(timeout)
                if thread.is_alive():
                    sys.stdout.write('Terminating process')
                    self.process.terminate()
                    thread.join()
                    sys.stdout.write(self.process.returncode)

        command = Command(exect_str)
        # command=Command(exect_str)
        command.run(timeout=2000)                       ###since hhblits requres more than 10 minutes to finish, maybe later we should consider using qsub to the server
        # command.run(timeout=1)
        # command=mtutils.Command(exect_str)
        # command.run(timeout=120)
        logging.info("Output file: %s\n" % output_file_loc)
        prediction_output_file_handle.close()


def predict_test_dataset_with_THOIPA_DEPRECATED(train_setname, test_setname, s, logging):
    """ Predict the interface of residues within one set (e.g. set03)
    with a trained model from another set (e.g. set04).

    The training set is currently the "run" set according to the settings file.
    The test set is an option in the settings file.

    IMPORTANT. CURRENTLY THERE IS NO REDUNDANCY CHECK.
     - the same proteins could be in the test and train datasets
     - homologues of the tested protein could be in the training dataset

    Parameters
    ----------
    train_setname : str
        Name of the dataset used for training. E.g. "set04".
    test_setname : str
        Name of the dataset used for testing. E.g. "set03".
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    THOIPA_pred_csv : csv
        csv file with prediction results
        columns = ["interface", "interface_score", "THOIPA", etc]
        rows = range(0, number of AA in test set)
    """

    model_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_ML_model.lpkl".format(train_setname))
    test_data_csv = os.path.join(s["Results_folder"], test_setname, "{}_train_data.csv".format(test_setname))
    THOIPA_pred_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "trainset{}_testset{}_predictions.csv".format(train_setname[-2:], test_setname[-2:]))

    fit = joblib.load(model_pkl)

    df_data = pd.read_csv(test_data_csv, index_col=0)
    df_testdata = drop_cols_not_used_in_ML(logging, df_data, s["excel_file_with_settings"])
    tX = df_testdata

    tp = fit.predict_proba(tX)

    df_out = df_data[["acc_db", "residue_num", "residue_name", "interface", "interface_score", "n_homologues"]]
    # if "interface_score" in df_data.columns:
    #     df_out["interface_score"] = df_data["interface_score"]
    #     df_out["interface_score"] = df_data["interface_score"]

    df_out["THOIPA"] = tp[:, 1]  # tools.normalise_0_1(tp[:, 1])[0]

    df_out.to_csv(THOIPA_pred_csv)
    logging.info('finished predict_test_dataset_with_THOIPA_DEPRECATED ({})'.format(THOIPA_pred_csv))


def create_one_out_train_data(acc_db,set_path,s):
    df_train = pd.DataFrame()
    df_set04 = pd.read_excel(set_path, sheet_name='proteins')
    for j in df_set04.index:
        acc1 = df_set04.loc[j, "acc_db"]
        if not acc1 == acc_db:
            database = df_set04.loc[j, "database"]
            feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                 "{}.surr20.gaps5.combined_features.csv".format(acc1))

            df_features_new_protein1 = pd.read_csv(feature_combined_file, index_col=0)
            df_features_new_protein1["acc_db"] = "{}-{}".format(acc1, database)

            # reorder the columns
            df_features_new_protein1 = thoipapy.utils.reorder_dataframe_columns(df_features_new_protein1,
                                                                                ['acc_db', 'residue_num', 'residue_name',
                                                                           'n_homologues'])
            # for the first protein, replace the empty dataframe
            if df_train.empty:
                df_train = df_features_new_protein1
            else:
                # concatenate the growing dataframe of combined proteins and new dataframe
                df_train = pd.concat([df_train, df_features_new_protein1])

                # reset the index to be a range (0,...).
    df_train.index = range(df_train.shape[0])
    return df_train


def calc_roc_each_tmd_separately_deprecated(s, df_set, logging):

    logging.info("start create_ROC_Curve_figs_THOIPA_PREDDIMER_TMDOCK_LIPS")
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])

    prediction_name_list = [pred_colname, "PREDDIMER", "TMDOCK", "LIPS_surface_ranked"]
    testsetname, trainsetname = get_testsetname_trainsetname_from_run_settings(s)
    if s["setname"] == testsetname:
        prediction_name_list.append(f"thoipa.train{trainsetname}")

    roc_each_tmd_separate_deprecated_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/{s['setname']}_roc_each_tmd_separate_deprecated.csv"
    roc_each_tmd_separate_deprecated_png = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/{s['setname']}_roc_each_tmd_separate_deprecated.png"

    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    mean_tpr_list=[]
    for n, predictor in enumerate(prediction_name_list):
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        for i in df_set.index:
            acc = df_set.loc[i, "acc"]
            database = df_set.loc[i, "database"]
            merged_data_csv_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/predictions/merged/{database}.{acc}.merged.csv"
            dfm = pd.read_csv(merged_data_csv_path, engine="python", index_col=0)
            dfm.dropna(inplace=True)
            interface = dfm["interface"].values
            if n ==0 or n == 3:
                predict = dfm[predictor].values
            else:
                predict = -1 * dfm[predictor].values
            fpr, tpr, thresholds = roc_curve(interface, predict, drop_intermediate=False)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0
        mean_tpr /= len(df_set.index)
        mean_tpr[-1] = 1.0
        mean_roc_auc = auc(mean_fpr, mean_tpr)
        mean_tpr_list.append(mean_tpr)
        ax.plot(mean_fpr, mean_tpr, lw=1,label="{} (area = {:.2f})".format(predictor, mean_roc_auc), alpha=0.8)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(roc_each_tmd_separate_deprecated_png, dpi=240)
    #fig.savefig(thoipapy.utils.pdf_subpath(ROC_4predictor_png))
    df_tpr = pd.DataFrame.from_records(list(map(list, zip(*mean_tpr_list))),
                                       columns=prediction_name_list)
    df_tpr.to_csv(roc_each_tmd_separate_deprecated_csv)
    logging.info(f"fig saved: {roc_each_tmd_separate_deprecated_png}")