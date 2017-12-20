import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import thoipapy
import sys
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import os
import subprocess, threading
from sklearn.externals import joblib
import pickle
import time
from korbinian.utils import convert_truelike_to_bool, convert_falselike_to_bool

# intersect function
def intersect(a, b):
     return list(set(a) & set(b))

def drop_cols_not_used_in_ML(df_data, s):
    """Remove columns not used in machine learning training or testing.

    This includes
        1) text and name columns (e.g. residue_name)
        2) bind columns (interface, interface_score).
        3) columns of non-normalised features, or features that we don't want to include anymore.

    Parameters
    ----------
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns = 'acc_db', "residue_num", "residue_name", etc
        rows = range of number of residues

    Returns
    -------
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns :
    """

    # columns in combined file that are not used in machine-learning training or testing
    # cols_to_drop = ["acc_db", "residue_num", "residue_name", "interface", "interface_score", "Entropy", "n_homologues", "LIPS_entropy", "LIPS_L*E", "LIPS_surface_ranked",
    #                 "CoevDImax", "CoevDI4", "CoevDI8",  "CoevMImax", "CoevMI4", "CoevMI8", "CumDI4", "CumDI8","CoevDI4_norm", "CoevMI4_norm",
    #                 "CumDI4_norm", "CumDI8_norm", "CumMI4_norm", "CumMI8_norm"]#, "RelPos_TMD", "RelPos_fullseq"
    #                 #"CoevDImax_norm", "CoevDI8_norm", "CoevMImax_norm", "CoevMI8_norm", ]#

    # cols_to_drop = ["acc_db", "residue_num", "residue_name", "interface", "interface_score", "Entropy",
    #                 #"CumDI4", "CumDI8", "CumMI4", "CumMI8", "CumDI4_norm", "CumDI8_norm", "CumMI4_norm", "CumMI8_norm",
    #                 #"Aromatic_sAA", "Polar_sAA", "Cbbranched_sAA", "Small_sAA",
    #                 #"H", "E", "K", "D", "N",
    #                 #"C", "S", "D", "E", "K", "R", "Q", "N",
    #                 "CS", "DE", "KR", "QN",
    #                 "RelPos_TMD", "RelPos_fullseq",
    #                 "CoevDImax", "CoevDI4", "CoevDI8",  "CoevMImax", "CoevMI4", "CoevMI8",
    #                 #"CoevDI4_norm", "CoevMI4_norm",
    #                 "CoevDImax_norm", 'CoevMImax_norm',
    #                 #"LIPS_entropy", "LIPS_L*E", "LIPS_lipo",
    #                 #"LIPS_surface",
    #                 #"LIPS_surface_ranked",
    #                 #"LIPS_surface_ranked_norm",
    #                 #"n_homologues",
    #                 "n_homol_norm",
    #                 #"GxxxG",
    #                 #"SmxxxSm",
    #                 #"PolarxxxPolar"
    #                 ]
    #                 #"n_homologues", , "LIPS_surface_ranked",
    #                 #"
    #                 #"CumDI4_norm", "CumDI8_norm", "CumMI4_norm", "CumMI8_norm"]#, "RelPos_TMD", "RelPos_fullseq"
    #                 #"CoevDImax_norm", "CoevDI8_norm", "CoevMImax_norm", "CoevMI8_norm", ]#

    # pd.Series(cols_to_drop).to_csv(r"D:\data_thoipapy\Results\set02\set02_dropped_features.csv")
    #
    #
    # # got only those that are actually in the columns
    # cols_to_drop = set(cols_to_drop).intersection(set(df_data.columns))
    # df_data = df_data.drop(cols_to_drop, axis=1)

    # read the features tab of the excel settings file
    features_df = pd.read_excel(s["excel_file_with_settings"], sheetname="features")
    # convert "WAHR" etc to true and false
    features_df["include"] = features_df["include"].apply(convert_truelike_to_bool, convert_nontrue=False)
    features_df["include"] = features_df["include"].apply(convert_falselike_to_bool)
    # drop any features that are not labeled TRUE for inclusion
    features_df = features_df.loc[features_df.include == True]
    # filter df_data to only keep the desired feature columns
    feature_list = features_df.feature.tolist()
    df_data = df_data.loc[:, feature_list]

    return df_data

def THOIPA_RF_classifier_with_settings(s):
    """ For tuning the RF parameters, they are always in one place, and determined by the settings file.

    Parameters
    ----------
    s

    Returns
    -------

    """

    # convert max_features to python None if "None"
    max_features = None if s["max_features"] == "None" else s["max_features"]

    forest = RandomForestClassifier(n_estimators=s["RF_number_of_estimators"], n_jobs=s["n_CPU_cores"], criterion=s["criterion"],
                                    min_samples_leaf=s["min_samples_leaf"],
                                    #max_depth=s["max_depth"],
                                    oob_score=True, max_features=max_features, bootstrap=bool(s["bootstrap"]),
                                    #random_state=s["random_state"]
                                    )
    return forest

def train_random_forest_model(s, logging):
    """Train the random forest model for a particular set.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Sawed Files
    -----------
    model_pkl : pickle
        Pickle containing the trained random forest model.

    """
    logging.info('starting to predict etra data with THOIPA prediction model')

    train_data_csv = os.path.join(s["set_results_folder"], "{}_train_data.csv".format(s["setname"]))
    train_data_used_for_model_csv = os.path.join(s["set_results_folder"], "{}_train_data_used_for_model.csv".format(s["setname"]))
    model_pkl = os.path.join(s["set_results_folder"], "{}_rfmodel.pkl".format(s["setname"]))

    df_data = pd.read_csv(train_data_csv, index_col=0)

    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    y = df_data["interface"]
    X = drop_cols_not_used_in_ML(df_data, s)

    X.to_csv(train_data_used_for_model_csv)

    if 1 not in y.tolist():
        raise ValueError("None of the residues are marked 1 for an interface residue!")

    forest = THOIPA_RF_classifier_with_settings(s)
    # save random forest model into local driver
    # pkl_file = r'D:\thoipapy\RandomForest\rfmodel.pkl'
    fit = forest.fit(X, y)
    joblib.dump(fit, model_pkl)

    logging.info('finished training random forest algorithm ({})'.format(model_pkl))

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

    model_pkl = os.path.join(s["set_results_folder"], "{}_rfmodel.pkl".format(train_setname))
    test_data_csv = os.path.join(s["Results_folder"], test_setname, "{}_train_data.csv".format(test_setname))
    THOIPA_pred_csv = os.path.join(s["set_results_folder"], "trainset{}_testset{}_predictions.csv".format(train_setname[-2:], test_setname[-2:]))

    fit = joblib.load(model_pkl)

    df_data = pd.read_csv(test_data_csv, index_col=0)
    df_testdata = drop_cols_not_used_in_ML(df_data, s)
    tX = df_testdata

    tp = fit.predict_proba(tX)

    df_out = df_data[["acc_db", "residue_num", "residue_name", "interface", "interface_score", "n_homologues"]]
    # if "interface_score" in df_data.columns:
    #     df_out["interface_score"] = df_data["interface_score"]
    #     df_out["interface_score"] = df_data["interface_score"]

    df_out["THOIPA"] = tp[:, 1]  # tools.normalise_0_1(tp[:, 1])[0]

    df_out.to_csv(THOIPA_pred_csv)
    logging.info('finished predict_test_dataset_with_THOIPA_DEPRECATED ({})'.format(THOIPA_pred_csv))

def run_Rscipt_random_forest(s, output_file_loc, logging):
    logging.info('begining to run random forest R code')
    Rscript_loc = s["Rscript_dir"]
    Random_Forest_R_code_file=s["Rcode"]
    train_data_file=os.path.join(s["RF_loc"],"NoRedundPro/TRAINDATA68.csv")
    acc = s["tm_protein_name"]
    tmp_protein_test_data = os.path.join(s["RF_loc"], "TestData/%s/%s.mem.2gap.physipara.testdata.csv") % (s["Datatype"],acc)
    #out_put_file_loc_handle=open(output_file_loc,"w")
    if os.path.isfile(tmp_protein_test_data):
        prediction_output_file = os.path.join(s["RF_loc"],"%s.pred.out") % acc
        prediction_output_file = os.path.join("/home/students/zeng/workspace/test2/out", "%s.pred.out") % acc
        prediction_output_file_handle=open(prediction_output_file,"w")
        exect_str = "{Rscript} {Random_Forest_R_code} {train_data} {test_data} {output}".format(Rscript=Rscript_loc, Random_Forest_R_code=Random_Forest_R_code_file,train_data=train_data_file,test_data=tmp_protein_test_data,output=output_file_loc)
        print(exect_str)
        class Command(object):
            def __init__(self, cmd):
                self.cmd = cmd
                self.process = None

            def run(self, timeout):
                def target():
                    print('Thread started')
                    self.process = subprocess.Popen(self.cmd,shell=True,stdout=subprocess.PIPE)
                    #subprocess.call(self.cmd,shell=True,stdout=prediction_output_file_handle)
                    subprocess.call(self.cmd, shell=True)
                    #self.process.communicate()
                    # print(self.process.communicate())
                    print('Thread finished')

                thread = threading.Thread(target=target)
                thread.start()

                thread.join(timeout)
                if thread.is_alive():
                    print('Terminating process')
                    self.process.terminate()
                    thread.join()
                print(self.process.returncode)

        command = Command(exect_str)
        # command=Command(exect_str)
        command.run(timeout=2000)                       ###since hhblits requres more than 10 minutes to finish, maybe later we should consider using qsub to the server
        # command.run(timeout=1)
        # command=mtutils.Command(exect_str)
        # command.run(timeout=120)
        logging.info("Output file: %s\n" % output_file_loc)
        prediction_output_file_handle.close()



def run_10fold_cross_validation(s, logging):
    """Run 10-fold cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.

    The number of folds is determined by "cross_validation_number_of_splits" in the settings file.

    The number of estimators in the random forest algorithm is determined by "RF_number_of_estimators" in the settings file,
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
    train_data_csv = os.path.join(s["set_results_folder"], "{}_train_data.csv".format(s["setname"]))
    #crossvalidation_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_10F_data.csv".format(s["setname"]))
    crossvalidation_pkl = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_10F_data.pkl".format(s["setname"]))
    features_csv = os.path.join(s["set_results_folder"], "{}_test_features.csv".format(s["setname"]))

    thoipapy.utils.make_sure_path_exists(crossvalidation_pkl, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    X = drop_cols_not_used_in_ML(df_data, s)
    y = df_data["interface"]

    skf = StratifiedKFold(n_splits=s["cross_validation_number_of_splits"])
    cv = list(skf.split(X, y))

    forest = THOIPA_RF_classifier_with_settings(s)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    # save all outputs to a cross-validation dictionary, to be saved as a pickle file
    xv_dict = {}

    start = time.clock()

    for i, (train, test) in enumerate(cv):
        sys.stdout.write("f{}.".format(i+1)), sys.stdout.flush()
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1])
        xv_dict["fpr{}".format(i)] = fpr
        xv_dict["tpr{}".format(i)] = tpr
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    sys.stdout.write("\n"), sys.stdout.flush()

    logging.info("tree depths : {}".format([estimator.tree_.max_depth for estimator in forest.estimators_]))

    duration = time.clock() - start

    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0

    mean_auc = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["mean_auc"] = mean_auc

    # save dict as pickle
    with open(crossvalidation_pkl, "wb") as f:
        pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    features_ser = pd.Series(X.columns)
    features_ser.to_csv(features_csv)
    logging.info('{} 10-fold validation. AUC({:.3f}). Time taken = {:.2f}.\nFeatures: {}'.format(s["setname"], mean_auc, duration, X.columns.tolist()))

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
    crossvalidation_png = os.path.join(s["set_results_folder"], "crossvalidation", "{}_10F_ROC.png".format(s["setname"]))
    crossvalidation_pkl = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_10F_data.pkl".format(s["setname"]))

    # open pickle file
    with open(crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for i in range(s["cross_validation_number_of_splits"]):
        roc_auc = auc(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)])
        ax.plot(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)], lw=1, label='fold %d (area = %0.2f)' % (i, roc_auc), alpha=0.8)

    mean_auc = xv_dict["mean_auc"]

    ax.plot(xv_dict["false_positive_rate_mean"], xv_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_auc, lw=1.5)
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
    """Run Leave-One-Out cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.
    Each protein is uniquely identified by the acc and the database (acc_db).
    The test and training datasets are based on the following:
        1) the df_set derived from the set of protein sequences, e.g. set03
            - created manually
        2) the train_data csv, e.g."D:\data_thoipapy\Results\set03\set03_train_data.csv".
            - filtered according to redundancy etc by combine_all_train_data_for_random_forest()
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
    names_excel_path = os.path.join(os.path.dirname(s["sets_folder"]), "ETRA_NMR_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    train_data_csv = os.path.join(s["set_results_folder"], "{}_train_data.csv".format(s["setname"]))
    crossvalidation_folder = os.path.join(s["set_results_folder"], "crossvalidation")
    LOO_crossvalidation_pkl = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))
    BO_curve_folder = os.path.join(s["set_results_folder"], "crossvalidation")
    BO_data_excel = os.path.join(BO_curve_folder, "data", "BO_curve_data.xlsx")
    BO_linechart_png = os.path.join(BO_curve_folder, "BO_linechart.png")
    BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")
    thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

    df_data = pd.read_csv(train_data_csv)

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    acc_db_list = df_data.acc_db.unique()
    xv_dict = {}
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    start = time.clock()
    BO_all_df = pd.DataFrame()

    for i in df_set.index:
        acc, acc_db, database  = df_set.loc[i, "acc"], df_set.loc[i, "acc_db"], df_set.loc[i, "database"]
        prediction_csv = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.LOO.prediction.csv".format(acc, s["setname"]))
        thoipapy.utils.make_sure_path_exists(prediction_csv, isfile=True)

        if not acc_db in acc_db_list:
            logging.warning("{} is in protein set, but not found in training data".format(acc_db))
            # skip protein
            continue
        df_train = df_data.loc[df_data.acc_db != acc_db]
        # df_train = df_data.loc[df_data.interface_score < 5.5]
        df_test = df_data.loc[df_data.acc_db == acc_db]
        X_train = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(df_train, s)
        y_train = df_train["interface"]
        X_test = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(df_test, s)
        y_test = df_test["interface"]
        forest = thoipapy.RF_features.RF_Train_Test.THOIPA_RF_classifier_with_settings(s)
        prediction = forest.fit(X_train, y_train).predict_proba(X_test)[:, 1]

        # add the prediction to the combined file
        pred_colname = "THOIPA_{}_LOO".format(s["set_number"])
        df_test[pred_colname] = prediction
        # save just the prediction alone to csv
        prediction_df = df_test[["residue_num", "residue_name", pred_colname]]
        prediction_df.to_csv(prediction_csv, index=False)

        fpr, tpr, thresholds = roc_curve(y_test, prediction)
        roc_auc = auc(fpr, tpr)
        xv_dict[acc_db] = {"fpr": fpr, "tpr": tpr, "auc": roc_auc}
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

    logging.info([estimator.tree_.max_depth for estimator in forest.estimators_])

    duration = time.clock() - start

    mean_tpr /= df_set.shape[0]
    mean_tpr[-1] = 1.0

    mean_auc = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["mean_auc"] = mean_auc

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
    names_excel_path = os.path.join(os.path.dirname(s["sets_folder"]), "ETRA_NMR_names.xlsx")

    #linechart_mean_obs_and_rand = thoipapy.figs.Create_Bo_Curve_files.analyse_bo_curve_underlying_data(BO_all_data_csv, crossvalidation_folder, names_excel_path)
    thoipapy.figs.Create_Bo_Curve_files.parse_BO_data_csv_to_excel(BO_all_data_csv, BO_data_excel, logging)
    AUBOC10 = thoipapy.figs.Create_Bo_Curve_files.save_BO_linegraph_and_barchart(BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging)
    if "you_want_more_details" == "TRUE":
        other_figs_path = os.path.join(BO_curve_folder, "other_figs")
        thoipapy.figs.Create_Bo_Curve_files.save_extra_BO_figs(BO_data_excel, other_figs_path)

    #sys.stdout.write("\nBO curve data analysed ({})".format(linechart_mean_obs_and_rand))

    logging.info('{} LOO crossvalidation. AUC({:.2f}), AUBOC10({:.2f}). Time taken = {:.2f}.'.format(s["setname"], mean_auc, AUBOC10, duration))


def create_LOO_validation_fig(s, df_set, logging):
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
    # drop redundant proteins according to CD-HIT
    df_set = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    #plt.rcParams.update({'font.size': 7})
    LOO_crossvalidation_pkl = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_LOO_crossvalidation.pkl".format(s["setname"]))
    LOO_crossvalidation_ROC_png = os.path.join(s["set_results_folder"], "crossvalidation", "{}_LOO_crossvalidation_ROC.png".format(s["setname"]))
    LOO_crossvalidation_AUC_bar_png = os.path.join(s["set_results_folder"], "crossvalidation", "{}_LOO_crossvalidation_AUC_bar.png".format(s["setname"]))
    AUC_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_LOO_AUC.csv".format(s["setname"]))
    BO_all_data_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "{}_LOO_BO_data.csv".format(s["setname"]))

    # open pickle file
    with open(LOO_crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    # due to problems on Bo's computer, set the figsize to double what we should be using for the publication?
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    auc_dict = {}

    for acc_db in df_set.acc_db:
        if acc_db in xv_dict:
            roc_auc = xv_dict[acc_db]["auc"]
            auc_dict[acc_db] = roc_auc
            ax.plot(xv_dict[acc_db]["fpr"], xv_dict[acc_db]["tpr"], lw=1, label='{} ({:.2f})'.format(acc_db, roc_auc), alpha=0.8)
        else:
            logging.warning("{} not in xv_dict after LOO validation".format(acc_db))

    mean_auc = xv_dict["mean_auc"]

    ax.plot(xv_dict["false_positive_rate_mean"], xv_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_auc, lw=1.5)
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
    AUC_ser.sort_values(inplace=True, ascending=False)
    AUC_ser.to_csv(AUC_csv)
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    AUC_ser.plot(ax=ax, kind="bar")
    ax.set_ylabel("performance (AUBOC10)")
    fig.tight_layout()
    fig.savefig(LOO_crossvalidation_AUC_bar_png, dpi=240)
    fig.savefig(thoipapy.utils.pdf_subpath(LOO_crossvalidation_AUC_bar_png))

    logging.info("create_LOO_validation_fig finished ({})".format(LOO_crossvalidation_AUC_bar_png))

    #
    #
    # plt.close("all")
    # fig, ax = plt.subplots()
    # ax2 = ax.twinx()
    # df_o_minus_r.T.mean().plot(ax=ax, color="#0f7d9b", linestyle="-", label="new method (o-r)")
    # ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="new method random", alpha=0.5)
    #
    # df_o_over_r.T.mean().plot(ax=ax2, color="#9b2d0f", linestyle="-", label="old method (o/r)")
    # ax2.plot([1, 10], [1, 1], color="#9b2d0f", linestyle="--", label="old method random", alpha=0.5)
    #
    # # ax.set_ylim(0)
    # ax.grid(False)
    # ax.set_ylabel("performance value\n(observed - random)", color="#0f7d9b")
    # ax2.set_ylabel("performance value\n (observed / random)", color="#9b2d0f")
    # ax.tick_params('y', colors="#0f7d9b")
    # ax2.tick_params('y', colors="#9b2d0f")
    # ax2.spines['left'].set_color("#0f7d9b")
    # ax2.spines['right'].set_color("#9b2d0f")
    #
    # ax.legend()
    # ax2.legend()
    # # fig.legend()
    # fig.tight_layout()
    # fig.savefig(linechart_BO_curve_single_dataset, dpi=140)


def calculate_RF_variable_importance(s, logging):
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
        Also includes the standard deviation supplied by the random forest algorithm
    """
    #logging.info('RF_variable_importance_calculate is running\n')
    train_data_csv = os.path.join(s["set_results_folder"], "{}_train_data.csv".format(s["setname"]))

    variable_importance_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "trainset{}_testset{}_variable_importance.csv".format(s["setname"][-2:], s["setname"][-2:]))
    thoipapy.utils.make_sure_path_exists(variable_importance_csv, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    X = drop_cols_not_used_in_ML(df_data, s)
    y = df_data["interface"]
    forest = THOIPA_RF_classifier_with_settings(s)
    forest.fit(X, y)
    importances_arr = forest.feature_importances_
    std_arr = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    indices_arr = np.argsort(importances_arr)[::-1]

    importances_text_list = X.columns.tolist()

    order_list = [importances_arr[indices_arr[f]] for f in range(X.shape[1])]

    logging.info("\nFeature ranking:")

    nested_dict = {}

    for f in range(X.shape[1]):
        if f < 10:
            logging.info("%d. feature %d (%f) %s" % (f + 1, indices_arr[f], importances_arr[indices_arr[f]], importances_text_list[indices_arr[f]]))
        single_feature_dict = {"feature_order" : indices_arr[f], "mean_decrease_gini" : importances_arr[indices_arr[f]], "feature" : importances_text_list[indices_arr[f]],  "std" : std_arr[f]}
        nested_dict[f + 1] = single_feature_dict

    sys.stdout.write("\n\n"), sys.stdout.flush()

    df_imp = pd.DataFrame(nested_dict).T
    df_imp["orig_order"] = df_imp.index
    df_imp.set_index("feature", inplace=True)

    df_imp.to_csv(variable_importance_csv)

def fig_variable_importance(s, logging):
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
    plt.rcParams.update({'font.size': 8})
    from korbinian.utils import create_colour_lists
    colour_dict = create_colour_lists()

    variable_importance_csv = os.path.join(s["set_results_folder"], "crossvalidation", "data", "trainset{}_testset{}_variable_importance.csv".format(s["setname"][-2:], s["setname"][-2:]))
    variable_importance_all_png = os.path.join(s["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance_all.png".format(s["setname"][-2:], s["setname"][-2:]))
    variable_importance_top_png = os.path.join(s["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance_top.png".format(s["setname"][-2:], s["setname"][-2:]))

    df_imp = pd.read_csv(variable_importance_csv, index_col = 0)

    create_var_imp_plot(df_imp, colour_dict, variable_importance_all_png, df_imp.shape[0])
    create_var_imp_plot(df_imp, colour_dict, variable_importance_top_png, 30)


def create_var_imp_plot(df_imp, colour_dict, variable_importance_png, n_features_in_plot):
    """Plot function for fig_variable_importance, allowing a variable number of features.
    """
    df_sel = df_imp.iloc[:n_features_in_plot, :].copy()
    df_sel.sort_values("mean_decrease_gini", ascending=True, inplace=True)
    min_ = df_sel.mean_decrease_gini.min()

    # determine the plot height by the number of features
    # currently set for 30
    plot_height = 4 * n_features_in_plot / 30
    figsize = np.array([3.42, plot_height])
    fig, ax = plt.subplots(figsize=figsize)

    TUMblue = colour_dict["TUM_colours"]['TUMBlue']
    df_sel.mean_decrease_gini.plot(kind="barh", color=TUMblue, ax=ax)# xerr=df_sel["std"]
    ax.errorbar(df_sel.mean_decrease_gini, range(len(df_sel.index)), xerr=df_sel["std"], fmt="none", ecolor="k", ls="none", capthick=0.5, elinewidth=0.5, capsize=1, label=None)

    ax.set_ylabel("")
    ax.set_xlabel("variable importance\n(mean decrease gini)")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(variable_importance_png, dpi=240)
    #fig.savefig(variable_importance_png[:-4] + ".pdf")
    fig.savefig(thoipapy.utils.pdf_subpath(variable_importance_png))

