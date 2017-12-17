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
# import sys
# import glob
import numpy as np
# import eccpy.tools as tools
from sklearn.externals import joblib
import pickle
import time

# intersect function
def intersect(a, b):
     return list(set(a) & set(b))

def drop_cols_not_used_in_ML(df_data):
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


    cols_to_drop = ["acc_db", "residue_num", "residue_name", "interface", "interface_score", "Entropy",
                    #"CumDI4", "CumDI8", "CumMI4", "CumMI8", "CumDI4_norm", "CumDI8_norm", "CumMI4_norm", "CumMI8_norm",
                    #"Aromatic_sAA", "Polar_sAA", "Cbbranched_sAA", "Small_sAA",
                    #"H", "E", "K", "D", "N",
                    #"C", "S", "D", "E", "K", "R", "Q", "N",
                    "CS", "DE", "KR", "QN",
                    "RelPos_TMD", "RelPos_fullseq",
                    "CoevDImax", "CoevDI4", "CoevDI8",  "CoevMImax", "CoevMI4", "CoevMI8", "CoevDI4_norm", "CoevMI4_norm",
                    "CoevDImax_norm", 'CoevMImax_norm',
                    #"LIPS_entropy", "LIPS_L*E", "LIPS_lipo",
                    "LIPS_surface",
                    "LIPS_surface_ranked",
                    #"n_homologues",
                    "n_homol_norm"
                    ]
                    #"n_homologues", , "LIPS_surface_ranked",
                    #"
                    #"CumDI4_norm", "CumDI8_norm", "CumMI4_norm", "CumMI8_norm"]#, "RelPos_TMD", "RelPos_fullseq"
                    #"CoevDImax_norm", "CoevDI8_norm", "CoevMImax_norm", "CoevMI8_norm", ]#

    # got only those that are actually in the columns
    cols_to_drop = set(cols_to_drop).intersection(set(df_data.columns))
    df_data = df_data.drop(cols_to_drop, axis=1)
    return df_data

def THOIPA_RF_classifier_with_settings(set_):
    """ For tuning the RF parameters, they are always in one place, and determined by the settings file.

    Parameters
    ----------
    set_

    Returns
    -------

    """

    # convert max_features to python None if "None"
    max_features = None if set_["max_features"] == "None" else set_["max_features"]

    forest = RandomForestClassifier(n_estimators=set_["RF_number_of_estimators"], n_jobs=set_["n_jobs"], criterion=set_["criterion"],
                                    min_samples_leaf=set_["min_samples_leaf"],
                                    #max_depth=set_["max_depth"],
                                    oob_score=True, max_features=max_features, bootstrap=bool(set_["bootstrap"]),
                                    #random_state=set_["random_state"]
                                    )
    return forest

def train_random_forest_model(set_, logging):
    """Train the random forest model for a particular set.

    Parameters
    ----------
    set_ : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Sawed Files
    -----------
    model_pkl : pickle
        Pickle containing the trained random forest model.

    """
    logging.info('starting to predict etra data with THOIPA prediction model')

    train_data_csv = os.path.join(set_["set_results_folder"], "{}_train_data.csv".format(set_["setname"]))
    train_data_used_for_model_csv = train_data_csv[:-4] + "used_for_model.csv"
    model_pkl = os.path.join(set_["set_results_folder"], "{}_rfmodel.pkl".format(set_["setname"]))

    df_data = pd.read_csv(train_data_csv, index_col=0)

    df_data = df_data.loc[df_data.n_homologues >= set_["min_n_homol_training"]]

    y = df_data["interface"]
    X = drop_cols_not_used_in_ML(df_data)

    X.to_csv(train_data_used_for_model_csv)

    if 1 not in y.tolist():
        raise ValueError("None of the residues are marked 1 for an interface residue!")

    forest = THOIPA_RF_classifier_with_settings(set_)
    # save random forest model into local driver
    # pkl_file = r'D:\thoipapy\RandomForest\rfmodel.pkl'
    fit = forest.fit(X, y)
    joblib.dump(fit, model_pkl)

    logging.info('finished training random forest algorithm ({})'.format(model_pkl))

def predict_test_dataset_with_THOIPA(train_setname, test_setname, set_, logging):
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
    set_ : dict
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

    model_pkl = os.path.join(set_["set_results_folder"], "{}_rfmodel.pkl".format(train_setname))
    test_data_csv = os.path.join(set_["Results_folder"], test_setname, "{}_train_data.csv".format(test_setname))
    THOIPA_pred_csv = os.path.join(set_["set_results_folder"], "trainset{}_testset{}_predictions.csv".format(train_setname[-2:], test_setname[-2:]))

    fit = joblib.load(model_pkl)

    df_data = pd.read_csv(test_data_csv, index_col=0)
    df_testdata = drop_cols_not_used_in_ML(df_data)
    tX = df_testdata

    tp = fit.predict_proba(tX)

    df_out = df_data[["acc_db", "residue_num", "residue_name", "interface", "interface_score", "n_homologues"]]
    # if "interface_score" in df_data.columns:
    #     df_out["interface_score"] = df_data["interface_score"]
    #     df_out["interface_score"] = df_data["interface_score"]

    df_out["THOIPA"] = tp[:, 1]  # tools.normalise_0_1(tp[:, 1])[0]

    df_out.to_csv(THOIPA_pred_csv)
    logging.info('finished predict_test_dataset_with_THOIPA ({})'.format(THOIPA_pred_csv))

    # # fit = joblib.load(pkl_file)
    #
    # # test etra data
    # testdata_list = glob.glob(os.path.join(set_["thoipapy_data_folder"],"Features", "combined/etra", "*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))
    #
    #
    # i = 0
    # for test_data in testdata_list:
    #     acc = test_data.split('\\')[-1][0:6]
    #     # if acc == "O75460":
    #     dirupt_path = os.path.join(set_["base_dir"],"data_xy","Figure","Show_interface","Interface_xlsx", "{}.xlsx".format(acc))
    #     ddf = pd.read_excel(dirupt_path, index_col=0)
    #     disruption = ddf.Disruption
    #     thoipa_out = os.path.join(set_["thoipapy_data_folder"],"Features","combined/etra", "{}.thoipa_pred.csv".format(acc))
    #     tdf = pd.read_csv(test_data, sep=',', engine='python', index_col=0)
    #     tdf.index = tdf.index.astype(int) + 1
    #     aa = tdf.residue_name
    #
    #     coev_colname_list = ["CoevDImax", "CoevDI4", "CoevDI8", "CumDI4", "CumDI8", "CoevMImax", "CoevMI4", "CoevMI8", "CumMI4", "CumMI8"]
    #     list_cols_not_used_in_ML = ['acc_db', "residue_num", "residue_name",  "n_homologues", "Entropy"] + coev_colname_list
    #
    #     tdf = tdf.drop(list_cols_not_used_in_ML, axis = 1)
    #     tX = tdf[tdf.columns]
    #     tp = fit.predict_proba(tX)
    #     odf = pd.DataFrame()
    #     odf["AA"] = aa
    #     odf["thoipa"] = tools.normalise_0_1(tp[:, 1])[0]
    #     odf["disruption"] = tools.normalise_0_1(disruption)[0]
    #     odf.to_csv(thoipa_out)



def run_Rscipt_random_forest(set_, output_file_loc, logging):
    logging.info('begining to run random forest R code')
    Rscript_loc = set_["Rscript_dir"]
    Random_Forest_R_code_file=set_["Rcode"]
    train_data_file=os.path.join(set_["RF_loc"],"NoRedundPro/TRAINDATA68.csv")
    acc = set_["tm_protein_name"]
    tmp_protein_test_data = os.path.join(set_["RF_loc"], "TestData/%s/%s.mem.2gap.physipara.testdata.csv") % (set_["Datatype"],acc)
    #out_put_file_loc_handle=open(output_file_loc,"w")
    if os.path.isfile(tmp_protein_test_data):
        prediction_output_file = os.path.join(set_["RF_loc"],"%s.pred.out") % acc
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



def run_10fold_cross_validation(set_, logging):
    """Run 10-fold cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.

    The number of folds is determined by "cross_validation_number_of_splits" in the settings file.

    The number of estimators in the random forest algorithm is determined by "RF_number_of_estimators" in the settings file,
    therefore it should match the full training result of the whole dataset.

    IMPORTANT. CURRENTLY THERE IS NO AUTOMATIC REDUNDANCY CHECK.
     - homologues of the tested protein could be in the training dataset

    Parameters
    ----------
    set_ : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    crossvalidation_pkl : pickle
        Pickled dictionary (xv_dict) containing the results for each fold of validation.
        Also contains the mean ROC curve, and the mean AUC.
    """
    logging.info('10-fold cross validation is running\n')
    train_data_csv = os.path.join(set_["set_results_folder"], "{}_train_data.csv".format(set_["setname"]))
    crossvalidation_csv = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_crossvalidation.csv".format(set_["setname"][-2:], set_["setname"][-2:]))
    crossvalidation_pkl = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_crossvalidation.pkl".format(set_["setname"][-2:], set_["setname"][-2:]))

    thoipapy.utils.make_sure_path_exists(crossvalidation_csv, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)

    # drop training data (full protein) that don't have enough homologues
    df_data = df_data.loc[df_data.n_homologues >= set_["min_n_homol_training"]]

    #data = pd.read_csv('/scratch2/zeng/homotypic_data/data/RandomForest/PsEnCo/TrainData2',delimiter="\s",engine='python')
    # del data["Residue_id"]
    # del data["Residue_name"]
    #print(data.as_matrix(data.columns))
    # features=data.columns[0:28]
    # X=data[features]
    # y=data["Bind"]
    X = drop_cols_not_used_in_ML(df_data)
    y = df_data["interface"]
    #n_samples, n_features = X.shape
    #random_state = np.random.RandomState(0)
    #X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]
    #print(X.iloc[[20,22]])
    #StratifiedKFold(n_splits=2, random_state=None, shuffle=False)
    #cv = StratifiedKFold(y, n_folds=6)

    skf = StratifiedKFold(n_splits=set_["cross_validation_number_of_splits"])
    cv = list(skf.split(X, y))

    # # convert max_features to python None if "None"
    # max_features = None if set_["max_features"] == "None" else set_["max_features"]
    #
    # forest = RandomForestClassifier(n_estimators=set_["RF_number_of_estimators"], n_jobs=set_["n_jobs"], criterion=set_["criterion"],
    #                                 min_samples_leaf=set_["min_samples_leaf"],
    #                                 #max_depth=set_["max_depth"],
    #                                 oob_score=True, max_features=max_features, bootstrap=bool(set_["bootstrap"]),
    #                                 #random_state=set_["random_state"]
    #                                 )

    forest = THOIPA_RF_classifier_with_settings(set_)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []
    df_xv = pd.DataFrame()
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
        roc_auc = auc(fpr, tpr)
        #plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

    print([estimator.tree_.max_depth for estimator in forest.estimators_])

    duration = time.clock() - start
    logging.info("time taken for cross validation = {:.2f} (s in Windows, min in Linux?)".format(duration))

    #plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0

    mean_auc = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["mean_auc"] = mean_auc

    # save dict as pickle
    with open(crossvalidation_pkl, "wb") as f:
        pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    # df_xv.loc["mean_auc"] = mean_auc
    # df_xv.to_csv(crossvalidation_csv)
    logging.info('10-fold cross validation is finished. Mean AUC = {:.3f}\nFeatures included:\n{}'.format(mean_auc, X.columns.tolist()))

def fig_10fold_cross_validation(set_, logging):
    """Create figure showing ROC curve for each fold in a 10-fold validation.

    The underlying data is created by run_10fold_cross_validation. If this has not been run,
    it will return a file-not-found error.

    Parameters
    ----------
    set_ : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    plt.rcParams.update({'font.size': 7})
    crossvalidation_csv = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_crossvalidation.csv".format(set_["setname"][-2:], set_["setname"][-2:]))
    crossvalidation_png = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_ROC.png".format(set_["setname"][-2:], set_["setname"][-2:]))
    crossvalidation_pkl = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_crossvalidation.pkl".format(set_["setname"][-2:], set_["setname"][-2:]))

    # open pickle file
    with open(crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    for i in range(set_["cross_validation_number_of_splits"]):
        roc_auc = auc(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)])
        ax.plot(xv_dict["fpr{}".format(i)], xv_dict["tpr{}".format(i)], lw=1, label='fold %d (area = %0.2f)' % (i, roc_auc), alpha=0.8)

    #mean_auc = auc(df_xv["false_positive_rate"], df_xv["true_positive_rate"])
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
    fig.savefig(crossvalidation_png[:-4] + ".pdf")

    # df_xv = pd.read_csv(crossvalidation_csv)
    #
    # #mean_auc = auc(df_xv["false_positive_rate"], df_xv["true_positive_rate"])
    # mean_auc = df_xv.loc[0, "mean_auc"]
    #
    # fig, ax = plt.subplots(figsize=(4.2, 4.2))
    #
    # ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    # ax.plot(xv_dict["false_positive_rate_mean"], xv_dict["true_positive_rate_mean"], 'k--',
    #          label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
    # ax.set_xlim([-0.05, 1.05])
    # ax.set_ylim([-0.05, 1.05])
    # ax.set_xlabel("False positive rate")
    # ax.set_ylabel("True positive rate")
    # ax.legend(loc="lower_right")
    # fig.tight_layout()
    # fig.savefig(crossvalidation_png)

    # plt.plot(df_xv["false_positive_rate"], df_xv["true_positive_rate"], 'k--',
    #          label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
    # plt.xlim([-0.05, 1.05])
    # plt.ylim([-0.05, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('Receiver operating characteristic example')
    # plt.legend(loc="lower right")
    # plt.savefig(crossvalidation_png)
    #plt.show()
    # ADD SAVED PLOT

def calculate_RF_variable_importance(set_, logging):
    """Calculate the variable importance (mean decrease gini) for all variables in THOIPA.

    Parameters
    ----------
    set_ : dict
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
    train_data_csv = os.path.join(set_["set_results_folder"], "{}_train_data.csv".format(set_["setname"]))

    variable_importance_csv = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance.csv".format(set_["setname"][-2:], set_["setname"][-2:]))
    thoipapy.utils.make_sure_path_exists(variable_importance_csv, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    X = drop_cols_not_used_in_ML(df_data)
    y = df_data["interface"]
    forest = THOIPA_RF_classifier_with_settings(set_)
    forest.fit(X, y)
    importances_arr = forest.feature_importances_
    std_arr = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    indices_arr = np.argsort(importances_arr)[::-1]

    importances_text_list = X.columns.tolist()

    order_list = [importances_arr[indices_arr[f]] for f in range(X.shape[1])]

    sys.stdout.write("\nFeature ranking:")

    nested_dict = {}

    for f in range(X.shape[1]):
        if f < 10:
            sys.stdout.write("\n%d. feature %d (%f) %s" % (f + 1, indices_arr[f], importances_arr[indices_arr[f]], importances_text_list[indices_arr[f]]))
            sys.stdout.flush()
        single_feature_dict = {"feature_order" : indices_arr[f], "mean_decrease_gini" : importances_arr[indices_arr[f]], "feature" : importances_text_list[indices_arr[f]],  "std" : std_arr[f]}
        nested_dict[f + 1] = single_feature_dict

    df_imp = pd.DataFrame(nested_dict).T
    df_imp["orig_order"] = df_imp.index
    df_imp.set_index("feature", inplace=True)

    df_imp.to_csv(variable_importance_csv)


def fig_variable_importance(set_, logging):
    """Create figures showing ML feature importance.

    Fig1 : Barchart all features
    Fig2 : Barchart top features (currently set at 30)

    Parameters
    ----------
    set_ : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    plt.style.use('seaborn-whitegrid')
    plt.rcParams['errorbar.capsize'] = 1
    plt.rcParams.update({'font.size': 8})
    from korbinian.utils import create_colour_lists
    colour_dict = create_colour_lists()

    variable_importance_csv = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance.csv".format(set_["setname"][-2:], set_["setname"][-2:]))
    variable_importance_all_png = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance_all.png".format(set_["setname"][-2:], set_["setname"][-2:]))
    variable_importance_top_png = os.path.join(set_["set_results_folder"], "crossvalidation", "trainset{}_testset{}_variable_importance_top.png".format(set_["setname"][-2:], set_["setname"][-2:]))

    df_imp = pd.read_csv(variable_importance_csv, index_col = 0)

    create_var_imp_plot(df_imp, colour_dict, variable_importance_all_png, df_imp.shape[0])
    create_var_imp_plot(df_imp, colour_dict, variable_importance_top_png, 30)

    # UPDATE AS ABOVE
    # data = pd.read_csv('/scratch2/zeng/homotypic_data/data/RandomForest/PsEnCo/TrainData2',delimiter="\s")
    # del data["Residue_num"]
    # del data["Residue_name"]
    # #print(data.as_matrix(data.columns))
    # features=data.columns[0:28]
    # X=data[features]
    # y=data["Bind"]
    # forest = RandomForestClassifier(n_estimators=100)
    # forest.fit(X, y)
    # importances_arr = forest.feature_importances_
    # std_arr = np.std_arr([tree.feature_importances_ for tree in forest.estimators_],
    #              axis=0)
    # indices = np.argsort(importances_arr)[::-1]
    # print("Feature ranking:")
    #
    # for f in range(X.shape[1]):
    #     print("%d. feature %d (%f)" % (f + 1, indices[f], importances_arr[indices[f]]))
    # plt.figure()
    # plt.title("Feature importances_arr")
    # plt.bar(range(X.shape[1]), importances_arr[indices],
    #        color="r", yerr=std_arr[indices], align="center")
    # #plt.xticks(range(x.shape[1]), indices)
    # plt.xticks(range(X.shape[1]), indices)
    # plt.xticks(range(X.shape[1]), indices)
    # plt.xlim([-1, X.shape[1]])
    # plt.show()

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
    fig.savefig(variable_importance_png[:-4] + ".pdf")

