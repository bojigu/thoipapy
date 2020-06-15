import warnings
from pathlib import Path

import thoipapy.validation.feature_selection

warnings.simplefilter(action='ignore', category=FutureWarning)
import thoipapy
import pandas as pd
import os
import numpy as np
import joblib
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
import pickle
import warnings
from scipy.stats import linregress
#from eccpy.tools import normalise_0_1
from thoipapy.utils import normalise_0_1
warnings.filterwarnings("ignore")
#
# def Test_Etra_deprecated(s):
#     testsetname = "set03"
#     train_set_list = s["train_datasets"].split(",")
#     for train_set in train_set_list:
#         trainsetname = "set{:02d}".format(int(train_set))
#         traindata_set = os.path.join(s["Results_folder"], trainsetname, "{}_train_data.csv".format(trainsetname))
#         train_df =  pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)
#
#         BO_curve_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Train{}_Test{}.bocurve.csv".format(trainsetname, testsetname))
#
#         testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["sets_folder"])
#
#         testdataset_df = pd.read_excel(testset_path,sheet_name="proteins")
#         acc_list = testdataset_df.acc.tolist()
#         database = testdataset_df.database[0]
#         dfc = pd.DataFrame()
#         for i in testdataset_df.index:
#             acc = testdataset_df.loc[i, "acc"]
#             database = testdataset_df.loc[i, "database"]
#             testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
#             test_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)
#
#             #odf = save_THOIPA_pred_indiv_prot(acc, train_df, test_df)
#             odf = save_THOIPA_pred_indiv_prot(s, model_pkl, testdata_combined_file, THOIPA_pred_csv, combined_incl_THOIPA_csv, logging)
#             if dfc.empty:
#                 dfc = odf
#             else:
#                 dfc = pd.concat([dfc, odf], axis=1, join="outer")
#
#         dfc.to_csv(BO_curve_csv)


def run_testset_trainset_validation(s, logging):
    if not os.path.exists(os.path.join(s["thoipapy_data_folder"], "Results", "Bo_curve")):
        os.makedirs(os.path.join(s["thoipapy_data_folder"], "Results", "Bo_curve"))

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset
    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    validate_LIPS_for_testset(s, logging)
    validate_LIPS_for_testset(s, logging, LIPS_name="LIPS_surface_ranked", pred_col="LIPS_surface_ranked")

    validate_THOIPA_for_testset_trainset_combination(s, test_set_list, train_set_list, logging)


def validate_THOIPA_for_testset_trainset_combination(s, test_set_list, train_set_list, logging):
    """ Creates ROC and BO-curve for a particular testset-trainset combination.

    Parameters
    ----------
	s : dict
        Settings dictionary for figures.
    test_set_list : list
        List of test datasets in selection
        E.g. ["set03", "set31"]
    train_set_list : list
        List of training datasets in selection
        E.g. ["set02", "set04"]

    Saved Files
    -----------
    THOIPA_pred_csv : csv
        THOIPA result for this testset-trainset combination
        Columns = "residue_num", "residue_name", "THOIPA"
        Index = range index of residues
    combined_incl_THOIPA_csv : csv
        The combined file with all features. THOIPA prediction is added as a new column
    THOIPA_ROC_pkl : pickle
        Pickled output dictionary with ROC curves
        keys = accessions
        values = dictionary with fpr, tpr etc for each protein
        Could not be saved easily as a dataframe, because the number of residues is different for each protein

    """
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    for n, train_set in enumerate(train_set_list):
        trainsetname = "set{:02d}".format(int(train_set))
        model_pkl = os.path.join(s["Results_folder"], trainsetname, "{}_ML_model.lpkl".format(trainsetname))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            THOIPA_BO_curve_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "data", "Test{}_Train{}.THOIPA.best_overlap_data.csv".format(testsetname, trainsetname))
            THOIPA_ROC_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "data", "Test{}_Train{}.THOIPA.ROC_data.pkl".format(testsetname, trainsetname))
            BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname))
            BO_data_excel = os.path.join(BO_curve_folder, "data", "BO_curve_data.xlsx")
            BO_linechart_png = os.path.join(BO_curve_folder, "BO_linechart.png")
            BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")

            thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

            testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["sets_folder"])

            testdataset_df = pd.read_excel(testset_path)
            THOIPA_BO_data_df = pd.DataFrame()
            #LIPS_BO_data_df = pd.DataFrame()

            # save all outputs to a cross-validation dictionary, to be saved as a pickle file
            xv_dict_THOIPA = {}
            mean_tpr = 0.0
            mean_fpr = np.linspace(0, 1, 100)

            for i in testdataset_df.index:
                acc = testdataset_df.loc[i, "acc"]
                database = testdataset_df.loc[i, "database"]
                acc_db = acc + "-" + database
                testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                      "{}.surr20.gaps5.combined_features.csv".format(acc))
                THOIPA_pred_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/data/{database}/{acc}.THOIPA.train{trainsetname}.csv"
                combined_incl_THOIPA_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/data/{database}/{acc}.THOIPA_incl_combined.train{trainsetname}.csv"

                #THOIPA_pred_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database,
                #                                      "{}.THOIPA.train{}.csv".format(acc, trainsetname))
                #combined_incl_THOIPA_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database,
                #                                      "{}.THOIPA_incl_combined.train{}.csv".format(acc, trainsetname))

                thoipapy.utils.make_sure_path_exists(combined_incl_THOIPA_csv, isfile=True)

                combined_incl_THOIPA_df = save_THOIPA_pred_indiv_prot(s, model_pkl, testdata_combined_file, THOIPA_pred_csv, combined_incl_THOIPA_csv, logging)

                #######################################################################################################
                #                                                                                                     #
                #                           Processing BO curve data for each single protein                          #
                #                                                                                                     #
                #######################################################################################################

                combined_incl_THOIPA_df["LIPS_L*E"] = -1 * combined_incl_THOIPA_df["LIPS_L*E"]

                if database == "crystal" or database == "NMR":
                    # (it is closest distance and low value means high propencity of interfacial)
                    combined_incl_THOIPA_df["interface_score"] = -1 * combined_incl_THOIPA_df["interface_score"]

                THOIPA_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, combined_incl_THOIPA_df)

                if THOIPA_BO_data_df.empty:
                    THOIPA_BO_data_df = THOIPA_BO_single_prot_df
                else:
                    THOIPA_BO_data_df = pd.concat([THOIPA_BO_data_df, THOIPA_BO_single_prot_df], axis=1, join="outer")

                #######################################################################################################
                #                                                                                                     #
                #                     Processing ROC data for each single protein, saving to nested dict              #
                #                                                                                                     #
                #######################################################################################################

                df_for_roc = combined_incl_THOIPA_df.dropna(subset=["interface_score"])

                predictor_name = "THOIPA"

                fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name], drop_intermediate=False)
                roc_auc = auc(fpr, tpr)
                mean_tpr += interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0

                xv_dict_THOIPA[acc_db] = {"fpr" : fpr, "tpr" : tpr, "roc_auc" : roc_auc}

            #######################################################################################################
            #                                                                                                     #
            #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
            #                                                                                                     #
            #######################################################################################################

            THOIPA_BO_data_df.to_csv(THOIPA_BO_curve_data_csv)

            #THOIPA_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(THOIPA_BO_curve_data_csv, BO_curve_folder, names_excel_path)
            parse_BO_data_csv_to_excel(THOIPA_BO_curve_data_csv, BO_data_excel, logging)
            AUC_ser = pd.Series(xv_dict_THOIPA[acc_db]["roc_auc"])
            AUBOC10 = save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser)

            if "you_want_more_details" == "TRUE":
                other_figs_path = os.path.join(BO_curve_folder, "other_figs")
                save_extra_BO_figs(BO_data_excel, other_figs_path)

            #######################################################################################################
            #                                                                                                     #
            #                     Processing dictionary with ROC data, saving to pickle                           #
            #                                                                                                     #
            #######################################################################################################
            mean_tpr /= testdataset_df.shape[0]
            mean_tpr[-1] = 1.0

            mean_roc_auc = auc(mean_fpr, mean_tpr)

            ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_THOIPA}
            ROC_out_dict["true_positive_rate_mean"] = mean_tpr
            ROC_out_dict["false_positive_rate_mean"] = mean_fpr
            ROC_out_dict["mean_roc_auc"] = mean_roc_auc

            # save dict as pickle
            with open(THOIPA_ROC_pkl, "wb") as f:
                pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

            create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl)

            logging.info("Test{}_Train{} AUC({:.03f}), AUBOC10({:.2f}). ({})".format(testsetname, trainsetname, mean_roc_auc, AUBOC10, BO_barchart_png))


def validate_LIPS_for_testset(s, logging, LIPS_name = "LIPS_LE", pred_col="LIPS_L*E"):

    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset
    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    for test_set in test_set_list:
        testsetname = "set{:02d}".format(int(test_set))
        LIPS_BO_curve_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name), "Test{}.{}.best_overlap_data.csv".format(testsetname, LIPS_name))
        LIPS_ROC_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name), "data", "Test{}.{}.ROC_data.pkl".format(testsetname, LIPS_name))
        BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name))
        thoipapy.utils.make_sure_path_exists(LIPS_BO_curve_data_csv, isfile=True)
        BO_data_excel = os.path.join(BO_curve_folder, "data", "BO_curve_data.xlsx")
        BO_linechart_png = os.path.join(BO_curve_folder, "BO_linechart.png")
        BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")
        thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

        testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["sets_folder"])

        testdataset_df = pd.read_excel(testset_path)
        LIPS_BO_data_df = pd.DataFrame()

        # save all outputs to a cross-validation dictionary, to be saved as a pickle file
        xv_dict_LIPS = {}
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        for i in testdataset_df.index:
            acc = testdataset_df.loc[i, "acc"]
            database = testdataset_df.loc[i, "database"]
            acc_db = acc + "-" + database

            testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                  "{}.surr20.gaps5.combined_features.csv".format(acc))

            combined_df = pd.read_csv(testdata_combined_file, index_col=0)

            #######################################################################################################
            #                                                                                                     #
            #                           Processing BO curve data for each single protein                          #
            #                                                                                                     #
            #######################################################################################################
            # SAVE LIPS PREDICTION DATA
            # this is somewhat inefficient, as it is conducted for every test dataset
            LIPS_pred_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database, "{}.LIPS_pred.csv".format(acc, testsetname))
            LIPS_pred_df = combined_df[["residue_name", "residue_num", "LIPS_polarity", "LIPS_entropy", "LIPS_L*E", "LIPS_surface", "LIPS_surface_ranked"]]
            thoipapy.utils.make_sure_path_exists(LIPS_pred_csv, isfile=True)
            LIPS_pred_df.to_csv(LIPS_pred_csv)

            if pred_col == "LIPS_L*E":
                combined_df[pred_col] = -1 * combined_df[pred_col]

            if database == "crystal" or database == "NMR":
                # (it is closest distance and low value means high propencity of interfacial)
                combined_df["interface_score"] = -1 * combined_df["interface_score"]

            LIPS_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, combined_df, experiment_col="interface_score", pred_col=pred_col)

            if LIPS_BO_single_prot_df.empty:
                LIPS_BO_data_df = LIPS_BO_single_prot_df
            else:
                LIPS_BO_data_df = pd.concat([LIPS_BO_data_df, LIPS_BO_single_prot_df], axis=1, join="outer")

            #######################################################################################################
            #                                                                                                     #
            #                     Processing ROC data for each single protein, saving to nested dict              #
            #                                                                                                     #
            #######################################################################################################

            df_for_roc = combined_df.dropna(subset=["interface_score"])

            fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[pred_col], drop_intermediate=False)
            roc_auc = auc(fpr, tpr)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0

            xv_dict_LIPS[acc_db] = {"fpr" : fpr, "tpr" : tpr, "roc_auc" : roc_auc}

        #######################################################################################################
        #                                                                                                     #
        #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
        #                                                                                                     #
        #######################################################################################################

        LIPS_BO_data_df.to_csv(LIPS_BO_curve_data_csv)
        names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

        #LIPS_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(LIPS_BO_curve_data_csv, BO_curve_folder, names_excel_path)

        #parse_BO_data_csv_to_excel(LIPS_BO_curve_data_csv, BO_curve_folder, names_excel_path)
        parse_BO_data_csv_to_excel(LIPS_BO_curve_data_csv, BO_data_excel, logging)
        AUC_ser = pd.Series(xv_dict_LIPS[acc_db]["roc_auc"])
        AUBOC10 = save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser)

        if "you_want_more_details" == "TRUE":
            other_figs_path = os.path.join(BO_curve_folder, "other_figs")
            save_extra_BO_figs(BO_data_excel, other_figs_path)

        #######################################################################################################
        #                                                                                                     #
        #                     Processing dictionary with ROC data, saving to pickle                           #
        #                                                                                                     #
        #######################################################################################################
        mean_tpr /= testdataset_df.shape[0]
        mean_tpr[-1] = 1.0

        mean_roc_auc = auc(mean_fpr, mean_tpr)

        ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_LIPS}
        ROC_out_dict["true_positive_rate_mean"] = mean_tpr
        ROC_out_dict["false_positive_rate_mean"] = mean_fpr
        ROC_out_dict["mean_roc_auc"] = mean_roc_auc

        # save dict as pickle
        with open(LIPS_ROC_pkl, "wb") as f:
            pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

        create_ROC_fig_for_testset_trainset_combination(LIPS_ROC_pkl)

        logging.info("Test{}.{} AUC({:.03f}), AUBOC10({:.2f}). ({})".format(testsetname, LIPS_name, mean_roc_auc, AUBOC10, BO_barchart_png))

def create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl):

    #plt.rcParams.update({'font.size': 6})
    ROC_pkl_basename = os.path.basename(THOIPA_ROC_pkl)[:-4]
    ROC_pkl_dir = os.path.dirname(THOIPA_ROC_pkl)

    ROC_png = os.path.join(ROC_pkl_dir, "{}.ROC.png".format(ROC_pkl_basename))
    thoipapy.utils.make_sure_path_exists(ROC_png, isfile=True)

    # open pickle file
    with open(THOIPA_ROC_pkl, "rb") as f:
        ROC_out_dict = pickle.load(f)

    xv_dict_THOIPA = ROC_out_dict["xv_dict_THOIPA"]

    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for acc_db in xv_dict_THOIPA:
        roc_auc = xv_dict_THOIPA[acc_db]["roc_auc"]
        ax.plot(xv_dict_THOIPA[acc_db]["fpr"], xv_dict_THOIPA[acc_db]["tpr"], lw=1, label='{} ({:0.2f})'.format(acc_db, roc_auc), alpha=0.8)

    #mean_roc_auc = auc(df_xv["false_positive_rate"], df_xv["true_positive_rate"])
    mean_roc_auc = ROC_out_dict["mean_roc_auc"]

    ax.plot(ROC_out_dict["false_positive_rate_mean"], ROC_out_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_roc_auc, lw=1.5)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_png, dpi=240)
    fig.savefig(thoipapy.utils.pdf_subpath(ROC_png))
    #fig.savefig(ROC_png[:-4] + ".pdf")


def save_THOIPA_pred_indiv_prot(s, model_pkl, testdata_combined_file, THOIPA_pred_csv, test_combined_incl_pred, logging):

    combined_incl_THOIPA_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)
    #combined_incl_THOIPA_df = combined_incl_THOIPA_df.dropna()

    #drop_cols_not_used_in_ML
    #X=train_df.drop(train_features_del,axis=1)
    # X = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, train_df, s["excel_file_with_settings"])
    # y = train_df["interface"]
    # clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    # fit = clf.fit(X,y)

    fit = joblib.load(model_pkl)

    # if database == "ETRA":
    #     dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
    #     ddf = pd.read_excel(dirupt_path, index_col=0)
    #     interface_score = ddf.Disruption
    #     interface_score = interface_score      #(it is closest experimental disruption and high value means high propencity of interfacial)

    #Lips_score = test_df.LIPS_polarity * test_df.LIPS_entropy

    #tX=test_df.drop(test_features_del,axis=1)
    test_X = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, combined_incl_THOIPA_df, s["excel_file_with_settings"])

    try:
        prob_arr = fit.predict_proba(test_X)[:, 1]
    except ValueError:
        train_data_used_for_model_csv = model_pkl[:-11] + "train_data_used_for_model.csv"
        df_train_data_used_for_model_csv = pd.read_csv(train_data_used_for_model_csv)
        train_cols = set(df_train_data_used_for_model_csv.columns)
        test_cols = set(test_X.columns)
        cols_not_in_train = test_cols - train_cols
        cols_not_in_test = train_cols - test_cols
        logging.warning("cols_not_in_train : {}\ncols_not_in_test : {}".format(cols_not_in_train, cols_not_in_test))
        prob_arr = fit.predict_proba(test_X)[:, 1]

    # if hasattr(clf,'predict_proba'):
    #     prob_arr = fit.predict_proba(tX)[:,1]
    # else:
    #     prob_arr = fit.decision_function(tX)
    #     prob_arr = (prob_arr - prob_arr.min())/(prob_arr.max() - prob_arr.min())

    combined_incl_THOIPA_df["THOIPA"] = prob_arr
    # save full combined file with THOIPA prediction
    combined_incl_THOIPA_df.to_csv(test_combined_incl_pred)
    # save only THOIPA prediction
    THOIPA_pred_df = combined_incl_THOIPA_df[["residue_num", "residue_name", "THOIPA"]]
    THOIPA_pred_df.to_csv(THOIPA_pred_csv)
    return combined_incl_THOIPA_df

#
# def save_THOIPA_pred_indiv_prot(acc, train_df, test_df, database):
#     #drop_cols_not_used_in_ML
#     #X=train_df.drop(train_features_del,axis=1)
#     X = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, train_df, s["excel_file_with_settings"])
#     y = train_df["interface"]
#     clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
#     fit = clf.fit(X,y)
#
#     interface_score = test_df.interface_score
#
#     if database == "crystal" or database == "NMR":
#         interface_score = test_df.interface_score
#         interface_score = -1 * interface_score  #(it is closest distance and low value means high propencity of interfacial)
#     # if database == "ETRA":
#     #     dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
#     #     ddf = pd.read_excel(dirupt_path, index_col=0)
#     #     interface_score = ddf.Disruption
#     #     interface_score = interface_score      #(it is closest experimental disruption and high value means high propencity of interfacial)
#     test_df.index = test_df.index.astype(int) + 1
#     Lips_score = test_df.LIPS_polarity * test_df.LIPS_entropy
#
#     #tX=test_df.drop(test_features_del,axis=1)
#     tX = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, test_df, s["excel_file_with_settings"])
#
#     if hasattr(clf,'predict_proba'):
#         prob_pos = fit.predict_proba(tX)[:,1]
#     else:
#         prob_pos = fit.decision_function(tX)
#         prob_pos = (prob_pos - prob_pos.min())/(prob_pos.max() - prob_pos.min())
#
#     if database == "crystal" or database == "NMR":
#         interface_score = -1 * interface_score  # (it is closest distance and low value means high propencity of interfacial)
#     elif database == "ETRA":
#         pass  # (it is closest experimental disruption and high value means high propencity of interfacial)
#     else:
#         raise ValueError()
#
#     THOIPA_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc, prob_pos, interface_score, database)
#
#     return THOIPA_BO_single_prot_df



# def create_extra_BO_curve_figs(bo_data_csv, out_folder, names_excel_path):
#     """Analyse the Bo-curve underlying data.
#
#     Parses data into more easily manipulated dataframes.
#     Creates figures if desired.
#     All figures are saved in a subfolder with the same name as the original file.
#
#     Parameters
#     ----------
#     bo_data_csv : str
#         Path to the csv created by bo, with the underlying data.
#         Index : Top1, Top2, Ono, Pno, Rno etc
#     names_excel_path : str
#         Path to the excel file with the protein short names and reference.
#
#     Usage
#     -----
#     import datoxr
#     bo_data_csv = r"D:\drive\TMD_homodimer\figs\SuppDataX02-best_overlap_data\SuppDataX02.csv"
#     names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
#     datoxr.figs.bo_curve_analysis.analyse_bo_curve_underlying_data(bo_data_csv, names_excel_path)
#     """
#
#     # change to empty list if you don't want to create figures. Or [7] if you only want to process figure 7, for example.
#     list_figs_to_create = range(1, 11)
#     #list_figs_to_create = [7]
#
#     other_figs_path = os.path.join(out_folder, "other_figs")
#     if not os.path.exists(other_figs_path):
#         os.makedirs(other_figs_path)
#
#     BO_data_excel = os.path.join(out_folder, "BO_curve_data.xlsx")
#     # linechart_mean_obs_and_rand = os.path.join(other_figs_path, "1_linechart_mean_obs_and_rand.png")
#     # linechart_obs_indiv = os.path.join(other_figs_path, "2_linechart_obs_indiv.png")
#     # linechart_p_indiv = os.path.join(other_figs_path, "3_linechart_p_indiv.png")
#     # linechart_o_minus_r = os.path.join(other_figs_path, "4_linechart_o_minus_r.png")
#     # linechart_o_over_r = os.path.join(other_figs_path, "5_linechart_o_over_r.png")
#     BO_linechart_png = os.path.join(out_folder, "BO_curve_single_dataset.png")
#     BO_barchart_png = os.path.join(out_folder, "performance_individual_proteins.png")
#
#     analyse_bo_curve_underlying_data(bo_data_csv, BO_data_excel)

    # dfb = pd.read_csv(bo_data_csv, index_col=0)
    #
    # """ORIGINAL BO DATA CSV LOOKS LIKE THIS
    # Top1 = sample size 1
    # Ono = overlap in data
    # Rno = random overlap based on that sequence length and sample size
    # Pono = p-value for finding that overlap
    #
    #      Unnamed: 1  O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286  Ratio (Average(Ono)/Average(Rno))
    # NaN         Ono    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    1.00    0.00    0.00    1.00    1.00    1.00    0.00    0.00    0.00                                NaN
    # Top1        Rno    0.05    0.04    0.05    0.04    0.05    0.06    0.04    0.04    0.05    0.04    0.04    0.05    0.06    0.05    0.06    0.04    0.05                               5.01
    # NaN        Pono    0.95    0.96    0.95    0.96    0.95    0.94    0.96    0.96    0.05    0.96    0.96    0.05    0.06    0.05    0.94    0.96    0.95                                NaN
    # NaN         Ono    1.00    1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    1.00    0.00    1.00    2.00    2.00    0.00    0.00    1.00                                NaN
    # Top2        Rno    0.19    0.17    0.21    0.15    0.19    0.24    0.17    0.17    0.19    0.17    0.16    0.19    0.25    0.20    0.24    0.17    0.20                               3.68
    # """
    #
    # # create an index based on sample size [1 1 1 2 2 2  etc..
    # ind = []
    # for i in range(1, int((len(dfb) / 3)) + 1):
    #     ind += list(np.array([1, 1, 1]) * i)
    # dfb.index = ind
    # dfb.index.name = "sample size"
    #
    # """NOW INDICES ARE BASED ON SAMPLE SIZE
    #
    #             Unnamed: 1  O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286  Ratio (Average(Ono)/Average(Rno))
    # sample size
    # 1                  Ono    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    1.00    0.00    0.00    1.00    1.00    1.00    0.00    0.00    0.00                                NaN
    # 1                  Rno    0.05    0.04    0.05    0.04    0.05    0.06    0.04    0.04    0.05    0.04    0.04    0.05    0.06    0.05    0.06    0.04    0.05                               5.01
    # 1                 Pono    0.95    0.96    0.95    0.96    0.95    0.94    0.96    0.96    0.05    0.96    0.96    0.05    0.06    0.05    0.94    0.96    0.95                                NaN
    # 2                  Ono    1.00    1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    1.00    0.00    1.00    2.00    2.00    0.00    0.00    1.00                                NaN
    # 2                  Rno    0.19    0.17    0.21    0.15    0.19    0.24    0.17    0.17    0.19    0.17    0.16    0.19    0.25    0.20    0.24    0.17    0.20
    # """
    #
    # # split into separate dataframes
    # # dataframe of observed overlaps
    # dfobs = dfb.iloc[::3, 1:-1].astype(int)
    # # dataframe of random calculated overlaps
    # dfrand = dfb.iloc[1::3, 1:-1]
    # # dataframe of p-values
    # dfp = dfb.iloc[2::3, 1:-1]
    #
    # """FOR EXAMPLE df_obs NOW LOOKS LIKE THIS, WITH A ROW FOR EACH SAMPLE SIZE:
    #
    #              O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286
    # sample size
    # 1                 0       0       0       0       0       0       0       0       1       0       0       1       1       1       0       0       0
    # 2                 1       1       0       1       0       0       1       0       1       1       0       1       2       2       0       0       1
    # 3                 1       1       1       2       1       0       1       0       2       2       0       2       2       2       0       0       2
    # 4                 2       2       1       2       1       0       2       1       3       2       1       3       2       2       0       1       2
    # 5                 3       3       2       3       2       1       2       2       4       4       2       4       3       2       1       2       3"""
    #
    # df_o_minus_r = dfobs - dfrand
    # df_o_over_r = dfobs / dfrand
    #
    # """df_o_minus_r is negative where the result is lower than random
    #
    #              O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286
    # sample size
    # 1             -0.05   -0.04   -0.05   -0.04   -0.05   -0.06   -0.04   -0.04    0.95   -0.04   -0.04    0.95    0.94    0.95   -0.06   -0.04   -0.05
    # 2              0.81    0.83   -0.21    0.85   -0.19   -0.24    0.83   -0.17    0.81    0.83   -0.16    0.81    1.75    1.80   -0.24   -0.17    0.80
    # 3              0.57    0.61    0.53    1.65    0.57   -0.53    0.61   -0.37    1.57    1.63   -0.36    1.57    1.44    1.55   -0.53   -0.37    1.55
    # 4              1.24    1.30    0.16    1.38    0.24   -0.94    1.30    0.33    2.24    1.33    0.36    2.24    1.00    1.20   -0.94    0.33    1.20
    # 5              1.81    1.91    0.68    2.04    0.81   -0.47    0.91    0.96    2.81    2.96    1.00    2.81    1.44    0.75   -0.47    0.96    1.75
    #
    #
    # df_o_over_r is where the result is lower than random.
    # This is not quite a fair comparison, as the zeros are caused by 0 overlap / signigicant random overlap
    #
    #                O75460    P02724    P05106    P06583    P08514    P0A6S5    P23470    P35590     Q08345    Q12983  Q16827     Q16832     Q6ZRP7     Q7L4S7    Q8NI60    Q92729    Q9Y286
    # sample size
    # 1            0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  20.000000  0.000000  0.0000  20.000000  16.666667  20.000000  0.000000  0.000000  0.000000
    # 2            5.263158  5.882353  0.000000  6.666667  0.000000  0.000000  5.882353  0.000000   5.263158  5.882353  0.0000   5.263158   8.000000  10.000000  0.000000  0.000000  5.000000
    # 3            2.325581  2.564103  2.127660  5.714286  2.325581  0.000000  2.564103  0.000000   4.651163  5.405405  0.0000   4.651163   3.571429   4.444444  0.000000  0.000000  4.444444
    # 4            2.631579  2.857143  1.190476  3.225806  1.315789  0.000000  2.857143  1.492537   3.947368  2.985075  1.5625   3.947368   2.000000   2.500000  0.000000  1.492537  2.500000
    # 5            2.521008  2.752294  1.515152  3.125000  1.680672  0.680272  1.834862  1.923077   3.361345  3.846154  2.0000   3.361345   1.923077   1.600000  0.680272  1.923077  2.400000
    # """
    #
    # #################################################################
    # #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    # #################################################################
    #
    # with pd.ExcelWriter(BO_data_excel) as writer:
    #     dfobs.to_excel(writer, sheet_name="dfobs")
    #     dfrand.to_excel(writer, sheet_name="dfrand")
    #     dfp.to_excel(writer, sheet_name="dfp")
    #     df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
    #     df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")


def parse_BO_data_csv_to_excel(bo_data_csv, BO_data_excel, logging, predictor_name=""):
    """

    Run using s["create_AUC_AUBOC_separate_database"]

    Parameters
    ----------
    bo_data_csv
    BO_data_excel
    logging
    predictor_name

    Returns
    -------

    """
    dfb = pd.read_csv(bo_data_csv, index_col=0)

    """ORIGINAL BO DATA CSV LOOKS LIKE THIS
    Top1 = sample size 1
    observed_overlap = overlap in data
    random_overlap = random overlap based on that sequence length and sample size
    p_value_from_obs_overlap = p-value for finding that overlap

                               parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    sample_size                                                                                                                                                                               ...                                                                                                                                       
    Top1                 observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1                   random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1         p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2                 observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2                   random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20
    """

    dfb.index = dfb.index + "_" + dfb.parameters

    """NOW INDICES ARE UNIQUE

                                             parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    parameters                                                                                                                                                                                                  ...                                                                                                                                       
    Top1_observed_overlap                  observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1_random_overlap                      random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1_p_value_from_obs_overlap  p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2_observed_overlap                  observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2_random_overlap                      random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20

    """
    cols_to_drop = set(["parameters", "Ratio (Average(Ono)/Average(Rno))"]).intersection(set(dfb.columns))
    dfb.drop(cols_to_drop, axis=1, inplace=True)

    # split into separate dataframes. Relabel index to match sample size.
    # dataframe of observed overlaps
    dfobs = dfb[dfb.index.str.contains("observed_overlap")].astype(int)
    dfobs.index = range(1, dfobs.shape[0] + 1)
    # dataframe of random calculated overlaps
    dfrand = dfb[dfb.index.str.contains("random_overlap")]
    dfrand.index = range(1, dfrand.shape[0] + 1)
    # dataframe of p-values
    dfp = dfb[dfb.index.str.contains("p_value_from_obs_overlap")]
    dfp.index = range(1, dfp.shape[0] + 1)

    """FOR EXAMPLE df_obs NOW LOOKS LIKE THIS, WITH A ROW FOR EACH SAMPLE SIZE:

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1               0               0               0               1               0               0               0               0               0               0     ...                 1            1            0            0            0            0            0            0            1            0
    2               0               0               0               2               0               1               0               0               1               0     ...                 1            1            0            1            0            2            2            1            1            2
    3               1               0               0               3               1               1               0               0               1               0     ...                 3            3            0            2            0            2            2            1            1            2
    4               1               0               0               3               1               1               0               0               3               0     ...                 3            4            1            4            0            2            3            2            2            2
    5               1               0               0               4               2               2               2               0               4               0     ...                 3            4            1            4            1            2            3            3            3            2
    """

    # convert the observed overlap (e.g. 4 residues from sample size 10) to a fraction (0.4 of 1)
    dfobs_frac = dfobs.divide(np.array(dfobs.index), axis=0)
    dfrand_frac = dfrand.divide(np.array(dfrand.index), axis=0)

    """If dfobs looks like this:
        ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1               0           0            0
    2               0           0            0
    3               0           0            1
    4               1           0            1
    5               1           1            2
    
    
    Then the dfobs_frac is simply the same data as a fraction (divided by the index)
    
    ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1            0.00         0.0     0.000000
    2            0.00         0.0     0.000000
    3            0.00         0.0     0.333333
    4            0.25         0.0     0.250000
    5            0.20         0.2     0.400000
    """

    # the observed minus random is now calculated for the fraction correct, rather than the original numbers
    df_o_minus_r = dfobs_frac - dfrand_frac
    #df_o_over_r = dfobs_frac / dfrand_frac

    """df_o_minus_r is negative where the result is lower than random

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1       -0.035714           -0.05           -0.05            0.95           -0.05       -0.045455       -0.043478       -0.034483       -0.043478       -0.034483     ...              0.95     0.958333        -0.04    -0.047619    -0.047619        -0.05    -0.047619    -0.041667     0.952381        -0.05
    2       -0.142857           -0.20           -0.20            1.80           -0.20        0.818182       -0.173913       -0.137931        0.826087       -0.137931     ...              0.80     0.833333        -0.16     0.809524    -0.190476         1.80     1.809524     0.833333     0.809524         1.80
    3        0.678571           -0.45           -0.45            2.55            0.55        0.590909       -0.391304       -0.310345        0.608696       -0.310345     ...              2.55     2.625000        -0.36     1.571429    -0.428571         1.55     1.571429     0.625000     0.571429         1.55
    4        0.428571           -0.80           -0.80            2.20            0.20        0.272727       -0.695652       -0.551724        2.304348       -0.551724     ...              2.20     3.333333         0.36     3.238095    -0.761905         1.20     2.238095     1.333333     1.238095         1.20
    5        0.107143           -1.25           -1.25            2.75            0.75        0.863636        0.913043       -0.862069        2.913043       -0.862069     ...              1.75     2.958333         0.00     2.809524    -0.190476         0.75     1.809524     1.958333     1.809524         0.75


    df_o_over_r is where the result is lower than random.
    This is probably not quite a fair comparison, as the zeros are caused by 0 overlap / signigicant random overlap

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1        0.000000             0.0             0.0       20.000000        0.000000        0.000000            0.00             0.0        0.000000             0.0     ...         20.000000        24.00       0.0000     0.000000         0.00     0.000000     0.000000     0.000000    21.000000     0.000000
    2        0.000000             0.0             0.0       10.000000        0.000000        5.500000            0.00             0.0        5.750000             0.0     ...          5.000000         6.00       0.0000     5.250000         0.00    10.000000    10.500000     6.000000     5.250000    10.000000
    3        3.111111             0.0             0.0        6.666667        2.222222        2.444444            0.00             0.0        2.555556             0.0     ...          6.666667         8.00       0.0000     4.666667         0.00     4.444444     4.666667     2.666667     2.333333     4.444444
    4        1.750000             0.0             0.0        3.750000        1.250000        1.375000            0.00             0.0        4.312500             0.0     ...          3.750000         6.00       1.5625     5.250000         0.00     2.500000     3.937500     3.000000     2.625000     2.500000
    5        1.120000             0.0             0.0        3.200000        1.600000        1.760000            1.84             0.0        3.680000             0.0     ...          2.400000         3.84       1.0000     3.360000         0.84     1.600000     2.520000     2.880000     2.520000     1.600000
    
    NOW df_o_minus_r EXPRESSED AS A FRACTION
    It's negative where the random value got a higher score.
    
        ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1       -0.034483   -0.043478     0.933333
    2        0.431034   -0.086957     0.366667
    3        0.229885   -0.130435     0.466667
    4        0.362069   -0.173913     0.483333
    5        0.227586    0.182609     0.466667
    
    """
    #################################################################
    #          CALCULATE MEAN VALUE FOR ALL TMDS FOR BO-CURVE       #
    #################################################################
    mean_o_minus_r_df = df_o_minus_r.mean(axis=1).to_frame(name="mean_o_minus_r")

    #######################################################################################################
    #                                                                                                     #
    #                CALCULATE AREA UNDER THE BO CURVE FOR SAMPLE SIZES 1-10 (AUBOC10)                    #
    #                                                                                                     #
    #######################################################################################################

    AUBOC10_df = pd.DataFrame()
    for acc_db in df_o_minus_r.columns:
        o_minus_r_ser = df_o_minus_r[acc_db]
        AUBOC10_df.loc[acc_db, "AUBOC10"] = np.trapz(o_minus_r_ser, o_minus_r_ser.index)

        AUBOC10_df.loc[acc_db, "mean_all_sample_sizes"] = o_minus_r_ser.mean()

    mean_AUBOC10 = AUBOC10_df["AUBOC10"].mean()
    logging.info("---{: >24} mean_AUBOC10({:.2f}) n={} ---".format(predictor_name, mean_AUBOC10, AUBOC10_df.shape[0]))
    #################################################################
    #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    #################################################################

    with pd.ExcelWriter(BO_data_excel) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        dfobs_frac.to_excel(writer, sheet_name="dfobs_frac")
        dfrand_frac.to_excel(writer, sheet_name="dfrand_frac")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        #df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")
        mean_o_minus_r_df.to_excel(writer, sheet_name="mean_o_minus_r")
        AUBOC10_df.to_excel(writer, sheet_name="AUBOC10")


def parse_BO_data_csv_to_excel_DEPRECATED_NONFRACTION_VERSION(bo_data_csv, BO_data_excel, logging, predictor_name=""):
    """
    """
    dfb = pd.read_csv(bo_data_csv, index_col=0)

    """ORIGINAL BO DATA CSV LOOKS LIKE THIS
    Top1 = sample size 1
    observed_overlap = overlap in data
    random_overlap = random overlap based on that sequence length and sample size
    p_value_from_obs_overlap = p-value for finding that overlap

                               parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    sample_size                                                                                                                                                                               ...                                                                                                                                       
    Top1                 observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1                   random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1         p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2                 observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2                   random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20
    """

    dfb.index = dfb.index + "_" + dfb.parameters

    """NOW INDICES ARE UNIQUE

                                             parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    parameters                                                                                                                                                                                                  ...                                                                                                                                       
    Top1_observed_overlap                  observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1_random_overlap                      random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1_p_value_from_obs_overlap  p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2_observed_overlap                  observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2_random_overlap                      random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20

    """
    cols_to_drop = set(["parameters", "Ratio (Average(Ono)/Average(Rno))"]).intersection(set(dfb.columns))
    dfb.drop(cols_to_drop, axis=1, inplace=True)

    # split into separate dataframes. Relabel index to match sample size.
    # dataframe of observed overlaps
    dfobs = dfb[dfb.index.str.contains("observed_overlap")].astype(int)
    dfobs.index = range(1, dfobs.shape[0] + 1)
    # dataframe of random calculated overlaps
    dfrand = dfb[dfb.index.str.contains("random_overlap")]
    dfrand.index = range(1, dfrand.shape[0] + 1)
    # dataframe of p-values
    dfp = dfb[dfb.index.str.contains("p_value_from_obs_overlap")]
    dfp.index = range(1, dfp.shape[0] + 1)

    """FOR EXAMPLE df_obs NOW LOOKS LIKE THIS, WITH A ROW FOR EACH SAMPLE SIZE:

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1               0               0               0               1               0               0               0               0               0               0     ...                 1            1            0            0            0            0            0            0            1            0
    2               0               0               0               2               0               1               0               0               1               0     ...                 1            1            0            1            0            2            2            1            1            2
    3               1               0               0               3               1               1               0               0               1               0     ...                 3            3            0            2            0            2            2            1            1            2
    4               1               0               0               3               1               1               0               0               3               0     ...                 3            4            1            4            0            2            3            2            2            2
    5               1               0               0               4               2               2               2               0               4               0     ...                 3            4            1            4            1            2            3            3            3            2
    """

    df_o_minus_r = dfobs - dfrand
    df_o_over_r = dfobs / dfrand

    """df_o_minus_r is negative where the result is lower than random

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1       -0.035714           -0.05           -0.05            0.95           -0.05       -0.045455       -0.043478       -0.034483       -0.043478       -0.034483     ...              0.95     0.958333        -0.04    -0.047619    -0.047619        -0.05    -0.047619    -0.041667     0.952381        -0.05
    2       -0.142857           -0.20           -0.20            1.80           -0.20        0.818182       -0.173913       -0.137931        0.826087       -0.137931     ...              0.80     0.833333        -0.16     0.809524    -0.190476         1.80     1.809524     0.833333     0.809524         1.80
    3        0.678571           -0.45           -0.45            2.55            0.55        0.590909       -0.391304       -0.310345        0.608696       -0.310345     ...              2.55     2.625000        -0.36     1.571429    -0.428571         1.55     1.571429     0.625000     0.571429         1.55
    4        0.428571           -0.80           -0.80            2.20            0.20        0.272727       -0.695652       -0.551724        2.304348       -0.551724     ...              2.20     3.333333         0.36     3.238095    -0.761905         1.20     2.238095     1.333333     1.238095         1.20
    5        0.107143           -1.25           -1.25            2.75            0.75        0.863636        0.913043       -0.862069        2.913043       -0.862069     ...              1.75     2.958333         0.00     2.809524    -0.190476         0.75     1.809524     1.958333     1.809524         0.75


    df_o_over_r is where the result is lower than random.
    This is probably not quite a fair comparison, as the zeros are caused by 0 overlap / signigicant random overlap

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1        0.000000             0.0             0.0       20.000000        0.000000        0.000000            0.00             0.0        0.000000             0.0     ...         20.000000        24.00       0.0000     0.000000         0.00     0.000000     0.000000     0.000000    21.000000     0.000000
    2        0.000000             0.0             0.0       10.000000        0.000000        5.500000            0.00             0.0        5.750000             0.0     ...          5.000000         6.00       0.0000     5.250000         0.00    10.000000    10.500000     6.000000     5.250000    10.000000
    3        3.111111             0.0             0.0        6.666667        2.222222        2.444444            0.00             0.0        2.555556             0.0     ...          6.666667         8.00       0.0000     4.666667         0.00     4.444444     4.666667     2.666667     2.333333     4.444444
    4        1.750000             0.0             0.0        3.750000        1.250000        1.375000            0.00             0.0        4.312500             0.0     ...          3.750000         6.00       1.5625     5.250000         0.00     2.500000     3.937500     3.000000     2.625000     2.500000
    5        1.120000             0.0             0.0        3.200000        1.600000        1.760000            1.84             0.0        3.680000             0.0     ...          2.400000         3.84       1.0000     3.360000         0.84     1.600000     2.520000     2.880000     2.520000     1.600000
    """
    #######################################################################################################
    #                                                                                                     #
    #                CALCULATE AREA UNDER THE BO CURVE FOR SAMPLE SIZES 1-10 (AUBOC10)                    #
    #                                                                                                     #
    #######################################################################################################

    AUBOC10_df = pd.DataFrame()
    for acc_db in df_o_minus_r.columns:
        o_minus_r_ser = df_o_minus_r[acc_db]
        AUBOC10_df.loc[acc_db, "AUBOC10"] = np.trapz(o_minus_r_ser, o_minus_r_ser.index)

    mean_AUBOC10 = AUBOC10_df["AUBOC10"].mean()
    logging.info("---{: >24} mean_AUBOC10({:.2f}) n={} ---".format(predictor_name, mean_AUBOC10, AUBOC10_df.shape[0]))
    #################################################################
    #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    #################################################################

    with pd.ExcelWriter(BO_data_excel) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")
        AUBOC10_df.to_excel(writer, sheet_name="AUBOC10")



def save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser, plot_o_over_r=False):

    df_o_minus_r = pd.read_excel(BO_data_excel, sheet_name="df_o_minus_r", index_col=0)

    #######################################################################################################
    #                                                                                                     #
    #               Create a dataframe with AUBOC10 and AUC for individual protein (df_valid_indiv)       #
    #                                                                                                     #
    #######################################################################################################
    # load AUBOC10 values as a series
    AUBOC10_ser = pd.read_excel(BO_data_excel, sheet_name="AUBOC10", index_col=0)["AUBOC10"].copy()
    # select sample sizes 5 and 10
    df_valid_indiv = df_o_minus_r.loc[[5, 10], :].T.copy()
    df_valid_indiv["AUBOC10"] = AUBOC10_ser
    df_valid_indiv["ROC AUC"] = AUC_ser
    df_valid_indiv.sort_values("AUBOC10", axis=0, ascending=False, inplace=True)

    """ df_valid_indiv should now have the results from BO curve and ROC for each protein
    
                      AUBOC10  sample size 5  sample size 10   ROC AUC
    3ij4_A-crystal  17.456522       1.913043        1.652174  0.714286
    4wit_A-crystal  16.620000       2.000000        2.000000  0.622807
    Q08345-ETRA     16.571429       2.809524        2.238095  0.842593
    P04626-ETRA     16.456522       1.913043        1.652174  0.916667
    P25189-ETRA     14.634615       2.038462        2.153846  0.812500
    """

    #######################################################################################################
    #                                                                                                     #
    #                                plot correlation between AUBOC10 and ROC                             #
    #                                                                                                     #
    #######################################################################################################
    # BO_barchart_png
    plt.close("all")
    #plt.rcParams.update({'font.size': 8})
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    #df_valid_indiv_scatter = df_valid_indiv[["AUBOC10", "ROC AUC"]]
    df_valid_indiv.plot(kind="scatter", ax=ax, x="AUBOC10", y="ROC AUC", alpha=0.7)

    # calculate linear regression for fitted line
    slope, intercept, r_value, p_value, std_err = linregress(df_valid_indiv["AUBOC10"], df_valid_indiv["ROC AUC"])
    #fit_fn = np.poly1d(linear_regression)

    #slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    x_first_last_dp = np.array([df_valid_indiv["AUBOC10"].min(), df_valid_indiv["AUBOC10"].max()])
    y_fitted = x_first_last_dp * slope + intercept
    ax.plot(x_first_last_dp, y_fitted, label="$R^2$ : {:.2f}".format(r_value**2))

    ax.set_xlabel("AUBOC10")
    ax.set_ylabel("ROC AUC")
    ax.legend()
    fig.tight_layout()
    ax.grid(False)
    #BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")

    fig.savefig(BO_barchart_png[:-12] + "scatter.png", dpi=240)

    # simply normalise all between 0 and 1
    for col in df_valid_indiv.columns:
        df_valid_indiv[col] = normalise_0_1(df_valid_indiv[col])[0] + 0.01

    BO_curve_folder = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / "crossvalidation"
    crossvalidation_data_dir = BO_curve_folder / "data"
    if not crossvalidation_data_dir.is_dir():
        crossvalidation_data_dir.mkdir(parents=True)
    BO_data_excel = os.path.join(BO_curve_folder, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))

    df_valid_indiv = df_valid_indiv.reindex(columns=["AUBOC10", 5, 10, "ROC AUC"])
    df_valid_indiv.columns = ["AUBOC10", "sample size 5", "sample size 10", "ROC AUC"]

    df_valid_indiv.to_csv(BO_data_excel[:-5] + "_valid_indiv.csv")

    """ df_valid_indiv is now normalised within each column, and sorted by AUBOC10
                          AUBOC10  sample size 5  sample size 10   ROC AUC
    3ij4_A-crystal       1.010000       0.789166        0.727758  0.724139
    4wit_A-crystal       0.980317       0.810587        0.793133  0.594927
    DDR1 [Q08345-ETRA]   0.978593       1.010000        0.837883  0.905371
    ErbB2 [P04626-ETRA]  0.974516       0.789166        0.727758  1.010000
    MPZ [P25189-ETRA]    0.909867       0.820061        0.822048  0.862866
    """

    #######################################################################################################
    #                                                                                                     #
    #                                       plot barchart                                                 #
    #                                                                                                     #
    #######################################################################################################
    # BO_barchart_png
    plt.close("all")
    #plt.rcParams.update({'font.size': 8})
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    # replace the protein names
    df_valid_indiv.index = pd.Series(df_valid_indiv.index).replace(namedict)
    df_valid_indiv.plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance value\n(observed overlap - random overlap)")
    ax.legend()#(["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(BO_barchart_png, dpi=240)


    #######################################################################################################
    #                                                                                                     #
    #                                plot linechart (combined data all proteins                           #
    #                                                                                                     #
    #######################################################################################################
    if plot_o_over_r:
        df_o_over_r = pd.read_excel(BO_data_excel, sheet_name="df_o_over_r", index_col=0)
        df_o_over_r_mean = df_o_over_r.T.mean()
    df_o_minus_r.columns = pd.Series(df_o_minus_r.columns).replace(namedict)
    df_o_minus_r_mean = df_o_minus_r.T.mean()
    # get the area under the curve
    AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)

    # BO_linechart_png
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    df_o_minus_r_mean.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC10 : {:0.2f}".format(AUBOC10))
    ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="random", alpha=0.5)

    if plot_o_over_r:
        ax2 = ax.twinx()
        df_o_over_r_mean.plot(ax=ax2, color="#9b2d0f", linestyle="-", label="old method (o/r)")
        ax2.plot([1, 10], [1, 1], color="#9b2d0f", linestyle="--", label="old method random", alpha=0.5)

    # ax.set_ylim(0)
    ax.grid(False)
    ax.set_ylabel("fraction of correctly predicted residues\n(observed - random)", color="#0f7d9b")
    ax.tick_params('y', colors="#0f7d9b")

    ax.spines['left'].set_color("#0f7d9b")
    ax.legend()
    if plot_o_over_r:
        ax2.tick_params('y', colors="#9b2d0f")
        ax2.spines['right'].set_color("#9b2d0f")
        #ax.set_ylabel("performance value\n (observed / random)", color="#9b2d0f")
        ax.set_ylabel("fraction of correctly predicted residues\n(observed / random)", color="#9b2d0f")
        ax2.legend()

    ax.set_xlabel("number of TMD residues\n(sample size)")
    fig.tight_layout()
    fig.savefig(BO_linechart_png, dpi=140)

    return AUBOC10


def save_extra_BO_figs(BO_data_excel, other_figs_path):
    linechart_mean_obs_and_rand = os.path.join(other_figs_path, "1_linechart_mean_obs_and_rand.png")
    linechart_obs_indiv = os.path.join(other_figs_path, "2_linechart_obs_indiv.png")
    linechart_p_indiv = os.path.join(other_figs_path, "3_linechart_p_indiv.png")
    linechart_o_minus_r = os.path.join(other_figs_path, "4_linechart_o_minus_r.png")
    linechart_o_over_r = os.path.join(other_figs_path, "5_linechart_o_over_r.png")

    dfrand = pd.read_excel(BO_data_excel, sheet_name="dfrand", index_col=0)
    dfobs = pd.read_excel(BO_data_excel, sheet_name="dfobs", index_col=0)
    df_o_minus_r = pd.read_excel(BO_data_excel, sheet_name="df_o_minus_r", index_col=0)
    # linechart_mean_obs_and_rand

    fig, ax = plt.subplots()
    dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
    dfobs.mean(axis=1).plot(ax=ax, color="k", label="mean observed")
    ax.grid(False)
    ax.set_ylabel("mean overlap")
    ax.legend()
    fig.savefig(linechart_mean_obs_and_rand, dpi=140)

    # linechart_obs_indiv

    plt.close("all")
    fig, ax = plt.subplots()
    dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
    dfobs.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper left", ncol=2)
    ax.set_ylabel("overlap")
    fig.savefig(linechart_obs_indiv, dpi=140)

    dfp = pd.read_excel(BO_data_excel, sheet_name="dfp", index_col=0)
    # linechart_p_indiv
    plt.close("all")
    fig, ax = plt.subplots()
    dfp.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper right", ncol=2)
    ax.set_ylabel("p-value of result")
    fig.savefig(linechart_p_indiv, dpi=140)

    # linechart_o_minus_r
    plt.close("all")
    fig, ax = plt.subplots()
    df_o_minus_r.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper left", ncol=2)
    ax.set_ylabel("observed - random")
    fig.savefig(linechart_o_minus_r, dpi=140)
