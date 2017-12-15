import sys
import thoipapy
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
import glob
import numpy as np
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
import pickle
import warnings
warnings.filterwarnings("ignore")

def Test_Etra(s):
    testsetname = "set03"
    train_set_list = s["train_datasets"].split(",")
    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))
        traindata_set = os.path.join(s["Result_folder"], trainsetname, "{}_train_data.csv".format(trainsetname))
        train_df =  pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)

        if not os.path.exists(s["Bo_Curve_path"]):
            os.makedirs(s["Bo_Curve_path"])

        BO_curve_csv = os.path.join(s["Bo_Curve_path"], "Train{}_Test{}.bocurve.csv".format(trainsetname, testsetname))

        train_features_del = ["residue_num","residue_name","acc_db","n_homologues","interface_score","bind"]
        test_features_del=["residue_num","residue_name","n_homologues","bind","Disruption"]

        # xlsx_list = glob.glob(os.path.join(s["set_path"], "{}*.xlsx".format(testsetname)))
        # if len(xlsx_list) == 1:
        #     testset_path = xlsx_list[0]
        # elif len(xlsx_list) == 0:
        #     raise FileNotFoundError(
        #     "Excel file with this test data set not found.\nsetname = {}\nexcel files in folder = {}".format(testsetname, xlsx_list))
        # elif len(xlsx_list) > 1:
        #     raise ValueError(
        #         "More than one excel file in set folder contains '{}' in the filename.\nexcel files in folder = {}".format(
        #             testsetname, xlsx_list))

        testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["set_path"])

        testdataset_df = pd.read_excel(testset_path,sheetname="proteins")
        acc_list = testdataset_df.acc.tolist()
        database = testdataset_df.database[0]
        dfc = pd.DataFrame()
        for acc in acc_list:
            testdata_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
            test_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)

            odf = save_THOIPA_pred_indiv_prot(acc, train_df, test_df)
            if dfc.empty:
                dfc = odf
            else:
                dfc = pd.concat([dfc, odf], axis=1, join="outer")

        dfc.to_csv(BO_curve_csv)



def pred_interf_single_prot_using_sel_train_datasets(s):
    if not os.path.exists(s["Bo_Curve_path"]):
        os.makedirs(s["Bo_Curve_path"])

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset

    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    for n, train_set in enumerate(train_set_list):
        trainsetname = "set{:02d}".format(int(train_set))
        #traindata_set = os.path.join(s["Result_folder"], trainsetname, "{}_train_data.csv".format(trainsetname))
        #train_df = pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)
        model_pkl = os.path.join(s["Result_folder"], trainsetname, "{}_rfmodel.pkl".format(trainsetname))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            THOIPA_BO_curve_data_csv = os.path.join(s["Bo_Curve_path"],"Test{}_Train{}.THOIPA.best_overlap_data.csv".format(testsetname, trainsetname))
            #LIPS_BO_curve_data_csv = os.path.join(s["Bo_Curve_path"], "Test{}.LIPS.best_overlap_data.csv".format(testsetname, trainsetname))
            THOIPA_ROC_pkl = os.path.join(s["Bo_Curve_path"], "Test{}_Train{}.THOIPA.ROC_data.pkl".format(testsetname, trainsetname))

            testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["set_path"])
            #train_features_del = ["residue_num", "residue_name", "acc_db", "n_homologues", "interface_score", "bind"]
            #test_features_del = ["residue_num", "residue_name", "n_homologues", "bind", "Disruption"]

            testdataset_df = pd.read_excel(testset_path, sheetname="proteins")
            acc_list = testdataset_df.acc.tolist()
            database = testdataset_df.database[0]
            THOIPA_BO_data_df = pd.DataFrame()
            #LIPS_BO_data_df = pd.DataFrame()

            # save all outputs to a cross-validation dictionary, to be saved as a pickle file
            xv_dict_THOIPA = {}
            mean_tpr = 0.0
            mean_fpr = np.linspace(0, 1, 100)

            for acc in acc_list:
                testdata_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database,
                                                      "{}.surr20.gaps5.combined_features.csv".format(acc))
                THOIPA_pred_csv = os.path.join(os.path.dirname(s["thoipapy_feature_folder"]), "Predictions", "testset_trainset", database,
                                                      "{}.THOIPA.train{}.csv".format(acc, trainsetname))
                combined_incl_THOIPA_csv = os.path.join(os.path.dirname(s["thoipapy_feature_folder"]), "Predictions", "testset_trainset", database,
                                                      "{}.THOIPA_incl_combined.train{}.csv".format(acc, trainsetname))
                thoipapy.utils.make_sure_path_exists(combined_incl_THOIPA_csv, isfile=True)

                combined_incl_THOIPA_df = save_THOIPA_pred_indiv_prot(model_pkl, testdata_combined_file, THOIPA_pred_csv, combined_incl_THOIPA_csv)

                #######################################################################################################
                #                                                                                                     #
                #                           Processing BO curve data for each single protein                          #
                #                                                                                                     #
                #######################################################################################################
                # SAVE LIPS PREDICTION DATA
                # this is somewhat inefficient, as it is conducted for every test dataset
                #LIPS_pred_csv = os.path.join(os.path.dirname(s["thoipapy_feature_folder"]), "Predictions", "testset_trainset", database,
                #                                      "{}.LIPS_pred.csv".format(acc, trainsetname))
                #LIPS_pred_df = combined_incl_THOIPA_df[["residue_name", "residue_num", "LIPS_lipo", "LIPS_entropy", "LIPS_L*E", "LIPS_surface"]]
                #LIPS_pred_df.to_csv(LIPS_pred_csv)

                combined_incl_THOIPA_df["LIPS_L*E"] = -1 * combined_incl_THOIPA_df["LIPS_L*E"]

                if database == "crystal" or database == "NMR":
                    # (it is closest distance and low value means high propencity of interfacial)
                    combined_incl_THOIPA_df["interface_score"] = -1 * combined_incl_THOIPA_df["interface_score"]

                THOIPA_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc, combined_incl_THOIPA_df)
                #LIPS_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc, combined_incl_THOIPA_df, experiment_col="interface_score", pred_col="LIPS_L*E")

                if THOIPA_BO_data_df.empty:
                    THOIPA_BO_data_df = THOIPA_BO_single_prot_df
                    #LIPS_BO_data_df = LIPS_BO_single_prot_df
                else:
                    THOIPA_BO_data_df = pd.concat([THOIPA_BO_data_df, THOIPA_BO_single_prot_df], axis=1, join="outer")
                    #LIPS_BO_data_df = pd.concat([LIPS_BO_data_df, LIPS_BO_single_prot_df], axis=1, join="outer")

                #######################################################################################################
                #                                                                                                     #
                #                     Processing ROC data for each single protein, saving to nested dict              #
                #                                                                                                     #
                #######################################################################################################

                df_for_roc = combined_incl_THOIPA_df.dropna(subset=["interface_score"])

                predictor_name = "THOIPA"

                fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name])
                auc_value = auc(fpr, tpr)
                mean_tpr += interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0

                xv_dict_THOIPA[acc] = {"fpr" : fpr, "tpr" : tpr, "auc" : auc_value}

            #######################################################################################################
            #                                                                                                     #
            #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
            #                                                                                                     #
            #######################################################################################################

            THOIPA_BO_data_df.to_csv(THOIPA_BO_curve_data_csv)
            #LIPS_BO_data_df.to_csv(LIPS_BO_curve_data_csv)
            names_excel_path = os.path.join(os.path.dirname(s["set_path"]), "ETRA_NMR_names.xlsx")

            THOIPA_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(THOIPA_BO_curve_data_csv, names_excel_path)
            #LIPS_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(LIPS_BO_curve_data_csv, names_excel_path)

            sys.stdout.write("\nBO curve data analysed ({})".format(THOIPA_linechart_mean_obs_and_rand))

            #######################################################################################################
            #                                                                                                     #
            #                     Processing dictionary with ROC data, saving to pickle                           #
            #                                                                                                     #
            #######################################################################################################
            mean_tpr /= len(acc_list)
            mean_tpr[-1] = 1.0

            mean_auc = auc(mean_fpr, mean_tpr)

            ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_THOIPA}
            ROC_out_dict["true_positive_rate_mean"] = mean_tpr
            ROC_out_dict["false_positive_rate_mean"] = mean_fpr
            ROC_out_dict["mean_auc"] = mean_auc

            # save dict as pickle
            with open(THOIPA_ROC_pkl, "wb") as f:
                pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

            sys.stdout.write("\nBO curve and ROC data collected. Mean AUC = {:.03f}. BO output={}".format(mean_auc, THOIPA_BO_curve_data_csv))

            create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl)

            sys.stdout.write("\nTest{}_Train{} finished\n".format(testsetname, trainsetname))
            sys.stdout.flush()

    validate_LIPS_for_testset(s)

    sys.stdout.write("\npred_interf_single_prot_using_sel_train_datasets finished (all train and test datasets)")
    sys.stdout.flush()



def validate_LIPS_for_testset(s):
    if not os.path.exists(s["Bo_Curve_path"]):
        os.makedirs(s["Bo_Curve_path"])

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset

    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    for test_set in test_set_list:
        testsetname = "set{:02d}".format(int(test_set))
        LIPS_BO_curve_data_csv = os.path.join(s["Bo_Curve_path"], "Test{}.LIPS.best_overlap_data.csv".format(testsetname))
        LIPS_ROC_pkl = os.path.join(s["Bo_Curve_path"], "Test{}.LIPS.ROC_data.pkl".format(testsetname))

        testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["set_path"])
        #train_features_del = ["residue_num", "residue_name", "acc_db", "n_homologues", "interface_score", "bind"]
        #test_features_del = ["residue_num", "residue_name", "n_homologues", "bind", "Disruption"]

        testdataset_df = pd.read_excel(testset_path, sheetname="proteins")
        acc_list = testdataset_df.acc.tolist()
        database = testdataset_df.database[0]
        LIPS_BO_data_df = pd.DataFrame()

        # save all outputs to a cross-validation dictionary, to be saved as a pickle file
        xv_dict_LIPS = {}
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        for acc in acc_list:
            testdata_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database,
                                                  "{}.surr20.gaps5.combined_features.csv".format(acc))

            combined_df = pd.read_csv(testdata_combined_file, index_col=0)

            #######################################################################################################
            #                                                                                                     #
            #                           Processing BO curve data for each single protein                          #
            #                                                                                                     #
            #######################################################################################################
            # SAVE LIPS PREDICTION DATA
            # this is somewhat inefficient, as it is conducted for every test dataset
            LIPS_pred_csv = os.path.join(os.path.dirname(s["thoipapy_feature_folder"]), "Predictions", "testset_trainset", database,
                                                  "{}.LIPS_pred.csv".format(acc, testsetname))
            LIPS_pred_df = combined_df[["residue_name", "residue_num", "LIPS_lipo", "LIPS_entropy", "LIPS_L*E", "LIPS_surface"]]
            LIPS_pred_df.to_csv(LIPS_pred_csv)

            combined_df["LIPS_L*E"] = -1 * combined_df["LIPS_L*E"]

            if database == "crystal" or database == "NMR":
                # (it is closest distance and low value means high propencity of interfacial)
                combined_df["interface_score"] = -1 * combined_df["interface_score"]

            LIPS_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc, combined_df, experiment_col="interface_score", pred_col="LIPS_L*E")

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

            predictor_name = "LIPS_L*E"

            fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name])
            auc_value = auc(fpr, tpr)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0

            xv_dict_LIPS[acc] = {"fpr" : fpr, "tpr" : tpr, "auc" : auc_value}

        #######################################################################################################
        #                                                                                                     #
        #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
        #                                                                                                     #
        #######################################################################################################

        LIPS_BO_data_df.to_csv(LIPS_BO_curve_data_csv)
        names_excel_path = os.path.join(os.path.dirname(s["set_path"]), "ETRA_NMR_names.xlsx")

        LIPS_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(LIPS_BO_curve_data_csv, names_excel_path)

        sys.stdout.write("\nBO curve data analysed ({})".format(LIPS_linechart_mean_obs_and_rand))

        #######################################################################################################
        #                                                                                                     #
        #                     Processing dictionary with ROC data, saving to pickle                           #
        #                                                                                                     #
        #######################################################################################################
        mean_tpr /= len(acc_list)
        mean_tpr[-1] = 1.0

        mean_auc = auc(mean_fpr, mean_tpr)

        ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_LIPS}
        ROC_out_dict["true_positive_rate_mean"] = mean_tpr
        ROC_out_dict["false_positive_rate_mean"] = mean_fpr
        ROC_out_dict["mean_auc"] = mean_auc

        # save dict as pickle
        with open(LIPS_ROC_pkl, "wb") as f:
            pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

        sys.stdout.write("\nBO curve and ROC data collected. Mean AUC = {:.03f}. BO output={}".format(mean_auc, LIPS_BO_curve_data_csv))

        create_ROC_fig_for_testset_trainset_combination(LIPS_ROC_pkl)

        sys.stdout.write("\nTest{} LIPS validation finished\n".format(testsetname))
        sys.stdout.flush()

    sys.stdout.write("\npred_interf_single_prot_using_sel_train_datasets finished (all train and test datasets)")
    sys.stdout.flush()



def create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl):

    plt.rcParams.update({'font.size': 6})
    ROC_pkl_foldername = os.path.basename(THOIPA_ROC_pkl)[:-4]
    ROC_pkl_dir = os.path.dirname(THOIPA_ROC_pkl)

    ROC_png = os.path.join(ROC_pkl_dir, ROC_pkl_foldername, "{}.ROC.png".format(ROC_pkl_foldername))
    thoipapy.utils.make_sure_path_exists(ROC_png, isfile=True)

    # open pickle file
    with open(THOIPA_ROC_pkl, "rb") as f:
        ROC_out_dict = pickle.load(f)

    xv_dict_THOIPA = ROC_out_dict["xv_dict_THOIPA"]

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    for acc in xv_dict_THOIPA:
        roc_auc = xv_dict_THOIPA[acc]["auc"]
        ax.plot(xv_dict_THOIPA[acc]["fpr"], xv_dict_THOIPA[acc]["tpr"], lw=1, label='{} ({:0.2f})'.format(acc, roc_auc), alpha=0.8)

    #mean_auc = auc(df_xv["false_positive_rate"], df_xv["true_positive_rate"])
    mean_auc = ROC_out_dict["mean_auc"]

    ax.plot(ROC_out_dict["false_positive_rate_mean"], ROC_out_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_auc, lw=1.5)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_png, dpi=240)
    fig.savefig(ROC_png[:-4] + ".pdf")

def save_THOIPA_pred_indiv_prot(model_pkl, testdata_combined_file, THOIPA_pred_csv, test_combined_incl_pred):

    combined_incl_THOIPA_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)

    #drop_cols_not_used_in_ML
    #X=train_df.drop(train_features_del,axis=1)
    # X = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(train_df)
    # y = train_df["interface"]
    # clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    # fit = clf.fit(X,y)

    fit = joblib.load(model_pkl)

    # if database == "ETRA":
    #     dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
    #     ddf = pd.read_excel(dirupt_path, index_col=0)
    #     interface_score = ddf.Disruption
    #     interface_score = interface_score      #(it is closest experimental disruption and high value means high propencity of interfacial)

    #Lips_score = test_df.LIPS_lipo * test_df.LIPS_entropy

    #tX=test_df.drop(test_features_del,axis=1)
    test_X = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(combined_incl_THOIPA_df)

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
#     X = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(train_df)
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
#     Lips_score = test_df.LIPS_lipo * test_df.LIPS_entropy
#
#     #tX=test_df.drop(test_features_del,axis=1)
#     tX = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(test_df)
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


def analyse_bo_curve_underlying_data(bo_data_csv, names_excel_path):
    """Analyse the Bo-curve underlying data.

    Parses data into more easily manipulated dataframes.
    Creates figures if desired.
    All figures are saved in a subfolder with the same name as the original file.

    Parameters
    ----------
    bo_data_csv : str
        Path to the csv created by bo, with the underlying data.
        Index : Top1, Top2, Ono, Pno, Rno etc
    names_excel_path : str
        Path to the excel file with the protein short names and reference.

    Usage
    -----
    import datoxr
    bo_data_csv = r"D:\drive\TMD_homodimer\figs\SuppDataX02-best_overlap_data\SuppDataX02.csv"
    names_excel_path = r"D:\drive\TMD_homodimer\data_xy\ETRA_NMR_names.xlsx"
    datoxr.figs.bo_curve_analysis.analyse_bo_curve_underlying_data(bo_data_csv, names_excel_path)
    """

    # change to empty list if you don't want to create figures. Or [7] if you only want to process figure 7, for example.
    list_figs_to_create = range(1, 11)
    #list_figs_to_create = [7]

    out_folder = os.path.join(os.path.dirname(bo_data_csv), os.path.basename(bo_data_csv[:-4]))
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # create output paths
    excel_out_path = os.path.join(out_folder, "bo_curve_underlying_data_indiv_df.xlsx")
    linechart_mean_obs_and_rand = os.path.join(out_folder, "1_linechart_mean_obs_and_rand.png")
    linechart_obs_indiv = os.path.join(out_folder, "2_linechart_obs_indiv.png")
    linechart_p_indiv = os.path.join(out_folder, "3_linechart_p_indiv.png")
    linechart_o_minus_r = os.path.join(out_folder, "4_linechart_o_minus_r.png")
    linechart_o_over_r = os.path.join(out_folder, "5_linechart_o_over_r.png")
    linechart_method_comparison_minus_vs_over = os.path.join(out_folder, "6_linechart_method_comparison_minus_vs_over.png")
    barchart_ss5_ss10_indiv_prot = os.path.join(out_folder, "7_barchart_ss5_ss10_indiv_prot.png")


    dfb = pd.read_csv(bo_data_csv, index_col=0)

    """ORIGINAL BO DATA CSV LOOKS LIKE THIS
    Top1 = sample size 1
    Ono = overlap in data
    Rno = random overlap based on that sequence length and sample size
    Pono = p-value for finding that overlap

         Unnamed: 1  O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286  Ratio (Average(Ono)/Average(Rno))
    NaN         Ono    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    1.00    0.00    0.00    1.00    1.00    1.00    0.00    0.00    0.00                                NaN
    Top1        Rno    0.05    0.04    0.05    0.04    0.05    0.06    0.04    0.04    0.05    0.04    0.04    0.05    0.06    0.05    0.06    0.04    0.05                               5.01
    NaN        Pono    0.95    0.96    0.95    0.96    0.95    0.94    0.96    0.96    0.05    0.96    0.96    0.05    0.06    0.05    0.94    0.96    0.95                                NaN
    NaN         Ono    1.00    1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    1.00    0.00    1.00    2.00    2.00    0.00    0.00    1.00                                NaN
    Top2        Rno    0.19    0.17    0.21    0.15    0.19    0.24    0.17    0.17    0.19    0.17    0.16    0.19    0.25    0.20    0.24    0.17    0.20                               3.68
    """

    # create an index based on sample size [1 1 1 2 2 2  etc..
    ind = []
    for i in range(1, int((len(dfb) / 3)) + 1):
        ind += list(np.array([1, 1, 1]) * i)
    dfb.index = ind
    dfb.index.name = "sample size"

    """NOW INDICES ARE BASED ON SAMPLE SIZE

                Unnamed: 1  O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286  Ratio (Average(Ono)/Average(Rno))
    sample size                                                                                                                                                                                      
    1                  Ono    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    1.00    0.00    0.00    1.00    1.00    1.00    0.00    0.00    0.00                                NaN
    1                  Rno    0.05    0.04    0.05    0.04    0.05    0.06    0.04    0.04    0.05    0.04    0.04    0.05    0.06    0.05    0.06    0.04    0.05                               5.01
    1                 Pono    0.95    0.96    0.95    0.96    0.95    0.94    0.96    0.96    0.05    0.96    0.96    0.05    0.06    0.05    0.94    0.96    0.95                                NaN
    2                  Ono    1.00    1.00    0.00    1.00    0.00    0.00    1.00    0.00    1.00    1.00    0.00    1.00    2.00    2.00    0.00    0.00    1.00                                NaN
    2                  Rno    0.19    0.17    0.21    0.15    0.19    0.24    0.17    0.17    0.19    0.17    0.16    0.19    0.25    0.20    0.24    0.17    0.20      
    """

    # split into separate dataframes
    # dataframe of observed overlaps
    dfobs = dfb.iloc[::3, 1:-1].astype(int)
    # dataframe of random calculated overlaps
    dfrand = dfb.iloc[1::3, 1:-1]
    # dataframe of p-values
    dfp = dfb.iloc[2::3, 1:-1]

    """FOR EXAMPLE df_obs NOW LOOKS LIKE THIS, WITH A ROW FOR EACH SAMPLE SIZE:

                 O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286
    sample size                                                                                                                                        
    1                 0       0       0       0       0       0       0       0       1       0       0       1       1       1       0       0       0
    2                 1       1       0       1       0       0       1       0       1       1       0       1       2       2       0       0       1
    3                 1       1       1       2       1       0       1       0       2       2       0       2       2       2       0       0       2
    4                 2       2       1       2       1       0       2       1       3       2       1       3       2       2       0       1       2
    5                 3       3       2       3       2       1       2       2       4       4       2       4       3       2       1       2       3"""

    df_o_minus_r = dfobs - dfrand
    df_o_over_r = dfobs / dfrand

    """df_o_minus_r is negative where the result is lower than random

                 O75460  P02724  P05106  P06583  P08514  P0A6S5  P23470  P35590  Q08345  Q12983  Q16827  Q16832  Q6ZRP7  Q7L4S7  Q8NI60  Q92729  Q9Y286
    sample size                                                                                                                                        
    1             -0.05   -0.04   -0.05   -0.04   -0.05   -0.06   -0.04   -0.04    0.95   -0.04   -0.04    0.95    0.94    0.95   -0.06   -0.04   -0.05
    2              0.81    0.83   -0.21    0.85   -0.19   -0.24    0.83   -0.17    0.81    0.83   -0.16    0.81    1.75    1.80   -0.24   -0.17    0.80
    3              0.57    0.61    0.53    1.65    0.57   -0.53    0.61   -0.37    1.57    1.63   -0.36    1.57    1.44    1.55   -0.53   -0.37    1.55
    4              1.24    1.30    0.16    1.38    0.24   -0.94    1.30    0.33    2.24    1.33    0.36    2.24    1.00    1.20   -0.94    0.33    1.20
    5              1.81    1.91    0.68    2.04    0.81   -0.47    0.91    0.96    2.81    2.96    1.00    2.81    1.44    0.75   -0.47    0.96    1.75


    df_o_over_r is where the result is lower than random.
    This is not quite a fair comparison, as the zeros are caused by 0 overlap / signigicant random overlap

                   O75460    P02724    P05106    P06583    P08514    P0A6S5    P23470    P35590     Q08345    Q12983  Q16827     Q16832     Q6ZRP7     Q7L4S7    Q8NI60    Q92729    Q9Y286
    sample size                                                                                                                                                                            
    1            0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  20.000000  0.000000  0.0000  20.000000  16.666667  20.000000  0.000000  0.000000  0.000000
    2            5.263158  5.882353  0.000000  6.666667  0.000000  0.000000  5.882353  0.000000   5.263158  5.882353  0.0000   5.263158   8.000000  10.000000  0.000000  0.000000  5.000000
    3            2.325581  2.564103  2.127660  5.714286  2.325581  0.000000  2.564103  0.000000   4.651163  5.405405  0.0000   4.651163   3.571429   4.444444  0.000000  0.000000  4.444444
    4            2.631579  2.857143  1.190476  3.225806  1.315789  0.000000  2.857143  1.492537   3.947368  2.985075  1.5625   3.947368   2.000000   2.500000  0.000000  1.492537  2.500000
    5            2.521008  2.752294  1.515152  3.125000  1.680672  0.680272  1.834862  1.923077   3.361345  3.846154  2.0000   3.361345   1.923077   1.600000  0.680272  1.923077  2.400000
    """

    #################################################################
    #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    #################################################################

    with pd.ExcelWriter(excel_out_path) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")

    #################################################################
    #             EXTRACT NAMES FROM NAMES EXCEL FILE               #
    #################################################################
    df_names = pd.read_excel(names_excel_path, index_col=0)
    # restrict names dict to only that database
    database = "ETRA"
    df_names = df_names.loc[df_names.database == database]
    df_names["label"] = df_names.shortname + " [" + df_names.index + "]"
    namedict = df_names["label"].to_dict()
    df_o_minus_r.columns = pd.Series(df_o_minus_r.columns).replace(namedict)

    # linechart_mean_obs_and_rand
    fignr = 1
    if fignr in list_figs_to_create:
        fig, ax = plt.subplots()
        dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
        dfobs.mean(axis=1).plot(ax=ax, color="k", label="mean observed")
        ax.grid(False)
        ax.set_ylabel("mean overlap")
        ax.legend()
        fig.savefig(linechart_mean_obs_and_rand, dpi=140)

    # linechart_obs_indiv
    fignr = 2
    if fignr in list_figs_to_create:
        plt.close("all")
        fig, ax = plt.subplots()
        dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
        dfobs.plot(ax=ax, alpha=0.7)
        ax.legend(loc="upper left", ncol=2)
        ax.set_ylabel("overlap")
        fig.savefig(linechart_obs_indiv, dpi=140)

    # linechart_p_indiv
    fignr = 3
    if fignr in list_figs_to_create:
        plt.close("all")
        fig, ax = plt.subplots()
        dfp.plot(ax=ax, alpha=0.7)
        ax.legend(loc="upper right", ncol=2)
        ax.set_ylabel("p-value of result")
        fig.savefig(linechart_p_indiv, dpi=140)

    # linechart_o_minus_r
    fignr = 4
    if fignr in list_figs_to_create:
        plt.close("all")
        fig, ax = plt.subplots()
        df_o_minus_r.plot(ax=ax, alpha=0.7)
        ax.legend(loc="upper left", ncol=2)
        ax.set_ylabel("observed - random")
        fig.savefig(linechart_o_minus_r, dpi=140)

    # linechart_o_over_r
    fignr = 5
    if fignr in list_figs_to_create:
        plt.close("all")
        fig, ax = plt.subplots()
        df_o_over_r.plot(ax=ax, alpha=0.7)
        ax.set_ylabel("observed / random")
        ax.legend(loc="upper left", ncol=2)
        fig.savefig(linechart_o_over_r, dpi=140)

    # linechart_method_comparison_minus_vs_over
    fignr = 6
    if fignr in list_figs_to_create:
        plt.close("all")
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        df_o_minus_r.T.mean().plot(ax=ax, color="#0f7d9b", linestyle="-", label="mean observed overlap - random overlap")
        df_o_over_r.T.mean().plot(ax=ax2, color="#9b2d0f", linestyle="--", label="mean observed overlap / random overlap")
        #ax.set_ylim(0)
        ax.grid(False)
        ax.set_ylabel("performance value", color="#0f7d9b")
        ax2.set_ylabel("performance value", color="#9b2d0f")
        ax.legend()
        fig.tight_layout()
        fig.savefig(linechart_method_comparison_minus_vs_over, dpi=140)

    # barchart_ss5_ss10_indiv_prot
    fignr = 7
    if fignr in list_figs_to_create:
        plt.close("all")
        plt.rcParams.update({'font.size': 8})
        fig, ax = plt.subplots(figsize=(3.42,3.42))
        df_o_minus_r_sel = df_o_minus_r.loc[[5, 10], :].T
        df_o_minus_r_sel.sort_values(5, axis=0, ascending=False, inplace=True)
        df_o_minus_r_sel.plot(kind="bar", ax=ax, alpha=0.7)
        ax.set_ylabel("performance value\n(observed overlap - random overlap)")
        ax.legend(["sample size = 5", "sample size = 10"])
        fig.tight_layout()
        ax.grid(False)
        fig.savefig(barchart_ss5_ss10_indiv_prot, dpi=240)

    return linechart_mean_obs_and_rand