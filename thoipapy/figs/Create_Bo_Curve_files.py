import sys
import thoipapy
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
import glob
import numpy as np
from sklearn.externals import joblib
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

            odf = add_THOIPA_pred_to_combined_file(acc, train_df, test_df)
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
    if isinstance(s["test_datasets"], int):
        test_set_list = [s["test_datasets"]]
    else:
        test_set_list = s["test_datasets"].split(",")

    if isinstance(s["train_datasets"], int):
        train_set_list = [s["train_datasets"]]
    else:
        train_set_list = s["train_datasets"].split(",")

    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))
        #traindata_set = os.path.join(s["Result_folder"], trainsetname, "{}_train_data.csv".format(trainsetname))
        #train_df = pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)
        model_pkl = os.path.join(s["Result_folder"], trainsetname, "{}_rfmodel.pkl".format(trainsetname))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            BO_curve_csv = os.path.join(s["Bo_Curve_path"],"Train{}_Test{}.bocurve.csv".format(trainsetname, testsetname))

            # xlsx_list = glob.glob(os.path.join(s["set_path"], "{}*.xlsx".format(testsetname)))
            # if len(xlsx_list) == 1:
            #     testset_path = xlsx_list[0]
            # elif len(xlsx_list) == 0:
            #     raise FileNotFoundError(
            #         "Excel file with this test data set not found.\nsetname = {}\nexcel files in folder = {}".format(
            #             testsetname, xlsx_list))
            # elif len(xlsx_list) > 1:
            #     raise ValueError(
            #         "More than one excel file in set folder contains '{}' in the filename.\nexcel files in folder = {}".format(
            #             testsetname, xlsx_list))

            testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["set_path"])
            #train_features_del = ["residue_num", "residue_name", "acc_db", "n_homologues", "interface_score", "bind"]
            #test_features_del = ["residue_num", "residue_name", "n_homologues", "bind", "Disruption"]

            testdataset_df = pd.read_excel(testset_path, sheetname="proteins")
            acc_list = testdataset_df.acc.tolist()
            database = testdataset_df.database[0]
            dfc = pd.DataFrame()
            for acc in acc_list:
                testdata_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database,
                                                      "{}.surr20.gaps5.combined_features.csv".format(acc))
                test_combined_incl_pred = os.path.join(os.path.dirname(s["thoipapy_feature_folder"]), "combined_incl_pred", database,
                                                      "{}.train{}.combined_features_incl_pred.csv".format(acc, trainsetname))
                thoipapy.utils.make_sure_path_exists(test_combined_incl_pred, isfile=True)

                add_THOIPA_pred_to_combined_file(model_pkl, testdata_combined_file, test_combined_incl_pred)

                test_df = pd.read_csv(test_combined_incl_pred)

                df_for_best_overlap = test_df[["residue_name", "interface_score", "THOIPA"]].copy()

                #test_df.index = test_df.index.astype(int) + 1

                #test_df["interface_score"] = -1 * test_df["interface_score"]

                if database == "crystal" or database == "NMR":
                    df_for_best_overlap["interface_score"] = -1 * df_for_best_overlap["interface_score"]
                    #interface_score = test_df.interface_score
                    #interface_score = -1 * interface_score  # (it is closest distance and low value means high propencity of interfacial)

                #if database == "crystal" or database == "NMR":
                #    interface_score = -1 * interface_score  # (it is closest distance and low value means high propencity of interfacial)
                #elif database == "ETRA":
                #    pass  # (it is closest experimental disruption and high value means high propencity of interfacial)
                #else:
                #    raise ValueError()

                #prob_arr = test_df["THOIPA"].as_matrix()
                #interface_score = test_df["interface_score"].as_matrix()


                odf = thoipapy.figs.fig_utils.calc_best_overlap(acc, df_for_best_overlap)

                if dfc.empty:
                    dfc = odf
                else:
                    dfc = pd.concat([dfc, odf], axis=1, join="outer")

            dfc.to_csv(BO_curve_csv)
            sys.stdout.write("\npred_interf_single_prot_using_sel_train_datasets finished ({})".format(BO_curve_csv))
            sys.stdout.flush()

    sys.stdout.write("\npred_interf_single_prot_using_sel_train_datasets finished (all train and test datasets)".format(BO_curve_csv))
    sys.stdout.flush()

def add_THOIPA_pred_to_combined_file(model_pkl, testdata_combined_file, test_combined_incl_pred):

    test_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)

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
    tX = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(test_df)

    prob_arr = fit.predict_proba(tX)[:, 1]

    # if hasattr(clf,'predict_proba'):
    #     prob_arr = fit.predict_proba(tX)[:,1]
    # else:
    #     prob_arr = fit.decision_function(tX)
    #     prob_arr = (prob_arr - prob_arr.min())/(prob_arr.max() - prob_arr.min())

    test_df["THOIPA"] = prob_arr
    test_df.to_csv(test_combined_incl_pred)


#
# def add_THOIPA_pred_to_combined_file(acc, train_df, test_df, database):
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
#     odf = thoipapy.figs.fig_utils.calc_best_overlap(acc, prob_pos, interface_score, database)
#
#     return odf