import sys
import thoipapy
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
import glob
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def Test_Etra(s):
    testsetname = "set03"
    train_set_list = s["train_datasets"].split(",")
    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))
        traindata_set = os.path.join(s["Result_folder"], "{}\{}_train_data.csv".format(trainsetname,trainsetname))
        traindata_set_df =  pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)

        if not os.path.exists(s["Bo_Curve_path"]):
            os.makedirs(s["Bo_Curve_path"])

        BoCurve_Output_path = os.path.join(s["Bo_Curve_path"],
                                               "Train{}_Test{}.bocurve.csv".format(trainsetname, testsetname))

        train_features_del = ["residue_num","residue_name","acc_db","n_homologues","interface_score","bind"]
        test_features_del=["residue_num","residue_name","n_homologues","bind","Disruption"]

        xlsx_list = glob.glob(os.path.join(s["set_path"], "{}*.xlsx".format(testsetname)))
        if len(xlsx_list) == 1:
            testset_path = xlsx_list[0]
        elif len(xlsx_list) == 0:
            raise FileNotFoundError(
            "Excel file with this test data set not found.\nsetname = {}\nexcel files in folder = {}".format(testsetname, xlsx_list))
        elif len(xlsx_list) > 1:
            raise ValueError(
                "More than one excel file in set folder contains '{}' in the filename.\nexcel files in folder = {}".format(
                    testsetname, xlsx_list))

        testdataset_df = pd.read_excel(testset_path,sheetname="proteins")
        acc_list = testdataset_df.acc.tolist()
        database = testdataset_df.database[0]
        dfc = pd.DataFrame()
        for acc in acc_list:
            testdata_combined_file = os.path.join(s["thoipapy_feature_folder"],"combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))
            test_data_df = pd.read_csv(testdata_combined_file,sep=',',engine='python',index_col=0)
            odf=Train_Input_Test_Pred_Out(acc, traindata_set_df, test_data_df, train_features_del, test_features_del, database, s)
            if dfc.empty:
                dfc = odf
            else:
                dfc = pd.concat([dfc, odf], axis=1, join="outer")

        dfc.to_csv(BoCurve_Output_path)



def pred_interf_single_prot_using_sel_train_datasets(s):
    if not os.path.exists(s["Bo_Curve_path"]):
        os.makedirs(s["Bo_Curve_path"])

    test_set_list = s["test_datasets"].split(",")
    train_set_list = s["train_datasets"].split(",")
    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))
        traindata_set = os.path.join(s["Result_folder"], "{}\{}_train_data.csv".format(trainsetname,trainsetname))
        traindata_set_df =  pd.read_csv(traindata_set, sep=',', engine='python', index_col=0)

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            BoCurve_Output_path = os.path.join(s["Bo_Curve_path"],
                                           "Train{}_Test{}.bocurve.csv".format(trainsetname, testsetname))

            xlsx_list = glob.glob(os.path.join(s["set_path"], "{}*.xlsx".format(testsetname)))
            if len(xlsx_list) == 1:
                testset_path = xlsx_list[0]
            elif len(xlsx_list) == 0:
                raise FileNotFoundError(
                    "Excel file with this test data set not found.\nsetname = {}\nexcel files in folder = {}".format(
                        testsetname, xlsx_list))
            elif len(xlsx_list) > 1:
                raise ValueError(
                    "More than one excel file in set folder contains '{}' in the filename.\nexcel files in folder = {}".format(
                        testsetname, xlsx_list))

            train_features_del = ["residue_num", "residue_name", "acc_db", "n_homologues", "interface_score", "bind"]
            test_features_del = ["residue_num", "residue_name", "n_homologues", "bind", "Disruption"]

            testdataset_df = pd.read_excel(testset_path, sheetname="proteins")
            acc_list = testdataset_df.acc.tolist()
            database = testdataset_df.database[0]
            dfc = pd.DataFrame()
            for acc in acc_list:
                testdata_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database,
                                                      "{}.surr20.gaps5.combined_features.csv".format(acc))
                test_data_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)
                odf = Train_Input_Test_Pred_Out(acc, traindata_set_df, test_data_df, train_features_del, test_features_del,
                                                database, s)
                if dfc.empty:
                    dfc = odf
                else:
                    dfc = pd.concat([dfc, odf], axis=1, join="outer")

            dfc.to_csv(BoCurve_Output_path)
    sys.stdout.write("pred_interf_single_prot_using_sel_train_datasets finished ({})".format(BoCurve_Output_path))
    sys.stdout.flush()




def Train_Input_Test_Pred_Out(acc, train_data_df, tdf, train_features_del, test_features_del, database, s):
    #drop_cols_not_used_in_ML
    #X=train_data_df.drop(train_features_del,axis=1)
    X = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(train_data_df)
    y = train_data_df["interface"]
    clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    fit = clf.fit(X,y)

    interface_score = tdf.interface_score

    if database == "crystal" or database == "NMR":
        interface_score = tdf.interface_score
        interface_score = -1 * interface_score  #(it is closest distance and low value means high propencity of interfacial)
    # if database == "ETRA":
    #     dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
    #     ddf = pd.read_excel(dirupt_path, index_col=0)
    #     interface_score = ddf.Disruption
    #     interface_score = interface_score      #(it is closest experimental disruption and high value means high propencity of interfacial)
    tdf.index = tdf.index.astype(int) +1
    Lips_score = tdf.LIPS_lipo * tdf.LIPS_entropy

    #tX=tdf.drop(test_features_del,axis=1)
    tX = thoipapy.RF_features.RF_Train_Test.drop_cols_not_used_in_ML(tdf)

    if hasattr(clf,'predict_proba'):
        prob_pos=fit.predict_proba(tX)[:,1]
    else:
        prob_pos=fit.decision_function(tX)
        prob_pos=(prob_pos - prob_pos.min())/(prob_pos.max() - prob_pos.min())
    odf=thoipapy.figs.fig_utils.Bo_Curve_Create(acc,prob_pos,Lips_score,interface_score,database)
    return odf