import thoipapy
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
import glob
import numpy as np



def Test_Etra(s):
    crystalnmr_train_data = os.path.join(s["Result_folder"], "set04\set04_train_data.csv")
    crystalnmr_train_data_df = pd.read_csv(crystalnmr_train_data, sep=',', engine='python', index_col=0)
    crystal_train_data = os.path.join(s["Result_folder"], "set01\set01_train_data.csv")
    crystal_train_data_df = pd.read_csv(crystal_train_data, sep=',', engine='python', index_col=0)
    nmr_train_data = os.path.join(s["Result_folder"], "set02\set02_train_data.csv")
    nmr_train_data_df = pd.read_csv(nmr_train_data, sep=',', engine='python', index_col=0)
    if not os.path.exists(s["Bo_Curve_path"]):
        os.makedirs(s["Bo_Curve_path"])

    BoCurve_CrystalNmrTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainCrystalNmr_TestEtra.bocurve.csv")
    BoCurve_CrystalTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainCrystal_TestEtra.bocurve.csv")
    BoCurve_NmrTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainNmr_TestEtra.bocurve.csv")

    train_features_del = ["residue_num","residue_name","acc_db","n_homologues","closedist","bind"]
    test_features_del=["residue_num","residue_name","n_homologues","bind","Disruption"]

    database = "ETRA"

    testdata_list = glob.glob(os.path.join(s["thoipapy_feature_folder"],"combined", database, "*.surr20.gaps5.combined_features.csv"))

    dfc = pd.DataFrame()
    dfn = pd.DataFrame()
    dfcn = pd.DataFrame()
    for test_data in testdata_list:
        acc = test_data.split('\\')[-1][0:6]
        test_data_df = pd.read_csv(test_data,sep=',',engine='python',index_col=0)
        odfcn=Train_Input_Test_Pred_Out(acc, crystalnmr_train_data_df, test_data_df, train_features_del, test_features_del, database, s)
        if dfc.empty:
            dfcn = odfcn
        else:
            dfcn = pd.concat([dfc, odfcn], axis=1, join="outer")

        odfc=Train_Input_Test_Pred_Out(acc, crystal_train_data_df, test_data_df, train_features_del, test_features_del, database, s)
        if dfc.empty:
            dfc = odfc
        else:
            dfc = pd.concat([dfc, odfc], axis=1, join="outer")

        odfn=Train_Input_Test_Pred_Out(acc, nmr_train_data_df, test_data_df, train_features_del, test_features_del, database, s)
        if dfn.empty:
            dfn = odfn
        else:
            dfn = pd.concat([dfn, odfn], axis=1, join="outer")

    dfcn.to_csv(BoCurve_CrystalNmrTrain_output_path)
    dfc.to_csv(BoCurve_CrystalTrain_output_path)
    dfn.to_csv(BoCurve_NmrTrain_output_path)




def Test_Crystal(s):
    if not os.path.exists(s["Bo_Curve_path"]):
        os.makedirs(s["Bo_Curve_path"])

    test_data_type_list=["crystal","NMR"]
    list_name_lists =["01","02"]
    for test_data_type, list_name in zip(test_data_type_list, list_name_lists):
        BoCurve_CrystalNmrTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainCrystalNmr_Test{}.bocurve.csv".format(test_data_type))
        BoCurve_CrystalTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainCrystal_Test{}.bocurve.csv".format(test_data_type))
        BoCurve_NmrTrain_output_path = os.path.join(s["Bo_Curve_path"], "TrainNmr_Test{}.bocurve.csv".format(test_data_type))

        train_features_del = ["residue_num","residue_name","acc_db","n_homologues","closedist","bind"]
        test_features_del=["residue_num","residue_name","acc_db","n_homologues","closedist","bind"]

        dfc = pd.DataFrame()
        dfn = pd.DataFrame()
        dfcn = pd.DataFrame()

        set_path = s["set_path"]

        df_set = pd.read_excel(os.path.join(set_path,"set{}_{}.xlsx".format(list_name,test_data_type)), sheetname='proteins')
        for i in df_set.index:
            acc = df_set.loc[i, "acc"]
            database = df_set.loc[i, "database"]
            feature_combined_file = os.path.join(s["thoipapy_feature_folder"], "combined", database, "{}.surr20.gaps5.combined_features.csv".format(acc))

            test_data_df = pd.read_csv(feature_combined_file, index_col=0)
            test_data_df["acc_db"] = "{}-{}".format(acc, database)

            # reorder the columns
            test_data_df = thoipapy.utils.reorder_dataframe_columns(test_data_df, ['acc_db', 'residue_num', 'residue_name', 'n_homologues'])

            set04_path= os.path.join(s["set_path"],"set04_crystal_NMR.xlsx")
            train04_df = thoipapy.figs.fig_utils.create_one_out_train_data(acc, set04_path ,s)
            set01_path= os.path.join(s["set_path"],"set01_crystal.xlsx")
            train01_df = thoipapy.figs.fig_utils.create_one_out_train_data(acc, set01_path ,s)
            set02_path= os.path.join(s["set_path"],"set02_NMR.xlsx")
            train02_df = thoipapy.figs.fig_utils.create_one_out_train_data(acc, set02_path ,s)

            odfcn=Train_Input_Test_Pred_Out(acc, train04_df, test_data_df, train_features_del, test_features_del, database, s)
            if dfc.empty:
                dfcn = odfcn
            else:
                dfcn = pd.concat([dfc, odfcn], axis=1, join="outer")

            odfc=Train_Input_Test_Pred_Out(acc, train01_df, test_data_df, train_features_del, test_features_del, database, s)
            if dfc.empty:
                dfc = odfc
            else:
                dfc = pd.concat([dfc, odfc], axis=1, join="outer")

            odfn=Train_Input_Test_Pred_Out(acc, train02_df, test_data_df, train_features_del, test_features_del, database, s)
            if dfn.empty:
                dfn = odfn
            else:
                dfn = pd.concat([dfn, odfn], axis=1, join="outer")
        dfcn.to_csv(BoCurve_CrystalNmrTrain_output_path)
        dfc.to_csv(BoCurve_CrystalTrain_output_path)
        dfn.to_csv(BoCurve_NmrTrain_output_path)




def Train_Input_Test_Pred_Out(acc, train_data_df, tdf, train_features_del, test_features_del, database, s):

    X=train_data_df.drop(train_features_del,axis=1)
    y=train_data_df["interface"]
    clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    fit = clf.fit(X,y)
    if database == "crystal" or database == "NMR":
        disrupt_or_closedist = tdf.closedist
        disrupt_or_closedist = -1 * disrupt_or_closedist  #(it is closest distance and low value means high propencity of interfacial)
    if database == "ETRA":
        dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
        ddf = pd.read_excel(dirupt_path, index_col=0)
        disrupt_or_closedist = ddf.Disruption
        disrupt_or_closedist = disrupt_or_closedist      #(it is closest experimental disruption and high value means high propencity of interfacial)
    tdf.index = tdf.index.astype(int) +1
    Lips_score = tdf.LIPS_lipo * tdf.LIPS_entropy
    tX=tdf.drop(test_features_del,axis=1)
    if hasattr(clf,'predict_proba'):
        prob_pos=fit.predict_proba(tX)[:,1]
    else:
        prob_pos=fit.decision_function(tX)
        prob_pos=(prob_pos - prob_pos.min())/(prob_pos.max() - prob_pos.min())
    odf=thoipapy.figs.fig_utils.Bo_Curve_Create(acc,prob_pos,Lips_score,disrupt_or_closedist,database)
    return odf