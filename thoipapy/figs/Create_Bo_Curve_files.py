import thoipapy
import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier
import glob
import numpy as np

def Test_Etra(s):
    crystalnmr_train_data = os.path.join(s["Result_folder"], "set04\set04_train_data.csv")
    crystal_train_data = os.path.join(s["Result_folder"], "set01\set01_train_data.csv")
    nmr_train_data = os.path.join(s["Result_folder"], "set02\set02_train_data.csv")
    BoCurvePath = os.path.join(s["Result_folder"], "Bo_Curve")
    if not os.path.exists(BoCurvePath):
        os.makedirs(BoCurvePath)
    BoCurve_CrystalNmrTrain_output_path = os.path.join(BoCurvePath, "TrainCrystalNmr_TestEtra.bocurve.csv")
    BoCurve_CrystalTrain_output_path = os.path.join(BoCurvePath, "TrainCrystal_TestEtra.bocurve.csv")
    BoCurve_NmrTrain_output_path = os.path.join(BoCurvePath, "TrainNmr_TestEtra.bocurve.csv")

    TrainCrystalNmr_OrTrainNmr_TestEtra(crystalnmr_train_data,BoCurve_CrystalNmrTrain_output_path,s)
    TrainCrystalNmr_OrTrainNmr_TestEtra(crystal_train_data, BoCurve_CrystalTrain_output_path, s)
    TrainCrystalNmr_OrTrainNmr_TestEtra(nmr_train_data,BoCurve_NmrTrain_output_path,s)


def TrainCrystalNmr_OrTrainNmr_TestEtra(train_data,bocurve_output,s):
    train_features_del = ["residue_num","residue_name","acc_db","n_homologues","closedist","bind"]
    test_features_del=["residue_num","residue_name","n_homologues"]
    dfm=pd.DataFrame()


    data = pd.read_csv(train_data,sep=',',engine='python',index_col=0)
    X=data.drop(train_features_del,axis=1)
    y=data["bind"]
    clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    fit = clf.fit(X,y)

    database = "ETRA"
    testdata_list = glob.glob(os.path.join(s["thoipapy_feature_folder"],"combined", database, "*.surr20.gaps5.combined_features_incl_phys_param.csv"))
    for test_data in testdata_list:
        acc = test_data.split('\\')[-1][0:6]
        dirupt_path = os.path.join(s["toxr_disruption_folder"],"{}_mul_scan_average_data.xlsx".format(acc))
        ddf = pd.read_excel(dirupt_path,index_col=0)
        disruption = ddf.Disruption
        tdf = pd.read_csv(test_data,sep=',',engine='python',index_col=0)
        tdf.index = tdf.index.astype(int) +1
        Lips_score = tdf.LIPS_lipo * tdf.LIPS_entropy
        tX=tdf.drop(test_features_del,axis=1)
        if hasattr(clf,'predict_proba'):
            prob_pos=fit.predict_proba(tX)[:,1]
        else:
            prob_pos=fit.decision_function(tX)
            prob_pos=(prob_pos - prob_pos.min())/(prob_pos.max() - prob_pos.min())
        odf=thoipapy.figs.fig_utils.Bo_Curve_Create(acc,prob_pos,Lips_score,disruption,database)
        if dfm.empty:
            dfm=odf
        else:
            dfm=pd.concat([dfm,odf],axis=1,join="outer")
    dfm.to_csv(bocurve_output)