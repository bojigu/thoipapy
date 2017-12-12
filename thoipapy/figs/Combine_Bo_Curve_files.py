import pandas as pd

def combine_BO_curve_files(s):

    Train04_Test01_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset04_Testset01.bocurve.csv"
    df41 = pd.read_csv(Train04_Test01_BoCurve_file,index_col=0)
    df41_ratio = df41[df41.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df41[df41.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df41_ratio_df = df41_ratio.to_frame(name="Tr4Te1Ratio")
    df41_LIPS_ratio = df41[df41.parameters=="LIPS_observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df41[df41.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df41_LIPS_ratio_df = df41_LIPS_ratio.to_frame(name="Tr4Te1LIPSRatio")


    Train04_Test02_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset04_Testset02.bocurve.csv"
    df42 = pd.read_csv(Train04_Test02_BoCurve_file,index_col=0)
    df42_ratio = df42[df42.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df42[df42.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df42_ratio_df = df42_ratio.to_frame(name="Tra4Tes2Ratio")
    df42_LIPS_ratio = df42[df42.parameters=="LIPS_observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df42[df42.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df42_LIPS_ratio_df = df42_LIPS_ratio.to_frame(name="Tr4Te2LIPSRatio")


    Train04_Test03_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset04_Testset03.bocurve.csv"
    df43 = pd.read_csv(Train04_Test03_BoCurve_file,index_col=0)
    df43_ratio = df43[df43.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df43[df43.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df43_ratio_df = df43_ratio.to_frame(name="Tra4Tes3Ratio")
    df43_LIPS_ratio = df43[df43.parameters=="LIPS_observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df43[df43.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df43_LIPS_ratio_df = df43_LIPS_ratio.to_frame(name="Tr4Te3LIPSRatio")


    Train02_Test01_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset02_Testset01.bocurve.csv"
    df21 = pd.read_csv(Train02_Test01_BoCurve_file,index_col=0)
    df21_ratio = df21[df21.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df21[df21.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df21_ratio_df = df21_ratio.to_frame(name="Tra2Te1Ratio")

    Train02_Test02_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset02_Testset02.bocurve.csv"
    df22 = pd.read_csv(Train02_Test02_BoCurve_file,index_col=0)
    df22_ratio = df22[df22.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df22[df22.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df22_ratio_df = df22_ratio.to_frame(name="Tra2Tes2Ratio")

    Train02_Test03_BoCurve_file = r"D:\THOIPA_data\Results\Bo_Curve\Trainset02_Testset03.bocurve.csv"
    df23 = pd.read_csv(Train02_Test03_BoCurve_file,index_col=0)
    df23_ratio = df23[df23.parameters=="observed_overlap"].drop("parameters",axis=1).mean(axis=1)/df23[df23.parameters=="random_overlap"].drop("parameters",axis=1).mean(axis=1)
    df23_ratio_df = df23_ratio.to_frame(name="Tra2Te3Ratio")

    dfc=pd.DataFrame()

    Combined_BoCurve_file= r"D:/THOIPA_data/Results/Bo_Curve/Combined_Bo_Curve_ratio_file.csv"
    dfc = pd.concat([df41_ratio_df, df42_ratio_df,df43_ratio_df,df21_ratio_df,df22_ratio_df,df23_ratio_df,df41_LIPS_ratio_df,df42_LIPS_ratio_df,df43_LIPS_ratio_df], axis=1, join="outer")
    dfc.index = range(1,11)
    dfc.index.name = "sample_size"
    dfc
    dfc.to_csv(Combined_BoCurve_file)
