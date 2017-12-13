import numpy as np
import pandas as pd
import thoipapy
import os
import matplotlib.pyplot as plt
import sys

def fig_plot_BO_curve_mult_train_datasets(s):

    plt.rcParams.update({'font.size': 7})

    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    test_dataset_str = "-".join([str(n) for n in test_set_list])
    train_dataset_str = "-".join([str(n) for n in train_set_list])

    mult_testname = "testsets({})_trainsets({})".format(test_dataset_str, train_dataset_str)
    mult_THOIPA_dir = os.path.join(s["Bo_Curve_path"], "mult_THOIPA", mult_testname)
    thoipapy.utils.make_sure_path_exists(mult_THOIPA_dir)

    plot_BO_curve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname)

    plot_BO_curve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheetname="df_o_over_r", suffix="_old_method_o_over_r")


def plot_BO_curve(s,train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheetname="df_o_minus_r", suffix=""):

    BO_curve_png = os.path.join(mult_THOIPA_dir, "{}{}.png".format(mult_testname, suffix))

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            #/media/mark/sindy/m_data/THOIPA_data/Results/Bo_Curve/Testset03_Trainset01.THOIPA.best_overlap_data/bo_curve_underlying_data_indiv_df.xlsx
            bo_curve_underlying_data_indiv_xlsx = os.path.join(s["Bo_Curve_path"], "Test{}_Train{}.THOIPA.best_overlap_data".format(testsetname, trainsetname), "bo_curve_underlying_data_indiv_df.xlsx")

            df = pd.read_excel(bo_curve_underlying_data_indiv_xlsx, sheetname=sheetname, index_col=0)

            df["mean_"] = df.mean(axis=1)

            # use the composite trapezoidal rule to get the area under the curve
            # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.trapz.html
            area_under_curve = np.trapz(y=df["mean_"], x=df.index)

            df["mean_"].plot(ax=ax, label="Test{}_Train{}(AUBOC10={:0.1f})".format(testsetname, trainsetname, area_under_curve))

    ax.set_xlabel("sample size")
    ax.set_ylabel("performance\n(observed overlap - random overlap)")
    ax.set_xticks(range(1, df.shape[0] + 1))
    ax.set_xticklabels(df.index)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BO_curve_png, dpi=240)
    fig.savefig(BO_curve_png[:-4] + ".pdf")
    sys.stdout.write("\nfig_plot_BO_curve_mult_train_datasets finished ({})".format(BO_curve_png))



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
