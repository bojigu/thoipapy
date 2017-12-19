import ast
import numpy as np
import pandas as pd
import thoipapy
import os
import matplotlib.pyplot as plt
import sys
import pickle
from korbinian.utils import convert_truelike_to_bool, convert_falselike_to_bool

def fig_plot_BO_curve_mult_train_datasets(s):
    """Plot the BO-curve for multiple training datasets.

    Takes the datasets listed in settings under "train_datasets" and "test_datasets"
    and plots the BO-curve of each combination in a single figure.

    The Area Under the BO Curve for a sample size of 0 to 10 (AUBOC10) is shown in the legend.

    Currently plots both the new and old performance method.

    NEW METHOD
    ----------
    Performance = overlap between experiment and predicted MINUS the overlap expected in random selections

    OLD METHOD
    ----------
    Performance = overlap between experiment and predicted DIVIDED BY the overlap expected in random selections

    Parameters
    ----------
    s : dict
        Settings dictionary for figures.

    """

    plt.rcParams.update({'font.size': 7})

    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    test_dataset_str = "-".join([str(n) for n in test_set_list])
    train_dataset_str = "-".join([str(n) for n in train_set_list])

    mult_testname = "testsets({})_trainsets({})".format(test_dataset_str, train_dataset_str)
    print(mult_testname)
    mult_THOIPA_dir = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "summaries", mult_testname)
    thoipapy.utils.make_sure_path_exists(mult_THOIPA_dir)

    plot_BO_curve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname)

    plot_BO_curve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheetname="df_o_over_r", suffix="_BO_curve_old_method")


def plot_BO_curve(s,train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheetname="df_o_minus_r", suffix="_BO_curve"):
    """ Separate function allowing a toggle of the OLD or NEW performance methods

    Parameters
    ----------
	s : dict
        Settings dictionary for figures.
    train_set_list : list
        List of training datasets in selection
        E.g. ["set02", "set04"]
    test_set_list : list
        List of test datasets in selection
        E.g. ["set03", "set31"]
    mult_THOIPA_dir : str
        Path to folder containing results for multiple THOIPA comparisons.
    mult_testname : str
        String denoting this combination of test and training datasets
        E.g. testsets(2)_trainsets(2)
    sheetname : str
        Excel sheetname
        This is the toggle deciding whether the OLD or NEW performance measure is used
        Default = new method ("df_o_minus_r"), where the overlap MINUS random_overlap is used
    suffix : str
        Suffix for figure
        E.g. "" or "_old_method_o_over_r"

    """

    BO_curve_png = os.path.join(mult_THOIPA_dir, "{}{}.png".format(mult_testname, suffix))

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            #/media/mark/sindy/m_data/THOIPA_data/Results/Bo_Curve/Testset03_Trainset01.THOIPA.validation/bo_curve_underlying_data_indiv_df.xlsx
            bo_curve_underlying_data_indiv_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "bo_curve_underlying_data_indiv_df.xlsx")

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


def compare_predictors(s):
    """Plot the BO-curve for multiple prediction methods

    Takes the datasets listed in settings under fig_plot_BO_curve_mult_predictors_list
    (e.g. ["Testset03_Trainset04.THOIPA","Testset03.LIPS"])
    and plots the BO-curves in a single figure.

    The Area Under the BO Curve for a sample size of 0 to 10 (AUBOC10) is shown in the legend.

    Currently plots both the new and old performance method.

    Performance is measured with the NEW METHOD:
    Performance = overlap between experiment and predicted MINUS the overlap expected in random selections

    Parameters
    ----------
    s : dict
        Settings dictionary for figures.

    """

    plt.rcParams.update({'font.size': 7})
    mult_pred_dir = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors")
    BO_curve_png = os.path.join(mult_pred_dir, "BO_curve_mult_pred.png")
    AUBOC10_bar_png = os.path.join(mult_pred_dir, "AUBOC10_barchart_mult_pred.png")
    ROC_png = os.path.join(mult_pred_dir, "ROC.png")

    thoipapy.utils.make_sure_path_exists(mult_pred_dir)

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    # list of predictors to compare, e.g. ["Testset03_Trainset04.THOIPA", "Testset03.LIPS"]
    #predictor_list = ast.literal_eval(s["fig_plot_BO_curve_mult_predictors_list"])

    predictors_df = pd.read_excel(s["excel_file_with_settings"], sheetname="predictors")
    predictors_df["include"] = predictors_df["include"].apply(convert_truelike_to_bool, convert_nontrue=False)
    predictors_df["include"] = predictors_df["include"].apply(convert_falselike_to_bool)
    predictors_df = predictors_df.loc[predictors_df.include == True]
    predictor_list = predictors_df.predictor.tolist()

    area_under_curve_dict = {}
    # create an empty dataframe to keep the pycharm IDE happy
    df = pd.DataFrame()

    for predictor_name in predictor_list:
        bo_curve_underlying_data_indiv_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "{}".format(predictor_name), "bo_curve_underlying_data_indiv_df.xlsx")

        if not os.path.isfile(bo_curve_underlying_data_indiv_xlsx):
            raise FileNotFoundError("bo_curve_underlying_data_indiv_xlsx does not exist ({}). Try running run_testset_trainset_validation in run_figs.py".format(bo_curve_underlying_data_indiv_xlsx))

        df = pd.read_excel(bo_curve_underlying_data_indiv_xlsx, sheetname="df_o_minus_r", index_col=0)

        df["mean_"] = df.mean(axis=1)

        # use the composite trapezoidal rule to get the area under the curve
        # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.trapz.html
        area_under_curve = np.trapz(y=df["mean_"], x=df.index)
        area_under_curve_dict[predictor_name] = area_under_curve

        df["mean_"].plot(ax=ax, label="{}(AUBOC10={:0.1f})".format(predictor_name, area_under_curve))

    ax.set_xlabel("sample size")
    ax.set_ylabel("performance\n(observed overlap - random overlap)")
    ax.set_xticks(range(1, df.shape[0] + 1))
    ax.set_xticklabels(df.index)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BO_curve_png, dpi=240)
    fig.savefig(BO_curve_png[:-4] + ".pdf")

    plt.close("all")
    AUBOC10_ser = pd.Series(area_under_curve_dict).sort_index()
    fig, ax = plt.subplots(figsize=(3.42, 3.42))
    AUBOC10_ser.plot(ax=ax, kind="bar")
    ax.set_ylabel("performance (AUBOC10)")
    fig.tight_layout()
    fig.savefig(AUBOC10_bar_png, dpi=240)
    fig.savefig(AUBOC10_bar_png[:-4] + ".pdf")

    plt.close("all")

    fig, ax = plt.subplots(figsize=(3.42, 3.42))

    for predictor_name in predictor_list:
        #"D:\data_thoipapy\Results\compare_testset_trainset\data\Testset03_Trainset04.THOIPA\Testset03_Trainset04.THOIPA.ROC_data.pkl"
        ROC_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", predictor_name, "{}.ROC_data.pkl".format(predictor_name))

        if os.path.isfile(ROC_pkl):
            with open(ROC_pkl, "rb") as f:
                ROC_out_dict = pickle.load(f)
                ax.plot(ROC_out_dict["false_positive_rate_mean"], ROC_out_dict["true_positive_rate_mean"], label='{} ({:0.2f})'.format(predictor_name, ROC_out_dict["mean_auc"]), lw=1.5)
        else:
            sys.stdout.write("PICKLE WITH ROC DATA NOT FOUND : {}".format(ROC_pkl))

    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_png, dpi=240)
    fig.savefig(ROC_png[:-4] + ".pdf")

    sys.stdout.write("\ncompare_predictors finished ({})".format(BO_curve_png))

def combine_BO_curve_files_HARDLINKED(s):

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
