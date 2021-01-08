import os
import pickle
import sys
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from thoipapy.utils import convert_truelike_to_bool, convert_falselike_to_bool
import thoipapy


def fig_plot_BOcurve_mult_train_datasets(s):
    """Plot the BO-curve for multiple training datasets.

    Takes the datasets listed in settings under "train_datasets" and "test_datasets"
    and plots the BO-curve of each combination in a single figure.

    The Area Under the BO Curve for a sample size of 0 to 10 (AUBOC) is shown in the legend.

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

    # plt.rcParams.update({'font.size': 7})

    test_set_list, train_set_list = thoipapy.utils.get_test_and_train_set_lists(s)

    test_dataset_str = "-".join([str(n) for n in test_set_list])
    train_dataset_str = "-".join([str(n) for n in train_set_list])

    mult_testname = "testsets({})_trainsets({})".format(test_dataset_str, train_dataset_str)
    sys.stdout.write(mult_testname)
    mult_THOIPA_dir = os.path.join(s["data_dir"], "results", "compare_testset_trainset", "summaries", mult_testname)
    thoipapy.utils.make_sure_path_exists(mult_THOIPA_dir)

    plot_BOcurve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname)

    plot_BOcurve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheet_name="df_o_over_r", suffix="_BO_curve_old_method")


def plot_BOcurve(s, train_set_list, test_set_list, mult_THOIPA_dir, mult_testname, sheet_name="df_o_minus_r", suffix="_BO_curve"):
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
    sheet_name : str
        Excel sheet_name
        This is the toggle deciding whether the OLD or NEW performance measure is used
        Default = new method ("df_o_minus_r"), where the overlap MINUS random_overlap is used
    suffix : str
        Suffix for figure
        E.g. "" or "_old_method_o_over_r"

    """

    BO_curve_png = os.path.join(mult_THOIPA_dir, "{}{}.png".format(mult_testname, suffix))

    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for train_set in train_set_list:
        trainsetname = "set{:02d}".format(int(train_set))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))
            # /media/mark/sindy/m_data/THOIPA_data/results/Bo_Curve/Testset03_Trainset01.THOIPA.validation/bocurve_data.xlsx
            bocurve_data_xlsx = os.path.join(s["data_dir"], "results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "data", "bocurve_data.xlsx")

            df = pd.read_excel(bocurve_data_xlsx, sheet_name=sheet_name, index_col=0)

            df["mean_"] = df.mean(axis=1)

            # apply cutoff (e.g. 5 residues for AUBOC5)
            auboc_ser = df["mean_"].iloc[:s["n_residues_AUBOC_validation"]]

            # use the composite trapezoidal rule to get the area under the curve
            # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.trapz.html

            AUBOC = np.trapz(y=auboc_ser, x=auboc_ser.index)

            df["mean_"].plot(ax=ax, label="Test{}_Train{}(AUBOC={:0.1f})".format(testsetname, trainsetname, AUBOC))

    ax.set_xlabel("sample size")
    ax.set_ylabel("performance\n(observed overlap - random overlap)")
    ax.set_xticks(range(1, df.shape[0] + 1))
    ax.set_xticklabels(df.index)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BO_curve_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(BO_curve_png))
    sys.stdout.write("\nfig_plot_BO_curve_mult_train_datasets finished ({})".format(BO_curve_png))


def compare_selected_predictors(s, logging):
    """Plot the BO-curve for multiple prediction methods

    Takes the datasets listed in settings under the "selected_predictors" tab
    (e.g. ["Testset03_Trainset04.THOIPA","Testset03.LIPS"])
    and plots the BO-curves in a single figure.

    The Area Under the BO Curve for a sample size of 0 to 10 (AUBOC) is shown in the legend.

    Currently plots both the new and old performance method.

    Performance is measured with the NEW METHOD:
    Performance = overlap between experiment and predicted MINUS the overlap expected in random selections

    Parameters
    ----------
    s : dict
        Settings dictionary for figures.

    """
    # if s["set_number"] != s["test_datasets"]:
    #    raise Exception("set_number and test_datasets are not identical in settings file. This is recommended for test/train validation.")

    # plt.rcParams.update({'font.size': 7})
    logging.info("\n--------------- starting compare_selected_predictors ---------------\n")
    BO_curve_png: Union[Path, str] = Path(s["data_dir"]) / f"results/{s['setname']}/blindvalidation/compare_selected_predictors_BO_curve.png"
    AUBOC_bar_png: Union[Path, str] = Path(s["data_dir"]) / f"results/{s['setname']}/blindvalidation/compare_selected_predictors_AUBOC_barchart.png"
    ROC_png: Union[Path, str] = Path(s["data_dir"]) / f"results/{s['setname']}/blindvalidation/compare_selected_predictors_ROC.png"

    thoipapy.utils.make_sure_path_exists(BO_curve_png, isfile=True)

    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    predictors_df = pd.read_excel(s["excel_file_with_settings"], sheet_name="selected_predictors")
    predictors_df["include"] = predictors_df["include"].apply(convert_truelike_to_bool, convert_nontrue=False)
    predictors_df["include"] = predictors_df["include"].apply(convert_falselike_to_bool)
    predictors_df = predictors_df.loc[predictors_df.include == True]
    predictor_list = predictors_df.predictor.tolist()

    area_under_curve_dict = {}
    # create an empty dataframe to keep the pycharm IDE happy
    df = pd.DataFrame()

    for predictor_name in predictor_list:
        bocurve_data_xlsx: Union[Path, str] = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/data/{s['setname']}_thoipa_loo_bo_curve_data.xlsx"

        if not os.path.isfile(bocurve_data_xlsx):
            raise FileNotFoundError("bocurve_data_xlsx does not exist ({}). Try running run_testset_trainset_validation".format(bocurve_data_xlsx))

        df = pd.read_excel(bocurve_data_xlsx, sheet_name="df_o_minus_r", index_col=0)

        df["mean_"] = df.mean(axis=1)

        # apply cutoff (e.g. 5 residues for AUBOC5)
        auboc_ser = df["mean_"].iloc[:s["n_residues_AUBOC_validation"]]

        # use the composite trapezoidal rule to get the area under the curve
        # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.trapz.html
        AUBOC = np.trapz(y=auboc_ser, x=auboc_ser.index)

        area_under_curve_dict[predictor_name] = AUBOC

        df["mean_"].plot(ax=ax, label="{}(AUBOC={:0.1f})".format(predictor_name, AUBOC))

    ax.set_xlabel("sample size")
    ax.set_ylabel("performance\n(observed overlap - random overlap)")
    ax.set_xticks(range(1, df.shape[0] + 1))
    ax.set_xticklabels(df.index)
    ax.legend()
    fig.tight_layout()
    fig.savefig(BO_curve_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(BO_curve_png))

    plt.close("all")
    AUBOC_ser = pd.Series(area_under_curve_dict).sort_index()
    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    AUBOC_ser.plot(ax=ax, kind="bar")
    ax.set_ylabel("performance (AUBOC)")
    fig.tight_layout()
    fig.savefig(AUBOC_bar_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(AUBOC_bar_png))

    plt.close("all")

    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for predictor_name in predictor_list:
        # "D:\data_thoipapy\results\compare_testset_trainset\data\Testset03_Trainset04.THOIPA\Testset03_Trainset04.THOIPA.ROC_data.pkl"
        # ROC_pkl = os.path.join(s["data_dir"], "results", "compare_testset_trainset", "data", predictor_name, "data", "{}.ROC_data.pkl".format(predictor_name))
        testsetname = "set{:02d}".format(int(s['test_datasets']))
        ROC_pkl = Path(s["data_dir"]) / "results" / testsetname / f"blindvalidation/{predictor_name}/ROC_data.pkl"

        if os.path.isfile(ROC_pkl):
            with open(ROC_pkl, "rb") as f:
                ROC_out_dict = pickle.load(f)
                ax.plot(ROC_out_dict["false_positive_rate_mean"], ROC_out_dict["true_positive_rate_mean"], label='{} ({:0.2f})'.format(predictor_name, ROC_out_dict["mean_roc_auc"]), lw=1.5)
        else:
            sys.stdout.write("\nPICKLE WITH ROC DATA NOT FOUND : {}".format(ROC_pkl))
            continue

    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(ROC_png))

    sys.stdout.write("\nBO_curve_png ({})\n".format(BO_curve_png))
    logging.info("\n--------------- finished compare_selected_predictors ---------------\n")


def combine_BOcurve_files_hardlinked(s):
    Train04_Test01_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset04_Testset01.bocurve.csv"
    df41 = pd.read_csv(Train04_Test01_BoCurve_file, index_col=0)
    df41_ratio = df41[df41.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df41[df41.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df41_ratio_df = df41_ratio.to_frame(name="Tr4Te1Ratio")
    df41_LIPS_ratio = df41[df41.parameters == "LIPS_observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df41[df41.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df41_LIPS_ratio_df = df41_LIPS_ratio.to_frame(name="Tr4Te1LIPSRatio")

    Train04_Test02_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset04_Testset02.bocurve.csv"
    df42 = pd.read_csv(Train04_Test02_BoCurve_file, index_col=0)
    df42_ratio = df42[df42.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df42[df42.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df42_ratio_df = df42_ratio.to_frame(name="Tra4Tes2Ratio")
    df42_LIPS_ratio = df42[df42.parameters == "LIPS_observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df42[df42.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df42_LIPS_ratio_df = df42_LIPS_ratio.to_frame(name="Tr4Te2LIPSRatio")

    Train04_Test03_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset04_Testset03.bocurve.csv"
    df43 = pd.read_csv(Train04_Test03_BoCurve_file, index_col=0)
    df43_ratio = df43[df43.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df43[df43.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df43_ratio_df = df43_ratio.to_frame(name="Tra4Tes3Ratio")
    df43_LIPS_ratio = df43[df43.parameters == "LIPS_observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df43[df43.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df43_LIPS_ratio_df = df43_LIPS_ratio.to_frame(name="Tr4Te3LIPSRatio")

    Train02_Test01_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset02_Testset01.bocurve.csv"
    df21 = pd.read_csv(Train02_Test01_BoCurve_file, index_col=0)
    df21_ratio = df21[df21.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df21[df21.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df21_ratio_df = df21_ratio.to_frame(name="Tra2Te1Ratio")

    Train02_Test02_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset02_Testset02.bocurve.csv"
    df22 = pd.read_csv(Train02_Test02_BoCurve_file, index_col=0)
    df22_ratio = df22[df22.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df22[df22.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df22_ratio_df = df22_ratio.to_frame(name="Tra2Tes2Ratio")

    Train02_Test03_BoCurve_file = r"D:\THOIPA_data\results\Bo_Curve\Trainset02_Testset03.bocurve.csv"
    df23 = pd.read_csv(Train02_Test03_BoCurve_file, index_col=0)
    df23_ratio = df23[df23.parameters == "observed_overlap"].drop("parameters", axis=1).mean(axis=1) / df23[df23.parameters == "random_overlap"].drop("parameters", axis=1).mean(axis=1)
    df23_ratio_df = df23_ratio.to_frame(name="Tra2Te3Ratio")

    dfc = pd.DataFrame()

    Combined_BoCurve_file = r"D:/THOIPA_data/results/Bo_Curve/Combined_Bo_Curve_ratio_file.csv"
    dfc = pd.concat([df41_ratio_df, df42_ratio_df, df43_ratio_df, df21_ratio_df, df22_ratio_df, df23_ratio_df, df41_LIPS_ratio_df, df42_LIPS_ratio_df, df43_LIPS_ratio_df], axis=1, join="outer")
    dfc.index = range(1, 11)
    dfc.index.name = "sample_size"
    dfc
    dfc.to_csv(Combined_BoCurve_file)
