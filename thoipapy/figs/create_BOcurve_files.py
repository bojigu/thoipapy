import warnings
from pathlib import Path
from typing import Union

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.stats import linregress
from thoipapy.utils import normalise_0_1, make_sure_path_exists

warnings.filterwarnings("ignore")


def save_BO_linegraph_and_barchart(s, bocurve_data_xlsx, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser, plot_o_over_r=False):

    df_o_minus_r = pd.read_excel(bocurve_data_xlsx, sheet_name="df_o_minus_r", index_col=0)
    BO_scatter_png = str(BO_barchart_png)[:-12] + "scatter.png"

    #######################################################################################################
    #                                                                                                     #
    #               Create a dataframe with AUBOC10 and AUC for individual protein (df_valid_indiv)       #
    #                                                                                                     #
    #######################################################################################################
    # load AUBOC10 values as a series
    mean_o_minus_r_by_sample_ser = pd.read_excel(bocurve_data_xlsx, sheet_name="mean_o_minus_r_by_sample", index_col=0)["mean_o_minus_r_by_sample"]
    # select sample sizes 5 and 10
    df_valid_indiv = df_o_minus_r.loc[[5, 10], :].T.copy()
    df_valid_indiv["AUBOC10"] = mean_o_minus_r_by_sample_ser
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

    fig.savefig(BO_scatter_png, dpi=240)

    # simply normalise all between 0 and 1
    for col in df_valid_indiv.columns:
        df_valid_indiv[col] = normalise_0_1(df_valid_indiv[col])[0] + 0.01

    bocurve_data_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/data/{s['setname']}_thoipa_loo_bo_curve_data.xlsx"
    BO_data_valid_indiv_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/data/{s['setname']}_BO_curve_data_valid_indiv.csv"
    make_sure_path_exists(bocurve_data_xlsx, isfile=True)

    df_valid_indiv = df_valid_indiv.reindex(columns=["AUBOC10", 5, 10, "ROC AUC"])
    df_valid_indiv.columns = ["AUBOC10", "sample size 5", "sample size 10", "ROC AUC"]

    df_valid_indiv.to_csv(BO_data_valid_indiv_csv)

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
        df_o_over_r = pd.read_excel(bocurve_data_xlsx, sheet_name="df_o_over_r", index_col=0)
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


def save_extra_BO_figs(bocurve_data_xlsx, other_figs_path):
    linechart_mean_obs_and_rand = os.path.join(other_figs_path, "1_linechart_mean_obs_and_rand.png")
    linechart_obs_indiv = os.path.join(other_figs_path, "2_linechart_obs_indiv.png")
    linechart_p_indiv = os.path.join(other_figs_path, "3_linechart_p_indiv.png")
    linechart_o_minus_r = os.path.join(other_figs_path, "4_linechart_o_minus_r.png")
    linechart_o_over_r = os.path.join(other_figs_path, "5_linechart_o_over_r.png")

    dfrand = pd.read_excel(bocurve_data_xlsx, sheet_name="dfrand", index_col=0)
    dfobs = pd.read_excel(bocurve_data_xlsx, sheet_name="dfobs", index_col=0)
    df_o_minus_r = pd.read_excel(bocurve_data_xlsx, sheet_name="df_o_minus_r", index_col=0)
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

    dfp = pd.read_excel(bocurve_data_xlsx, sheet_name="dfp", index_col=0)
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
