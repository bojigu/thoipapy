import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp

import thoipapy
from thoipapy.utils import make_sure_path_exists


def validate_multiple_predictors_and_subsets_auboc(s, df_set, logging):
    logging.info("start create_AUBOC_43databases_figs")

    predictors = ["THOIPA_{}_LOO".format(s["set_number"]), "PREDDIMER", "TMDOCK", "LIPS_surface_ranked", "random"]  # "LIPS_L*E",
    subsets = ["crystal", "NMR", "ETRA"]
    for subset in subsets:
        df_o_minus_r_mean_df = pd.DataFrame()
        AUBOC_list = []
        mean_AUBOC_file = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/data/{s['setname']}.{subset}.4predictors_mean_AUBOC.csv"
        mean_AUBOC_barplot_png = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/figs/{s['setname']}.{subset}.4predictors_mean_AUBOC.png"
        BOCURVE_linechart_csv = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/data/{s['setname']}.{subset}.4predictors_BOCURVE_linechart.csv"
        BOCURVE_linechart_png = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/figs/{s['setname']}.{subset}.4predictors_BOCURVE_linechart.png"
        make_sure_path_exists(BOCURVE_linechart_png, isfile=True)

        make_sure_path_exists(mean_AUBOC_file, isfile=True)

        for predictor in predictors:
            bocurve_data_xlsx = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/indiv_validation/bocurve/data/{predictor}/bocurve_data.xlsx"
            df_o_minus_r = pd.read_excel(bocurve_data_xlsx, sheet_name="df_o_minus_r", index_col=0)
            df_o_minus_r = df_o_minus_r.filter(regex=subset, axis=1)
            df_o_minus_r_mean = df_o_minus_r.T.mean()

            # apply cutoff (e.g. 5 residues for AUBOC5)
            auboc_ser = df_o_minus_r_mean.iloc[:s["n_residues_AUBOC_validation"]]

            AUBOC = np.trapz(y=auboc_ser, x=auboc_ser.index)

            AUBOC_list.append(AUBOC)
            df_o_minus_r_mean_df = pd.concat([df_o_minus_r_mean_df, df_o_minus_r_mean], axis=1, join="outer")

        auboc_mean_df = pd.DataFrame.from_records([AUBOC_list], columns=predictors)
        auboc_mean_df.to_csv(mean_AUBOC_file)
        df_o_minus_r_mean_df.columns = predictors
        df_o_minus_r_mean_df.to_csv(BOCURVE_linechart_csv)
        plt.close("all")
        figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        color_list = ["r", "g", "b", "k", "0.5", "0.25"]
        fig, ax = plt.subplots(figsize=figsize)
        for i, column in enumerate(df_o_minus_r_mean_df.columns):
            # df_o_minus_r_mean_df.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC : {:0.2f}".format(AUBOC))
            label_name = "{}(AUBOC:{:.2f})".format(predictors[i], AUBOC_list[i])
            df_o_minus_r_mean_df[column].plot(ax=ax, linestyle="-", label=label_name, color=color_list[i])
        ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="random", alpha=0.5)
        ax.grid(False)
        ax.set_ylabel("performance value\n(observed - random)", color="#0f7d9b")
        ax.tick_params('y', colors="#0f7d9b")

        ax.spines['left'].set_color("#0f7d9b")
        ax.legend()
        fig.tight_layout()
        fig.savefig(BOCURVE_linechart_png, dpi=140)

        plt.close("all")
        # plt.rcParams.update({'font.size': 8})
        # figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        fig, ax = plt.subplots(figsize=figsize)
        # replace the protein names

        auboc_mean_df.plot(kind="bar", ax=ax, alpha=0.7)

        ax.set_ylabel("performance value\n(AUBOC)")
        ax.legend()  # (["sample size = 5", "sample size = 10"])

        fig.tight_layout()
        ax.grid(False)
        fig.savefig(mean_AUBOC_barplot_png, dpi=240)


def validate_multiple_predictors_and_subsets_auc(s, df_set, logging):
    logging.info("start create_AUC_4predictors_3databases_figs")
    predictors = ["THOIPA_{}_LOO".format(s["set_number"]), "PREDDIMER", "TMDOCK", "LIPS_surface_ranked"]  # "LIPS_L*E",
    subsets = ["crystal", "NMR", "ETRA"]
    for subset in subsets:
        mean_roc_auc_list = []
        mean_tpr_list = []
        # outputs
        mean_AUC_file = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/data/{s['setname']}.{subset}.4predictors_mean_AUC.csv"
        mean_AUC_barplot_png = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/figs/{s['setname']}.{subset}.4predictors_mean_AUC.png"
        ROC_curve_csv = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/data/{s['setname']}.{subset}.4predictors_AUC_ROC.csv"
        AUC_ROC_png = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/figs/{s['setname']}.{subset}.4predictors_AUC_ROC.png"

        make_sure_path_exists(mean_AUC_file, isfile=True)

        plt.close("all")
        figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        fig, ax = plt.subplots(figsize=figsize)
        for predictor_name in predictors:
            mean_roc_auc = []
            mean_tpr = 0.0
            big_list_of_tprs = []
            mean_fpr = np.linspace(0, 1, 100)
            n = 0
            # input
            auc_pkl = Path(s["data_dir"]) / f"results/{s['setname']}/crossvalidation/indiv_validation/roc_auc/{predictor_name}/ROC_AUC_data.pkl"
            # open pickle file
            with open(auc_pkl, "rb") as f:
                xv_dict = pickle.load(f)
                for k, v in xv_dict.items():
                    if re.search(subset, k):
                        mean_roc_auc.append(v['roc_auc'])
                        # mean_roc_auc.append(v['auc'])
                        fpr = v['fpr']
                        tpr = v['tpr']
                        mean_tpr += interp(mean_fpr, fpr, tpr)
                        mean_tpr[0] = 0.0
                        n = n + 1

                        big_list_of_tprs.append(tpr)
            mean_tpr /= n
            mean_tpr[-1] = 1.0

            dfb = pd.DataFrame(big_list_of_tprs)

            mean_roc_auc = np.mean(mean_roc_auc)
            mean_roc_auc_list.append(mean_roc_auc)
            mean_tpr_list.append(mean_tpr)
            ax.plot(mean_fpr, mean_tpr, lw=1, label="{} (area = {:.2f})".format(predictor_name, mean_roc_auc), alpha=0.8)
        ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
        ax.set_xlim([-0.05, 1.05])
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(AUC_ROC_png, dpi=240)
        # fig.savefig(thoipapy.utils.pdf_subpath(AUC_ROC_png))

        df_tpr = pd.DataFrame.from_records(list(map(list, zip(*mean_tpr_list))),
                                           columns=predictors)
        df_tpr.to_csv(ROC_curve_csv)

        auc_mean_df = pd.DataFrame.from_records([mean_roc_auc_list], columns=predictors)
        auc_mean_df.to_csv(mean_AUC_file)

        plt.close("all")
        # plt.rcParams.update({'font.size': 8})
        # figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        fig, ax = plt.subplots(figsize=figsize)
        # replace the protein names

        auc_mean_df.plot(kind="bar", ax=ax, alpha=0.7)

        ax.set_ylabel("performance value\n(AUC)")
        ax.legend()  # (["sample size = 5", "sample size = 10"])

        fig.tight_layout()
        ax.grid(False)
        fig.savefig(mean_AUC_barplot_png, dpi=240)
