import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp

import thoipapy


def create_AUBOC10_4predictors_3databases_figs(s,df_set,logging):

    logging.info("start create_AUBOC10_43databases_figs")
    predictor_name_list = ["THOIPA_{}_LOO".format(s["set_number"]),"PREDDIMER", "TMDOCK", "LIPS_surface_ranked"] #"LIPS_L*E",
    databases = ["crystal", "NMR", "ETRA"]
    for database in databases:
        df_o_minus_r_mean_df = pd.DataFrame()
        AUBOC10_list = []
        mean_AUBOC_file = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_mean_AUBOC10.csv"
        mean_AUBOC_barplot_png = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_mean_AUBOC10.png"
        BOCURVE_linechart_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_BOCURVE_linechart.csv"
        BOCURVE_linechart_png = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_BOCURVE_linechart.png"

        for predictor_name in predictor_name_list:
            crossvalidation_dir = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "crossvalidation")
            BO_data_excel = os.path.join(crossvalidation_dir, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))
            df_o_minus_r = pd.read_excel(BO_data_excel,sheet_name="df_o_minus_r",index_col=0)
            df_o_minus_r = df_o_minus_r.filter(regex=database, axis=1)
            df_o_minus_r_mean = df_o_minus_r.T.mean()
            AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)
            AUBOC10_list.append(AUBOC10)
            df_o_minus_r_mean_df = pd.concat([df_o_minus_r_mean_df, df_o_minus_r_mean], axis=1, join="outer")

        auboc_mean_df = pd.DataFrame.from_records([AUBOC10_list], columns=predictor_name_list)
        auboc_mean_df.to_csv(mean_AUBOC_file)
        df_o_minus_r_mean_df.columns = predictor_name_list
        df_o_minus_r_mean_df.to_csv(BOCURVE_linechart_csv)
        plt.close("all")
        figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        color_list = 'rgbk'
        fig, ax = plt.subplots(figsize=figsize)
        for i,column in enumerate(df_o_minus_r_mean_df.columns):
            # df_o_minus_r_mean_df.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC10 : {:0.2f}".format(AUBOC10))
            label_name = "{}(AUBOC:{:.2f})".format(predictor_name_list[i] ,AUBOC10_list[i])
            df_o_minus_r_mean_df[column].plot(ax=ax,  linestyle="-",label=label_name, color = color_list[i])
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

        ax.set_ylabel("performance value\n(AUBOC10)")
        ax.legend()  # (["sample size = 5", "sample size = 10"])

        fig.tight_layout()
        ax.grid(False)
        fig.savefig(mean_AUBOC_barplot_png, dpi=240)


def create_AUC_4predictors_3databases_figs(s,df_set,logging):

    logging.info("start create_AUC_4predictors_3databases_figs")
    predictor_name_list = ["THOIPA_{}_LOO".format(s["set_number"]),"PREDDIMER", "TMDOCK", "LIPS_surface_ranked"] #"LIPS_L*E",
    databases = ["crystal", "NMR", "ETRA"]
    for database in databases:
        mean_roc_auc_list = []
        mean_tpr_list= []

        mean_AUC_file = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_mean_AUC.csv"
        mean_AUC_barplot_png = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_mean_AUC.png"
        ROC_curve_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_AUC_ROC.csv"
        AUC_ROC_png = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{database}.4predictors_AUC_ROC.png"

        plt.close("all")
        figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        fig, ax = plt.subplots(figsize=figsize)
        for predictor_name in predictor_name_list:
            mean_roc_auc = []
            mean_tpr = 0.0
            big_list_of_tprs = []
            mean_fpr = np.linspace(0, 1, 100)
            n=0
            #AUC_data_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
            #                             "{}.{}_AUC_data.pkl".format(s["setname"],predictor_name.replace('*',"")))
            AUC_data_pkl = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/indiv_validation/{predictor_name}/ROC_AUC_data.pkl"
            # open pickle file
            with open(AUC_data_pkl, "rb") as f:
                xv_dict = pickle.load(f)
                for k,v in xv_dict.items():
                    if re.search(database,k):
                        mean_roc_auc.append(v['roc_auc'])
                        #mean_roc_auc.append(v['auc'])
                        fpr = v['fpr']
                        tpr = v['tpr']
                        mean_tpr += interp(mean_fpr, fpr, tpr)
                        mean_tpr[0] = 0.0
                        n = n+1

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
        # fig.savefig(crossvalidation_png[:-4] + ".pdf")
        fig.savefig(thoipapy.utils.pdf_subpath(AUC_ROC_png))

        df_tpr = pd.DataFrame.from_records(list(map(list, zip(*mean_tpr_list))),
                                                    columns=predictor_name_list)
        df_tpr.to_csv(ROC_curve_csv)

        auc_mean_df = pd.DataFrame.from_records([mean_roc_auc_list], columns=predictor_name_list)
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