import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import pairwise2
from scipy import interp
from scipy.stats import linregress
from sklearn.metrics import roc_curve, auc, precision_recall_curve

import thoipapy


def collect_indiv_validation_data(s, df_set, logging, namedict, predictor_name_list, THOIPA_predictor_name, unique_database_labels):
    """

    Parameters
    ----------
    s
    df_set
    logging
    namedict
    predictor_name_list
    THOIPA_predictor_name

    Returns
    -------

    """
    logging.info("start collect_indiv_validation_data THOIPA_PREDDIMER_TMDOCK")
    ROC_AUC_df = pd.DataFrame()
    PR_AUC_df = pd.DataFrame()
    AUBOC10_df = pd.DataFrame()
    AUBOC10_from_complete_data_ser = pd.Series()

    AUC_AUBOC_name_list = []
    linechar_name_list = []
    AUBOC10_list = []
    df_o_minus_r_mean_df = pd.DataFrame()
    roc_auc_mean_list=[]
    roc_auc_std_list = []
    
    indiv_validation_dir = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "indiv_validation")
    indiv_validation_data_xlsx = os.path.join(indiv_validation_dir, "indiv_validation_data.xlsx")

    thoipapy.utils.make_sure_path_exists(indiv_validation_dir)
    #if not os.path.isdir(os.path.dirname(BOAUC10_barchart_pdf)):
    #    os.makedirs(os.path.dirname(BOAUC10_barchart_pdf))

    for predictor_name in predictor_name_list:
        BO_data_df = pd.DataFrame()
        
        xv_dict = {}
        ROC_AUC_dict = {}
        PR_AUC_dict = {}
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        predictor_dir = os.path.join(indiv_validation_dir, predictor_name)
        auc_pkl = os.path.join(predictor_dir, "ROC_AUC_data.pkl")
        BO_curve_data_csv = os.path.join(predictor_dir, "BO_Curve_data.csv")
        BO_data_excel = os.path.join(predictor_dir, "BO_curve_data.xlsx")
        BO_linechart_png = os.path.join(predictor_dir, "BO_linechart.png")
        BO_barchart_png = os.path.join(predictor_dir, "AUBOC10_barchart.png")
        df_o_minus_r_mean_csv = os.path.join(predictor_dir, "df_o_minus_r_mean.csv")
        thoipapy.utils.make_sure_path_exists(predictor_dir)

        for i in df_set.index:
            sys.stdout.write(".")
            sys.stdout.flush()
            acc = df_set.loc[i, "acc"]

            database = df_set.loc[i, "database"]
            acc_db = acc + "-" + database
            merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))
            merged_data_df = pd.read_csv(merged_data_csv_path,engine="python")
            merged_data_df["LIPS_L*E"] = -1 * merged_data_df["LIPS_L*E"]
            merged_data_df["PREDDIMER"] = -1 * merged_data_df["PREDDIMER"]
            merged_data_df["TMDOCK"] = -1 * merged_data_df["TMDOCK"]

            if database == "crystal" or database == "NMR":
                # (it is closest distance and low value means high propensity of interfacial)
                merged_data_df["interface_score"] = -1 * merged_data_df["interface_score"]
            experiment_col = "interface_score"
            BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, merged_data_df, experiment_col, predictor_name)
            if BO_data_df.empty:
                BO_data_df = BO_single_prot_df
            else:
                BO_data_df = pd.concat([BO_data_df, BO_single_prot_df], axis=1, join="outer")

            df_for_roc = merged_data_df.dropna(subset=["interface", predictor_name])
            fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name], drop_intermediate=False)

            precision, recall, thresholds_PRC = precision_recall_curve(df_for_roc.interface, df_for_roc[predictor_name])

            pr_auc = auc(recall, precision)
            PR_AUC_dict[acc_db] = pr_auc

            roc_auc = auc(fpr, tpr)
            ROC_AUC_dict[acc_db] = roc_auc

            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0

            xv_dict[acc_db] = {"fpr": fpr, "tpr": tpr, "roc_auc": roc_auc, "precision" : precision, "recall" : recall, "pr_auc" : pr_auc}

        # save dict as pickle
        with open(auc_pkl, "wb") as f:
            pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        BO_data_df.to_csv(BO_curve_data_csv)
        # parse BO data csv
        # print out mean values
        thoipapy.figs.create_BOcurve_files.parse_BO_data_csv_to_excel(BO_curve_data_csv, BO_data_excel, logging, predictor_name)

        # ROC AUC validation
        ROC_AUC_ser = pd.Series(ROC_AUC_dict)
        ROC_AUC_ser.sort_values(inplace=True, ascending=False)
        roc_auc_mean_list.append(ROC_AUC_ser.mean())
        roc_auc_std_list.append(ROC_AUC_ser.std())

        # precision-recall AUC validation
        PR_AUC_ser = pd.Series(PR_AUC_dict)
        PR_AUC_ser.sort_values(inplace=True, ascending=False)

        # BO curve AUBOC10 validation
        AUBOC10_ser = pd.read_excel(BO_data_excel, sheetname="AUBOC10", index_col=0)["AUBOC10"].copy()
        df_o_minus_r = pd.read_excel(BO_data_excel, sheetname="df_o_minus_r", index_col=0)
        df_o_minus_r.columns = pd.Series(df_o_minus_r.columns).replace(namedict)
        df_o_minus_r_mean = df_o_minus_r.T.mean()
        #df_o_minus_r_mean_df= pd.concat([df_o_minus_r_mean_df,df_o_minus_r_mean],axis=1, join="outer")
        df_o_minus_r_mean_df[predictor_name] = df_o_minus_r_mean
        AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)
        AUBOC10_list.append(AUBOC10)
        AUBOC10_from_complete_data_ser[predictor_name] = AUBOC10
        linechar_name_list.append(predictor_name)
        AUC_AUBOC_name_list.append("{}-AUC".format(predictor_name))
        AUC_AUBOC_name_list.append("{}-AUBOC10".format(predictor_name))
        thoipapy.figs.create_BOcurve_files.save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict,
                                                                          logging, ROC_AUC_ser)

        ROC_AUC_df[predictor_name] = ROC_AUC_ser
        PR_AUC_df[predictor_name] = PR_AUC_ser
        AUBOC10_df[predictor_name] = AUBOC10_ser

    means_df = pd.DataFrame()
    means_df["ROC_AUC"] = ROC_AUC_df.mean()
    means_df["PR_AUC"] = PR_AUC_df.mean()
    means_df["AUBOC10_mean_indiv"] = AUBOC10_df.mean()
    means_df["AUBOC10_from_complete_data"] = AUBOC10_from_complete_data_ser

    """ means_df looks like this:
                          ROC_AUC    PR_AUC   AUBOC10
    THOIPA_5_LOO         0.629557  0.505823  1.202355
    PREDDIMER            0.566582  0.416761  0.515193
    TMDOCK               0.598387  0.421462  0.666720
    """

    std_df = pd.DataFrame()
    std_df["ROC_AUC"] = ROC_AUC_df.std()
    std_df["PR_AUC"] = PR_AUC_df.std()
    std_df["AUBOC10"] = AUBOC10_df.std()

    SEM_df = pd.DataFrame()
    SEM_df["ROC_AUC"] = ROC_AUC_df.std() / np.sqrt(ROC_AUC_df.shape[0])
    SEM_df["PR_AUC"] = PR_AUC_df.std() / np.sqrt(PR_AUC_df.shape[0])
    SEM_df["AUBOC10"] = AUBOC10_df.std() / np.sqrt(AUBOC10_df.shape[0])

    with pd.ExcelWriter(indiv_validation_data_xlsx) as writer:

        means_df.to_excel(writer, sheet_name="means")
        std_df.to_excel(writer, sheet_name="std")
        SEM_df.to_excel(writer, sheet_name="SEM")

        ROC_AUC_df.to_excel(writer, sheet_name="ROC_AUC_indiv")
        PR_AUC_df.to_excel(writer, sheet_name="PR_AUC_indiv")
        AUBOC10_df.to_excel(writer, sheet_name="BO_AUBOC10_indiv")

        df_o_minus_r_mean_df.to_excel(writer, sheet_name="BO_o_minus_r")

        if "TMDOCK" in PR_AUC_df.columns and "PREDDIMER" in PR_AUC_df.columns:
            df_THOIPA_vs_others = pd.DataFrame()
            df_THOIPA_vs_others["THOIPA_better_TMDOCK"] = PR_AUC_df[THOIPA_predictor_name] > PR_AUC_df.TMDOCK
            df_THOIPA_vs_others["THOIPA_better_PREDDIMER"] = PR_AUC_df[THOIPA_predictor_name] > PR_AUC_df.PREDDIMER
            df_THOIPA_vs_others["THOIPA_better_both"] = df_THOIPA_vs_others[["THOIPA_better_TMDOCK", "THOIPA_better_PREDDIMER"]].sum(axis=1) == 2
            n_THOIPA_better_both = df_THOIPA_vs_others["THOIPA_better_both"].sum()
            logging.info("THOIPA has higher precision-recall AUC than both TMDOCK and PREDDIMER for {}/{} proteins in {}".format(n_THOIPA_better_both, PR_AUC_df.shape[0], s["setname"]))
            df_THOIPA_vs_others.to_excel(writer, sheet_name="THOIPA_vs_others")

        # #sys.stdout.write(roc_auc_mean_list)
        # AUBOC10_mean_df = pd.DataFrame.from_records([AUBOC10_list], columns=linechar_name_list)
        # #AUBOC10_mean_df.to_csv(mean_AUBOC_file)
        # AUBOC10_mean_df.to_excel(writer, sheet_name="AUBOC10_mean")
        # df_o_minus_r_mean_df.columns = linechar_name_list
        # #ROC_AUC_df.columns = AUC_AUBOC_name_list
        # ROC_AUC_df.index.name = "acc_db"
        # #ROC_AUC_df.to_csv(AUC_AUBOC_file)
        # THOIPA_best_set = s["THOIPA_best_set"]
        #
        # # AUC for barchart, 4 predictors, mean AUC of all proteins in dataset
        # #logging.info("_finder : {}".format(mean_roc_auc_barchart_csv))
        # AUC_4pred_mean_all_indiv_prot_df = pd.DataFrame(index = linechar_name_list)
        # #AUC_4pred_mean_all_indiv_prot_df = pd.DataFrame([roc_auc_mean_list, roc_auc_std_list], index = linechar_name_list, columns=["mean", "std"])
        # AUC_4pred_mean_all_indiv_prot_df["roc_auc_mean"] = roc_auc_mean_list
        # AUC_4pred_mean_all_indiv_prot_df["roc_auc_std"] = roc_auc_std_list
        # AUC_4pred_mean_all_indiv_prot_df["n"] = df_set.shape[0]
        # AUC_4pred_mean_all_indiv_prot_df["SEM"] = AUC_4pred_mean_all_indiv_prot_df.roc_auc_std / AUC_4pred_mean_all_indiv_prot_df["n"].apply(np.sqrt)
        # #AUC_4pred_mean_all_indiv_prot_df.to_csv(mean_roc_auc_barchart_csv)
        # AUC_4pred_mean_all_indiv_prot_df.to_excel(writer, sheet_name="ROC_AUC_mean_indiv")


def create_indiv_validation_figs(s, logging, namedict, predictor_name_list, THOIPA_predictor_name, unique_database_labels):
    indiv_validation_dir = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "indiv_validation")
    indiv_validation_data_xlsx = os.path.join(indiv_validation_dir, "indiv_validation_data.xlsx")
    indiv_ROC_AUC_barchart_png = os.path.join(indiv_validation_dir, "indiv_ROC_AUC_barchart.png")
    indiv_PR_AUC_barchart_png = os.path.join(indiv_validation_dir, "indiv_PR_AUC_barchart.png")
    AUBOC10_barchart_png = os.path.join(indiv_validation_dir, "indiv_AUBOC10_barchart.png")
    BOCURVE_linechart_png = os.path.join(indiv_validation_dir, "BOcurve_linechart.png")
    mean_ROC_AUC_barchart_png = os.path.join(indiv_validation_dir, "mean_ROC_AUC_barchart.png")
    mean_PR_AUC_barchart_png = os.path.join(indiv_validation_dir, "mean_PR_AUC_barchart.png")
    ROC_AUC_vs_PR_AUC_scatter_png = os.path.join(indiv_validation_dir, "ROC_AUC_vs_PR_AUC_scatter.png")
    perc_interf_vs_PR_cutoff_linechart_png = os.path.join(indiv_validation_dir, "perc_interf_vs_PR_cutoff_linechart.png")

    ROC_AUC_df = pd.read_excel(indiv_validation_data_xlsx, sheetname="ROC_AUC_indiv")
    PR_AUC_df = pd.read_excel(indiv_validation_data_xlsx, sheetname="PR_AUC_indiv")
    AUBOC10_df = pd.read_excel(indiv_validation_data_xlsx, sheetname="BO_AUBOC10_indiv")
    df_o_minus_r_mean_df = pd.read_excel(indiv_validation_data_xlsx, sheetname="BO_o_minus_r")

    create_ROC_AUC_barchart(ROC_AUC_df, indiv_ROC_AUC_barchart_png, namedict, THOIPA_predictor_name)
    create_PR_AUC_barchart(PR_AUC_df, indiv_PR_AUC_barchart_png, namedict, THOIPA_predictor_name)
    create_AUBOC10_barchart(AUBOC10_df, AUBOC10_barchart_png, namedict, THOIPA_predictor_name)

    create_BOcurve_linechart(df_o_minus_r_mean_df, BOCURVE_linechart_png)

    create_mean_ROC_AUC_barchart(ROC_AUC_df, mean_ROC_AUC_barchart_png)
    create_mean_PR_AUC_barchart(PR_AUC_df, mean_PR_AUC_barchart_png)

    create_scatter_ROC_AUC_vs_PR_AUC(s, predictor_name_list, ROC_AUC_vs_PR_AUC_scatter_png)

    # for the complete list of proteins
    create_linechart_perc_interf_vs_PR_cutoff(s, predictor_name_list, perc_interf_vs_PR_cutoff_linechart_png)

    # for each dataset(e.g. ETRA) separately. Saved in "each_dataset" subfolder
    for database in unique_database_labels:
        perc_interf_vs_PR_cutoff_linechart_single_database_png = os.path.join(indiv_validation_dir, "each_dataset", "{}_perc_interf_vs_PR_cutoff_linechart.png".format(database))
        create_linechart_perc_interf_vs_PR_cutoff(s, predictor_name_list, perc_interf_vs_PR_cutoff_linechart_single_database_png, database=database)

    logging.info("finished run_indiv_validation_THOIPA_PREDDIMER_TMDOCK")

def precision_recall_curve_rises_above_threshold(precision, recall, threshold=0.5):
    """Determines whether the PR curve rises above a threshold P and R value at any point.

    BOTH precision and recall need to simultaneously be above the threshold at some stage, for
    the result to be True.

    Validation inspired by Lensink, M. F., & Wodak, S. J. (2010). Blind predictions of protein interfaces
    by docking calculations in CAPRI. Proteins: Structure, Function and Bioinformatics, 78(15), 3085–3095.
    https://doi.org/10.1002/prot.22850

    Parameters
    ----------
    precision : np.ndarray
        array of precision values in precision-recall curve
    recall : np.ndarray
        array of recall values in precision-recall curve
    threshold : float
        minimum value that has to be attained by BOTH precision and recall

    Returns
    -------
    PR_ rises_above_threshold : boolean
        whether the thresholds were exceeded at any stage in the curve

    Usage
    -----
    precision, recall, thresholds_PRC = precision_recall_curve(true_interface, prediction)
    PR_rises_above_threshold = precision_recall_curve_rises_above_threshold(precision, recall)
    """
    df_pr = pd.DataFrame()
    df_pr["precision"] = precision
    df_pr["recall"] = recall
    df_pr = df_pr > threshold
    df_pr = df_pr.all(axis=1)
    if True in df_pr.tolist():
        PR_rises_above_threshold = True
    else:
        PR_rises_above_threshold = False
    return PR_rises_above_threshold


def create_scatter_ROC_AUC_vs_PR_AUC(s, predictor_name_list, ROC_AUC_vs_PR_AUC_scatter_png):

    fig, ax = plt.subplots(figsize=(8, 8))
    for predictor_name in predictor_name_list:
        auc_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "indiv_validation", predictor_name, "ROC_AUC_data.pkl")
        with open(auc_pkl, "rb") as f:
            xv_dict = pickle.load(f)
        roc_auc_list = []
        pr_auc_list = []
        for key in xv_dict:
            roc_auc_list.append(xv_dict[key]["roc_auc"])
            pr_auc_list.append(xv_dict[key]["pr_auc"])
        ax.scatter(roc_auc_list, pr_auc_list, alpha=0.75, label=predictor_name)
        lr = linregress(roc_auc_list, pr_auc_list)
        yfit = np.array(roc_auc_list) * lr[0] + lr[1]
        ax.plot(roc_auc_list, yfit)
    ax.legend()
    ax.set_xlabel("ROC_AUC")
    ax.set_ylabel("PR_AUC")
    fig.savefig(ROC_AUC_vs_PR_AUC_scatter_png)


def create_mean_ROC_AUC_barchart(ROC_AUC_df, mean_ROC_AUC_barchart_png):
    # plt.close("all")
    # # plt.rcParams.update({'font.size': 2})
    # mean_roc_auc_name = [linechar_name_list[0],'\n{}\n'.format(linechar_name_list[1]),linechar_name_list[2],'\n{}\n'.format(linechar_name_list[3])]
    # figsize = np.array([3.42, 3.42])   # DOUBLE the real size, due to problems on Bo computer with fontsizes
    # fig, ax = plt.subplots(figsize=figsize)
    # # replace the protein names
    # x = y_pos = np.arange(len(linechar_name_list))
    # plt.bar(x, roc_auc_mean_list, width=0.6, color = 'rgbk', alpha=0.5)
    # plt.xticks(y_pos, mean_roc_auc_name,fontsize=6)
    # plt.ylabel("performance value\n(mean auc)")
    # 
    # #ax.set_ylabel("performance value\n(auc)")
    # ax.set_ylim(0, 0.70)
    # ax.legend()  # (["sample size = 5", "sample size = 10"])
    # 
    # fig.tight_layout()
    # ax.grid(False)
    # fig.savefig(mean_ROC_AUC_barchart_png, dpi=240)

    figsize = np.array([3.42, 3.42])   # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    ROC_AUC_df.mean().plot(kind="bar", ax=ax)
    ax.set_ylabel("performance value\n(mean auc)")

    #ax.set_ylabel("performance value\n(auc)")
    ax.set_ylim(0, 0.70)
    ax.legend()  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(mean_ROC_AUC_barchart_png, dpi=240)
    
def create_mean_PR_AUC_barchart(PR_AUC_df, mean_PR_AUC_barchart_png):


    figsize = np.array([3.42, 3.42])   # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    PR_AUC_df.mean().plot(kind="bar", ax=ax)
    ax.set_ylabel("performance value\n(mean auc)")

    #ax.set_ylabel("performance value\n(auc)")
    ax.set_ylim(0, 0.70)
    ax.legend()  # (["sample size = 5", "sample size = 10"])
    ax.set_facecolor('white')

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(mean_PR_AUC_barchart_png, dpi=240)






def create_BOcurve_linechart(df_o_minus_r_mean_df, BOCURVE_linechart_png):
    # BO_linechart_png
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    color_list = 'rgbk'
    fig, ax = plt.subplots(figsize=figsize)
    # for i,column in enumerate(df_o_minus_r_mean_df.columns):
    #     # df_o_minus_r_mean_df.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC10 : {:0.2f}".format(AUBOC10))
    #     label_name = "{}(AUBOC:{:.2f})".format(linechar_name_list[i] ,AUBOC10_list[i])
    #     df_o_minus_r_mean_df[column].plot(ax=ax,  linestyle="-",label=label_name, color = color_list[i])

    df_o_minus_r_mean_df.plot(ax=ax)
    
    ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="random", alpha=0.5)
    ax.grid(False)
    ax.set_ylabel("fraction of correctly predicted residues\n(observed - random)", color="#0f7d9b")
    ax.set_xlabel("number of TMD residues\n(sample size)")
    ax.tick_params('y', colors="#0f7d9b")

    ax.spines['left'].set_color("#0f7d9b")
    ax.legend()
    fig.tight_layout()
    fig.savefig(BOCURVE_linechart_png, dpi=140)

def create_ROC_AUC_barchart(ROC_AUC_df, ROC_AUC_barchart_png, namedict, THOIPA_predictor_name):

    # plt.rcParams.update({'font.size': 8})
    #figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes

    fig, ax = plt.subplots(figsize=figsize)

    # replace the protein names based on the naming file
    ROC_AUC_df.index = pd.Series(ROC_AUC_df.index).replace(namedict)
    ROC_AUC_df.sort_values([THOIPA_predictor_name], ascending=False, inplace=True)
    ROC_AUC_df.plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance (ROC AUC)")
    ax.legend(loc="upper right")  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(ROC_AUC_barchart_png, dpi=240)
    #fig.savefig(ROC_AUC_barchart_png[:-4] + ".pdf")


def create_PR_AUC_barchart(PR_AUC_df, indiv_PR_AUC_barchart_png, namedict, THOIPA_predictor_name):

    # replace the protein names based on the naming file
    PR_AUC_df.index = pd.Series(PR_AUC_df.index).replace(namedict)

    # # plt.rcParams.update({'font.size': 8})
    # # figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    # figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    # fig, ax = plt.subplots(figsize=figsize)
    # PR_AUC_df.sort_values([THOIPA_predictor_name], ascending=False, inplace=True)
    # PR_AUC_df.plot(kind="bar", ax=ax, alpha=0.7)
    # ax.set_ylabel("performance (precision recall AUC)")
    # ax.legend(loc="upper right")  # (["sample size = 5", "sample size = 10"])
    # fig.tight_layout()
    # ax.grid(False)
    # fig.savefig(indiv_PR_AUC_barchart_png, dpi=240)


    Width = 0.2
    Fontsize = 14
    #file_location = r"I:\THOIPA_data\Results\set05\indiv_validation\indiv_validation_data.xlsx"
    #PR_AUC_df = pd.read_excel(file_location, sheet_name="PR_AUC_indiv")

    for i in PR_AUC_df.index:
        if "crystal" in i:
            PR_AUC_df.loc[i, "sort_list"] = 2
        if "NMR" in i:
            PR_AUC_df.loc[i, "sort_list"] = 1
        if "ETRA" in i:
            PR_AUC_df.loc[i, "sort_list"] = 0

    PR_AUC_df.sort_values(['sort_list', THOIPA_predictor_name], ascending=[True, False], inplace=True)

    PR_AUC_df.rename(columns=lambda x: x.replace(THOIPA_predictor_name, 'THOIPA'), inplace=True)
    PR_AUC_df.rename(columns=lambda x: x.replace('LIPS_surface_ranked', 'LIPS'), inplace=True)

    PR_AUC_df.drop(['LIPS', "sort_list"], axis=1, inplace=True)

    color_list = ["#E95D12", "#0065BD", "k", "#B4B3B3"]

    fig, ax = plt.subplots(figsize=(13, 7))
    x = list(range(0, len(PR_AUC_df), 1))

    m = 0
    for i, color in zip(PR_AUC_df.columns, color_list):
        y = PR_AUC_df[i].tolist()
        plt.bar(x, y, width=Width, label=i, color=color)
        x = [i + Width for i in x]

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, ncol=31, loc=2, fontsize=Fontsize, frameon=True, bbox_to_anchor=(-0.0090, 1.15), facecolor='white', edgecolor="k")

    plt.rcParams['xtick.labelsize'] = Fontsize
    plt.rcParams['ytick.labelsize'] = Fontsize
    ax.tick_params(axis='y', labelsize=Fontsize, pad=2)
    ax.tick_params(axis='x', labelsize=Fontsize, pad=2)
    x_label = list(range(0, len(PR_AUC_df), 1))
    x_label = [i + Width for i in x_label]

    ax.set_ylabel('performance (precision-recall AUC)', fontsize=Fontsize, labelpad=1)

    ax.set_xticks(x_label)

    ax.set_xticklabels(PR_AUC_df.index.tolist(), fontsize=Fontsize, rotation=90)
    plt.ylim(0, 1.1)
    fig.tight_layout()
    plt.xlim(-0.5, 54)

    plt.savefig(indiv_PR_AUC_barchart_png, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.savefig(indiv_PR_AUC_barchart_png[:-4] + ".pdf", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()



def create_AUBOC10_barchart(AUBOC10_df, AUBOC10_barchart_png, namedict, THOIPA_predictor_name):
    #AUC_AUBOC_df = AUBOC10_df.T.sort_values(by=[THOIPA_predictor_name], ascending=False)
    # plt.rcParams.update({'font.size': 8})
    # figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes

    fig, ax = plt.subplots(figsize=figsize)
    ## replace the protein names
    #AUBOC10_df.index = pd.Series(AUBOC10_df.index).replace(namedict)
    # replace old "crystal" references with "X-ray"
    AUBOC10_df.index = pd.Series(AUBOC10_df.index).replace("crystal", "X-ray")
    AUBOC10_df.sort_values([THOIPA_predictor_name], ascending=False, inplace=True)
    AUBOC10_df.plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance (AUBOC10))")
    ax.legend()  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(AUBOC10_barchart_png, dpi=240)
    #fig.savefig(AUBOC10_barchart_png[:-4] + ".pdf")


def create_ROC_AUC_barchart_DEPRECATED(ROC_AUC_df, ROC_AUC_barchart_png, ROC_AUC_barchart_pdf, BOAUC10_barchart_png,BOAUC10_barchart_pdf,namedict, THOIPA_best_set):
    bo_auc_list = "?"
    THOIPA_best_setnumber = int(THOIPA_best_set[3:])
    #colname = "THOIPA_5_LOOAUC"
    THOIPA_x_LOOAUC = "THOIPA_{}_LOO-AUC".format(THOIPA_best_setnumber)
    THOIPA_x_LOOAUBOC10 = "THOIPA_{}_LOO-AUBOC10".format(THOIPA_best_setnumber)
    # auc_list = AUC_AUBOC_df.columns[[0, 2, 4, 6]]
    # bo_auc_list = AUC_AUBOC_df.columns[[1, 3, 5, 7]]

    auc_list = ROC_AUC_df.columns[::2]
    #bo_auc_list = ROC_AUC_df.columns[1::2]
    AUC_AUBOC_df = ROC_AUC_df.sort_values(by=[THOIPA_x_LOOAUC], ascending=False)
    plt.close("all")
    # plt.rcParams.update({'font.size': 8})
    #figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    figsize = np.array([9, 6])  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    # replace the protein names

    AUC_AUBOC_df.index = pd.Series(AUC_AUBOC_df.index).replace(namedict)
    AUC_AUBOC_df[auc_list].plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance value\n(auc)")
    ax.legend()  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(ROC_AUC_barchart_png, dpi=240)
    fig.savefig(ROC_AUC_barchart_pdf, dpi=240)

    AUC_AUBOC_df = AUC_AUBOC_df.sort_values(by=[THOIPA_x_LOOAUBOC10], ascending=False)
    plt.close("all")
    # plt.rcParams.update({'font.size': 8})
    #figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    # replace the protein names

    AUC_AUBOC_df.index = pd.Series(AUC_AUBOC_df.index).replace(namedict)
    AUC_AUBOC_df[bo_auc_list].plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance value\n(observed overlap - random overlap)")
    ax.legend()  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(BOAUC10_barchart_png, dpi=240)
    fig.savefig(BOAUC10_barchart_pdf, dpi=240)


def merge_4_files_alignment_method_deprecated(acc, full_seq, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file):
    """Deprecated method to merge predictions, with lots of checks to ensure sequences are the same.

    Parameters
    ----------
    acc
    full_seq
    train_data_file
    THOIPA_prediction_file
    PREDDIMER_prediction_file
    TMDOCK_prediction_file
    merged_data_xlsx_path
    columns_kept_in_combined_file

    Returns
    -------

    """
    all_files_exist = True
    for path in [train_data_file, THOIPA_prediction_file,PREDDIMER_prediction_file,TMDOCK_prediction_file]:
        if not os.path.isfile(path):
            all_files_exist = False
            sys.stdout.write("{} does not exist".format(path))
            break
    if not all_files_exist:
        sys.stdout.write("\n{} skipped. Input file missing.".format(acc))
        sys.stdout.flush()
        # skip this protein
        return None, None, None

    df_train = pd.read_csv(train_data_file,index_col = 0)
    df_train.index = range(1,df_train.shape[0]+1)
    df_thoipa = pd.read_csv(THOIPA_prediction_file)
    df_preddimer = pd.read_csv(PREDDIMER_prediction_file, index_col = 0)
    df_tmdock = pd.read_csv(TMDOCK_prediction_file, index_col = 0)
    if "closedist" in df_preddimer.columns:
        df_preddimer.rename(columns={"closedist": "PREDDIMER"}, inplace=True)
    if "closedist" in df_tmdock.columns:
        df_tmdock.rename(columns={"closedist": "TMDOCK"}, inplace=True)
    df_train_seq = df_train["residue_name"].str.cat()
    df_thoipa_seq = df_thoipa["residue_name"].str.cat()
    df_preddimer_seq = df_preddimer["residue_name"].str.cat()
    df_tmdock_seq = df_tmdock["residue_name"].str.cat()

    seqlist = [df_train_seq, df_thoipa_seq, df_preddimer_seq, df_tmdock_seq]
    for seq in seqlist:
        if seq not in full_seq:
            sys.stdout.write("Sequence in residue_name column of dataframe is not found in the original df_set sequence."
                             "acc : {}\nTMD_seq : {}\nfull seq in df_set : {}\nall TM sequences in list : {}".format(acc, seq, full_seq, seqlist))
            return None, None, None

    df_train = thoipapy.utils.add_res_num_full_seq_to_df(acc, df_train, df_train_seq, full_seq)
    df_thoipa = thoipapy.utils.add_res_num_full_seq_to_df(acc, df_thoipa, df_thoipa_seq, full_seq)
    df_preddimer = thoipapy.utils.add_res_num_full_seq_to_df(acc, df_preddimer, df_preddimer_seq, full_seq)
    df_tmdock = thoipapy.utils.add_res_num_full_seq_to_df(acc, df_tmdock, df_tmdock_seq, full_seq)

    dfs = pd.DataFrame()

    d = {}
    d["df_train"] = df_train_seq
    d["df_thoipa"] = df_thoipa_seq
    d["df_preddimer"] = df_preddimer_seq
    d["df_tmdock"] = df_tmdock_seq
    dfs["seq"] = pd.Series(d)
    # get the length and length_order of the sequences
    dfs["length"] = dfs.seq.str.len()
    dfs["length_order"] = dfs.length.argsort()

    unique_seq_list = dfs.seq.unique()

    if unique_seq_list.shape[0] == 1:
        # if all sequences match, just use the original index in each separate csv
        df_train["IND"] = df_train["residue_name"] + df_train["residue_num"].astype(str)
        df_train.set_index("IND", inplace=True, drop=False)
        df_thoipa["IND"] = df_thoipa["residue_name"] + df_thoipa["residue_num"].astype(str)
        df_thoipa.set_index("IND", inplace=True)
        df_preddimer["IND"] = df_preddimer["residue_name"] + df_preddimer["residue_num"].astype(str)
        df_preddimer.set_index("IND", inplace=True)
        df_tmdock["IND"] = df_tmdock["residue_name"] + df_tmdock["residue_num"].astype(str)
        df_tmdock.set_index("IND", inplace=True)

    elif unique_seq_list.shape[0] > 2:
        # skip protein if there are more than 3 sequences
        # a multiple sequence alignment would be necessary
        sys.stdout.write("4 sequences has more than 2 unique sequences, alignment not possible. protein {} skipped.".format(acc))
        # skip protein
        return None, None, None

    elif unique_seq_list.shape[0] ==2:
        sys.stdout.write("\n\nstarting reindexing of different TM lengths. ")
        # create a pairwise alignment are reindex dataframes
        for n, seq in enumerate(unique_seq_list):
            # select the rows that match that particular unique sequence
            dfs_sel = dfs.loc[dfs.seq == seq]
            # give them a unique number (0 or 1
            dfs.loc[dfs_sel.index, "unique_num"] = n
        dfs["unique_num"] = dfs["unique_num"].astype(int)

        # Align the two unique sequences(adding gaps at start or end)
        aligned_seq1, aligned_seq2, _, _, _ = \
        pairwise2.align.globalxx(unique_seq_list[0], unique_seq_list[1], one_alignment_only=True)[0]
        seq1 = aligned_seq1.strip("-")
        seq2 = aligned_seq2.strip("-")
        seq1 = aligned_seq1.replace("-", "")
        seq2 = aligned_seq2.replace("-", "")
        # add the sequence (including gaps) to dfs
        dfs.loc[dfs.seq == seq1, "aligned_seq"] = aligned_seq1
        dfs.loc[dfs.seq == seq2, "aligned_seq"] = aligned_seq2

        sys.stdout.write("Check the 4 aligned sequences.")
        for seq in dfs.aligned_seq.tolist():
            sys.stdout.write("\n{}".format(seq))
        # simply count the gaps at the start and the end, for each seq
        dfs["n_of_gaps_at_start_and_end_of_seq"] = dfs.aligned_seq.apply(
            thoipapy.utils.get_n_of_gaps_at_start_and_end_of_seq)
        # number of gaps at the start
        df_train_n_at_start = dfs.loc["df_train", "n_of_gaps_at_start_and_end_of_seq"][0]
        df_thoipa_n_at_start = dfs.loc["df_thoipa", "n_of_gaps_at_start_and_end_of_seq"][0]
        df_preddimer_n_at_start = dfs.loc["df_preddimer", "n_of_gaps_at_start_and_end_of_seq"][0]
        df_tmdock_n_at_start = dfs.loc["df_tmdock", "n_of_gaps_at_start_and_end_of_seq"][0]

        # the final index will be based on the longest unique sequence
        longest_seq_len = dfs.length.max()
        # reindex to add the gaps at the start and the end
        df_train.index = range(df_train_n_at_start, df_train_n_at_start + df_train.shape[0])
        df_thoipa.index = range(df_thoipa_n_at_start, df_thoipa_n_at_start + df_thoipa.shape[0])
        df_preddimer.index = range(df_preddimer_n_at_start, df_preddimer_n_at_start + df_preddimer.shape[0])
        df_tmdock.index = range(df_tmdock_n_at_start, df_tmdock_n_at_start + df_tmdock.shape[0])

        # reindex so that index is now consistent between all 3 dataframes, and any missing rows are added
        df_train = df_train.reindex(index=range(longest_seq_len))
        df_thoipa = df_thoipa.reindex(index=range(longest_seq_len))
        df_preddimer = df_preddimer.reindex(index=range(longest_seq_len))
        df_tmdock = df_tmdock.reindex(index=range(longest_seq_len))

        df_train.dropna(axis=0, how="all", inplace=True)
        df_thoipa.dropna(axis=0, how="all", inplace=True)
        df_preddimer.dropna(axis=0, how="all", inplace=True)
        df_tmdock.dropna(axis=0, how="all", inplace=True)

        df_train["IND"] = df_train["residue_name"] + df_train.index.astype(str)
        df_thoipa["IND"] = df_thoipa["residue_name"] + df_thoipa.index.astype(str)
        df_preddimer["IND"] = df_preddimer["residue_name"] + df_preddimer.index.astype(str)
        df_tmdock["IND"] = df_tmdock["residue_name"] + df_tmdock.index.astype(str)

        df_train.set_index("IND", inplace=True, drop=False)  # keep this one in the final merged df
        df_thoipa.set_index("IND", inplace=True)
        df_preddimer.set_index("IND", inplace=True)
        df_tmdock.set_index("IND", inplace=True)

    dfm = pd.concat([df_train, df_thoipa, df_preddimer, df_tmdock], axis=1, join="outer")

    dfm["aa_pos_in_dfm"] = dfm.index.str[1:].astype(int)
    # start numbering at 1 instead of 0
    dfm["aa_pos_in_dfm"] = dfm["aa_pos_in_dfm"] + 1
    dfm.sort_values("aa_pos_in_dfm", inplace=True)
    #dfm["Polarity"] = dfm["lipo_Hessa"]
    # use -entropy, named as conservation
    #dfm["Conservation"] = -dfm["Entropy"]
    dfm = dfm.loc[:, columns_kept_in_combined_file]

    with pd.ExcelWriter(merged_data_xlsx_path) as writer:
        df_train.to_excel(writer, sheet_name="df_train")
        df_thoipa.to_excel(writer, sheet_name="df_thoipa")
        df_preddimer.to_excel(writer, sheet_name="df_preddimer")
        df_tmdock.to_excel(writer, sheet_name="df_tmdock")
        dfm.to_excel(writer, sheet_name="dfm")

    sys.stdout.write("\n{} finished. Merged data saved to {}".format(acc, merged_data_xlsx_path))
    sys.stdout.flush()

def create_ROC_comp_4predictors(s, df_set, logging):

    logging.info("start create_ROC_Curve_figs_THOIPA_PREDDIMER_TMDOCK_LIPS")
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])
    prediction_name_list = [pred_colname, "PREDDIMER", "TMDOCK", "LIPS_surface_ranked"]
    ROC_4predictor_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                      "{}_ROC_4predictors.csv".format(s['setname']))
    ROC_4predictor_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                               "{}_ROC_4predictors.png".format(s['setname']))
    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    mean_tpr_list=[]
    for n, predictor in enumerate(prediction_name_list):
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        for i in df_set.index:
            acc = df_set.loc[i, "acc"]
            database = df_set.loc[i, "database"]
            merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))
            dfm = pd.read_csv(merged_data_csv_path, engine="python", index_col=0)
            dfm.dropna(inplace=True)
            interface = dfm["interface"].values
            if n ==0 or n == 3:
                predict = dfm[predictor].values
            else:
                predict = -1 * dfm[predictor].values
            fpr, tpr, thresholds = roc_curve(interface, predict, drop_intermediate=False)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0
        mean_tpr /= len(df_set.index)
        mean_tpr[-1] = 1.0
        mean_roc_auc = auc(mean_fpr, mean_tpr)
        mean_tpr_list.append(mean_tpr)
        ax.plot(mean_fpr, mean_tpr, lw=1,label="{} (area = {:.2f})".format(predictor, mean_roc_auc), alpha=0.8)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_4predictor_png, dpi=240)
    # fig.savefig(crossvalidation_png[:-4] + ".pdf")
    fig.savefig(thoipapy.utils.pdf_subpath(ROC_4predictor_png))
    df_tpr = pd.DataFrame.from_records(list(map(list, zip(*mean_tpr_list))),
                                       columns=prediction_name_list)
    df_tpr.to_csv(ROC_4predictor_csv)

def create_linechart_perc_interf_vs_PR_cutoff(s, predictor_name_list, perc_interf_vs_PR_cutoff_linechart_png, database="all"):
    """ Create linechart (and barchart) showing percentage of interface residues correctly predicted, according
    to precision-recall cutoffs.

    Parameters
    ----------
    s : dict
        Settings dictionary
    predictor_name_list : list
        List of predictors to include in plot
    perc_interf_vs_PR_cutoff_linechart_png : str
        Linechart path
    database : str
        Database for which figures are processed (e.g. ETRA). Default is all.

    Returns
    -------

    """

    #######################################################################################################
    #                                                                                                     #
    #                            Linechart at various precision-recall cutoffs                            #
    #                                                                                                     #
    #######################################################################################################

    plt.rcParams.update(plt.rcParamsDefault)

    thoipapy.utils.make_sure_path_exists(perc_interf_vs_PR_cutoff_linechart_png, isfile=True)

    df_PR_cutoff_all = pd.DataFrame()

    # list of cutoffs to plot on x-axis
    cutoff_list = np.arange(0.1, 0.95, 0.05)

    fig, ax = plt.subplots(figsize=(3.42, 3.42))
    for predictor_name in predictor_name_list:

        result_dict = {}

        auc_pkl = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "indiv_validation", predictor_name, "ROC_AUC_data.pkl")
        with open(auc_pkl, "rb") as f:
            xv_dict = pickle.load(f)

        for acc_db in xv_dict:
            if database is not "all":
                if database not in acc_db:
                    # skip the TMDs that are not in the right database
                    continue

            precision = xv_dict[acc_db]["precision"]
            recall = xv_dict[acc_db]["recall"]
            result_list = []
            for cutoff in cutoff_list:
                PR_rises_above_threshold = precision_recall_curve_rises_above_threshold(precision, recall, threshold=cutoff)
                result_list.append(PR_rises_above_threshold)
            result_dict[acc_db] = result_list

        df_cutoffs = pd.DataFrame(result_dict, index = cutoff_list)
        # save percentage of TMDs that exceeds cutoff
        df_PR_cutoff_all[predictor_name] = df_cutoffs.sum(axis=1) / df_cutoffs.shape[1]


    THOIPA_name = [x for x in predictor_name_list if "THOIPA" in x][0]
    df_PR_cutoff_all.columns = pd.Series(df_PR_cutoff_all.columns).replace(THOIPA_name, "THOIPA")

    # round index so values at index 0.5 can be identified
    df_PR_cutoff_all.index = pd.Series(df_PR_cutoff_all.index).round(2)


    #######################################################################################################
    #                                                                                                     #
    #                   Linechart showing percentage TMDs above cutoff vs cutoffs                         #
    #                                                                                                     #
    #######################################################################################################
    colour_dict = {"THOIPA" : "#E95D12", "TMDOCK" : "#0065BD", "PREDDIMER" : "k", "random" : "grey"}
    colour_list = ["#E95D12","#0065BD","k","grey"]
    fontsize = 8
    linewidth = 0.7
    #
    cols = ['THOIPA', 'TMDOCK','PREDDIMER','random']
    #df_PR_cutoff_all = df_PR_cutoff_all.reindex(columns = cols)
    #df_PR_cutoff_all.plot(ax=ax, color=colour_list, linestyle=linestyle_list)
    for col in cols:
        # set dotted line for random
        if col == "random":
            linestyle = "--"
        else:
            linestyle = "-"

        ax.plot(df_PR_cutoff_all.index, df_PR_cutoff_all[col], linestyle=linestyle, color=colour_dict[col], linewidth=linewidth, label=col)


    # general properties
    ax.grid(False)
    ax.set_facecolor('white')
    plt.rcParams["axes.edgecolor"] = "k"
    plt.rcParams["axes.linewidth"] = 0.4
    fig.patch.set_visible(True)
    ax.legend(fontsize=fontsize)

    # axis limitos
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05,1.05)

    # tick properties
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.tick_params(direction='out', length=0, width=1, colors='k')
    ax.tick_params(axis='y', labelsize=fontsize,pad=2)
    ax.tick_params(axis='x', labelsize=fontsize,pad=2)

    ax.set_xlabel("recall and precision cut-off", fontsize=fontsize)
    ax.set_ylabel("fraction of interfaces that meets the cut-off", fontsize=fontsize, labelpad=1)

    fig.tight_layout()
    fig.savefig(perc_interf_vs_PR_cutoff_linechart_png, dpi=300)
    fig.savefig(perc_interf_vs_PR_cutoff_linechart_png[:-4] + ".pdf")

    df_PR_cutoff_all.to_csv(perc_interf_vs_PR_cutoff_linechart_png[:-4] + "_data.csv")

    #######################################################################################################
    #                                                                                                     #
    #                            Barchart at precision-recall cutoff of 0.5                               #
    #                                                                                                     #
    #######################################################################################################
    # bar positions and labels
    x = [0.05, 1.05, 2, 3]
    labels = ['THOIPA', 'TMDOCK','PREDDIMER', 'random']

    df = df_PR_cutoff_all.reindex(columns = labels)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(1.32,2.32))

    fontsize = 6.5
    linewidth = 0.7
    colour_list = ["#E95D12","#0065BD","k","0.5"]
    colour_list_anno = ["k", "white", "white", "k"]

    # plot only the row with a cutoff of 0.5
    df.loc[0.5, :].plot(ax=ax, kind="bar", color = colour_list, width = 0.6)

    # general properties
    ax.grid(False)
    plt.rcParams["axes.edgecolor"] = "k"
    plt.rcParams["axes.linewidth"] = 0.4

    # tick properties
    ax.set_xticks(x)
    ax.set_xticklabels("")
    ax.tick_params(axis='y', labelsize=fontsize,pad=2)
    ax.tick_params(direction='in', length=0, width=0, colors='k')

    #ax.set_ylim(0,0.55)

    ax.set_ylabel('fraction correct interfaces (at 0.5 cutoff)   ', fontsize=fontsize, labelpad=1)

    # text annotations in bars
    for i, txt in enumerate(labels):
        ax.annotate(txt, (x[i],0.005), size=fontsize, color=colour_list_anno[i], ha='center',rotation=90,va="bottom")

    fig.tight_layout()
    bar_png = os.path.join(perc_interf_vs_PR_cutoff_linechart_png[:-13] + "barchart.png")
    fig.savefig(bar_png, dpi=300)
    fig.savefig(bar_png[:-4] + ".pdf")


