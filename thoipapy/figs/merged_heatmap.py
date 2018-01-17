import os
import sys
import thoipapy
import eccpy
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def create_merged_heatmap(s, df_set, logging):
    """Create heatmap from merged disruption, combined_prediction, and features in  traindata.csv files.
        Parameters
        ----------
        df_set : dict
            dictionary contains the set file contents
        s : dict
            settings dictionary

        Output files
        ------------
        heatmap : png
            heatmap_path = os.path.join(p["BZ12_2file_heatmap"], "{}.png".format(savename))
        heatmap : pdf
            heatmap_pdf_path = os.path.join(p["BZ12_2file_heatmap_pdf"], "{}.pdf".format(savename))
        heatmap_data : xlsx
            heatmap_data_xlsx_path = os.path.join(p["BZ12_2file_data"], "{}_heatmap_data.xlsx".format(savename))
        """
    sys.stdout.write(
        '\n~~~~~~~~~~~~                 starting create_heatmap_from_merged_files              ~~~~~~~~~~~~\n')
    sys.stdout.flush()
    #################################################################
    #             EXTRACT NAMES FROM NAMES EXCEL FILE               #
    #################################################################
    names_excel_path = os.path.join(s["dropbox_dir"],"ETRA_NMR_names.xlsx")
    df_names = pd.read_excel(names_excel_path, index_col=0)
    # restrict names dict to only that database

    dfh_cols = ["res_num_full_seq", "residue_name","interface", "interface_score", "THOIPA_5_LOO", "PREDDIMER", "TMDOCK", "LIPS_surface_ranked", "Conservation", "polarity", "coev_i4_DI"]
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        #if acc =="2axtM1":
        database = df_set.loc[i, "database"]
        df_names = df_names.loc[df_names.database == database]
        if acc in df_names.index:
            savename = "{}_{}".format(acc, df_names.loc[acc, "shortname"])
            fig_label = df_names.loc[acc, "concatname"]
        else:
            savename = acc
            fig_label = acc
        create_single_merged_heatmap(s, acc, database,savename, fig_label,dfh_cols)

def create_single_merged_heatmap(s, acc, database,savename, fig_label, dfh_cols):
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Merged", database, "{}.merged.csv".format(acc))
        dfm = pd.read_csv(merged_data_csv_path, engine = "python")
        dfm["LIPS_L*E"] = -1 * dfm["LIPS_L*E"]
        dfm["PREDDIMER"] = -1 * dfm["PREDDIMER"]
        dfm["TMDOCK"] = -1 * dfm["TMDOCK"]
        if database == "crystal" or database == "NMR":
            # (it is closest distance and low value means high propencity of interfacial)
            dfm["interface_score"] = -1 * dfm["interface_score"]

        heatmap_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database, "{}.png".format(savename))
        heatmap_pdf_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database,"pdf","{}.pdf".format(savename))
        heatmap_data_xlsx_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database,"xlsx","{}_merged.xlsx".format(savename))
        thoipapy.utils.make_sure_path_exists(heatmap_pdf_path, isfile=True)
        thoipapy.utils.make_sure_path_exists(heatmap_data_xlsx_path, isfile=True)
        # create dfh, dataframe for heatmap
        dfh = dfm[dfh_cols].copy()
        dfh.dropna(inplace=True)
        dfh.reset_index(drop=True, inplace=True)
        # normalise all the data columns between 0 and 1
        cols_to_plot = dfh_cols[2:]
        for col in cols_to_plot:
            dfh[col] = eccpy.tools.normalise_0_1(dfh[col])[0]

        # transpose dataframe so that "disruption" etc is on the left
        dfh_to_plot = dfh[cols_to_plot].T
        df_labels = dfh_to_plot.isnull().replace(False, "")
        df_labels = df_labels.replace(True, "X")
        dfh_to_plot.fillna(0, inplace=True)

        fontsize = 10
        tum_blue4_as_python_color = np.array([0, 82, 147]) / 255
        cmap = sns.light_palette(tum_blue4_as_python_color, as_cmap=True)

        plt.close("all")
        # sns.set_context("paper", rc={"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10})
        plt.rcParams['font.size'] = 10
        # fig, ax = plt.subplots(figsize=(3.42, 1))
        fig, ax = plt.subplots(figsize=(7, 2))

        # SOMEWHAT INELEGANT: In order to create a second xticklabels, a second axis was created, and the heatmap rendered twice
        # in ax1, the aa numbers are used as the xticklabels at the bottom
        # in ax2, the amino acid letters are used as the xticklabels on the top
        ax2 = ax.twiny()
        # plot in ax and ax2
        sns.heatmap(dfh_to_plot, ax=ax, cbar=False, cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        sns.heatmap(dfh_to_plot, ax=ax2, cbar=False, cmap=cmap, annot=df_labels, fmt="s", annot_kws={"color": "k"})
        # set aa position and letter labels
        ax.set_xticklabels(dfh.index, fontsize=fontsize, rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize, rotation=0)

        ax.set_xlabel("")
        # ax2.set_xlim(ax.get_xlim())
        # ax2.set_ylim(ax.get_ylim())

        ax2.set_xlabel(fig_label)
        ax2.set_xticks(ax.get_xticks())
        ax2.xaxis.tick_top()
        ax2.set_xticklabels(dfh.residue_name, fontsize=fontsize)
        ax.tick_params(direction='out', pad=1, tick1On=False)
        ax2.tick_params(direction='out', pad=0.1, tick2On=False)
        fig.tight_layout()
        fig.savefig(heatmap_path, dpi=240)
        fig.savefig(heatmap_pdf_path)

        with pd.ExcelWriter(heatmap_data_xlsx_path) as writer:
            dfm.to_excel(writer, sheet_name="dfm")
            dfh.to_excel(writer, sheet_name="dfh")
            dfh_to_plot.to_excel(writer, sheet_name="dfh_to_plot")

        sys.stdout.write("\n{} heatmap finished.".format(savename))
        sys.stdout.flush()

