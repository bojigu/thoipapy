import os
import sys
from pathlib import Path
from typing import Union

import thoipapy
#import eccpy
import numpy as np
import pandas as pd
import seaborn as sns; sns.set()
from matplotlib import pyplot as plt
#from eccpy.tools import normalise_between_2_values
from thoipapy.utils import normalise_between_2_values

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

    THOIPA_col = "THOIPA_{}_LOO".format(s["set_number"])
    LIPS_col = "LIPS_surface"#"LIPS_surface_ranked"
    coev_col = "DI4mean"

    dfh_cols = ["res_num_full_seq", "residue_name", "interface", "interface_score", THOIPA_col, "PREDDIMER", "TMDOCK", LIPS_col, "conservation", "relative_polarity", coev_col]

    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    df_names = pd.read_excel(names_excel_path, index_col=0)
    df_names["acc_db"] = df_names.index + "_" + df_names["database"]
    df_names["acc"] = df_names.index
    df_names.set_index("acc_db", inplace=True)

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]

        if database == "crystal":
            acc_db = acc + "_X-ray"
        else:
            acc_db = acc + "_" + database
        shortname = df_names.loc[acc_db, "shortname"]
        uniprot = df_names.loc[acc_db, "uniprot"]

        if database == "ETRA":
            ref = "".join(df_names.loc[acc_db, "source":"date"].dropna().astype(str).tolist())
            savename = "{}_{}".format(acc, shortname)
            fig_label = "{shortname} [{subset} dataset, {acc}, {ref}]".format(shortname=shortname,
                                                                            subset=database, acc=acc, ref=ref)
        elif database == "NMR":
            ref = "".join(df_names.loc[acc_db, "source":"date"].dropna().astype(str).to_list())
            savename = "{}_{}".format(acc, shortname)
            fig_label = "{shortname} [{subset} dataset, {acc}, PDB:{pdb}, {ref}]".format(shortname=shortname,
                        subset=database, acc=acc, pdb=df_names.loc[acc_db, "PDB acc"],ref=ref)
        elif database == "crystal":
            savename = acc + "_".format(database)
            #fig_label = acc + " [{} subset, PDB:{}, chain:{}, TMD:{}]".format(database, acc[:-2], acc[-2], acc[-1])
            fig_label = "{} [X-ray dataset, {}, PDB:{}, chain:{}, TMD:{}]]".format(shortname, uniprot, acc[:-2], acc[-2], acc[-1])
        else:
            raise ValueError("database not recognised : {}".format(database))

        sys.stdout.write("\n{} {}".format(savename, fig_label))
        create_single_merged_heatmap(s, acc, database,savename, fig_label, dfh_cols, THOIPA_col, LIPS_col, coev_col)

def create_single_merged_heatmap(s, acc, database, savename, fig_label, dfh_cols, THOIPA_col, LIPS_col, coev_col):
        merged_data_csv_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/predictions/merged/{database}.{acc}.merged.csv"
        dfm = pd.read_csv(merged_data_csv_path, engine = "python")


        heatmap_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database, "{}.png".format(acc))
        heatmap_pdf_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database,"pdf","{}.pdf".format(acc))
        heatmap_data_xlsx_path = os.path.join(s["thoipapy_data_folder"], "heatmap",database,"xlsx","{}_merged.xlsx".format(acc))
        hetero_bind_file = os.path.join(s["thoipapy_data_folder"], "Features", "Structure",database, "{}.hetero.bind.csv").format(acc)
        thoipapy.utils.make_sure_path_exists(heatmap_pdf_path, isfile=True)
        thoipapy.utils.make_sure_path_exists(heatmap_data_xlsx_path, isfile=True)

        # create dfh, dataframe for heatmap
        dfh = dfm[dfh_cols].copy()

        # Drop only positions where there is no interface data
        dfh.dropna(subset=["interface"], inplace=True)
        # why reset index here???
        dfh.reset_index(drop=True, inplace=True)
        # set index to start with 1
        dfh.index = range(1, dfh.shape[0] + 1)
        # normalise all the data columns between 0 and 1
        #cols_to_plot = dfh_cols[2:]
        dfh["PREDDIMER_norm"] = normalise_between_2_values(dfh["PREDDIMER"], 2.5, 8, invert=True)
        #dfm["PREDDIMER"] = -1 * dfm["PREDDIMER"]
        dfh["TMDOCK_norm"] = normalise_between_2_values(dfh["TMDOCK"], 2.5, 8, invert=True)
        #dfm["TMDOCK"] = -1 * dfm["TMDOCK"]
        if database == "crystal" or database == "NMR":
            # normalize crystal and NMR closedistance to between 0 and 1 with invert, min and max values were set as 2 and 10 angstrom
            dfh["interface score_norm"] = normalise_between_2_values(dfh["interface_score"],2,10,invert=True)
            #dfm["interface_score"] = -1 * dfm["interface_score"]
        elif database == "ETRA":
            ###normalize ETRA experimental disruption value to the range of 0 to 1 without invert, the min and max values were set as -0.4 and 0.4
            dfh["interface score_norm"] = normalise_between_2_values(dfh["interface_score"], -0.4, 0.4)

        # norm conservation
        dfh["conservation_norm"] = normalise_between_2_values(dfh["conservation"], 1.25, 3, invert=False)
        # norm THOIPA
        dfh["THOIPA_norm"] = normalise_between_2_values(dfh[THOIPA_col], 0.15, 0.5, invert=False)
        # norm polarity
        dfh["relative polarity_norm"] = normalise_between_2_values(dfh["relative_polarity"], 0.5, 2.5, invert=False)
        # norm LIPS
        #dfh["LIPS_norm"] = normalise_between_2_values(dfh[LIPS_col], -0.4, 1, invert=False)
        dfh["LIPS_norm"] = dfh[LIPS_col]

        # currently coevolution doesn't need to be normalised (already norm in each TMD between 0 and 1
        dfh["coevolution_norm"] = dfh[coev_col]

        cols_to_plot = ['interface score_norm', 'interface', 'PREDDIMER_norm', 'TMDOCK_norm', 'THOIPA_norm', 'LIPS_norm', 'conservation_norm', 'relative polarity_norm', 'coevolution_norm']
        cols_to_plot_renamed = [x[:-5] if "_norm" in x else x for x in cols_to_plot]

        # transpose dataframe so that "interface" etc is on the left
        dfh_to_plot = dfh[cols_to_plot].T
        dfh_to_plot.index = cols_to_plot_renamed
        # create dataframe that has an "-" in areas of np.nan
        # this is used to label the missing residues in TMDOCK
        df_labels = dfh_to_plot.isnull().replace(False, "")
        df_labels = df_labels.replace(True, "-")

        #################################################################################
        #                                                                               #
        #          @BO RECOMMEND INSERTING CODE TO LABEL "FOLDING" residues here        #
        #                                                                               #
        #################################################################################

        """
        for res in wherever_your_hetero(folding)_residues_are:
            if res == folding:
                df_labels.loc["interface_score", res] = "*"
        
        """
        if os.path.exists(hetero_bind_file):
            df_hetero = pd.read_csv(hetero_bind_file,engine="python")
            for i in df_hetero.index:
                if df_hetero.iloc[i]["hetero_interface"] == 1:
                    df_labels.loc["interface score",i+1] = "*"

        # now replace np.nan with 0 in original shading dataframe (colour will look like 0, rather than white)
        dfh_to_plot.fillna(0, inplace=True)

        fontsize = 16
        tum_blue4_as_python_color = np.array([0, 82, 147]) / 255
        cmap = sns.light_palette(tum_blue4_as_python_color, as_cmap=True)


        plt.close("all")
        # sns.set_context("paper", rc={"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10})
        """
        IMPORTANT!!
        The default fontsize controls the spacing between the subplots, EVEN IF THERE ARE NO TITLES or XLABELS!      
        """
        plt.rcParams['font.size'] = fontsize/2
        # relative height of each subplot heatmap
        gridspec_kw = {"height_ratios": [4, 8, 6]}

        # create 3 subplots
        fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(10, 3), gridspec_kw=gridspec_kw)

        # SOMEWHAT INELEGANT: In order to create a second xticklabels, a second axis is created, and each heatmap rendered twice
        # in ax2_2, the aa numbers are used as the xticklabels at the bottom
        # in ax2_0, the amino acid letters are used as the xticklabels on the top
        ax2_0 = axes[0].twiny()
        ax2_1 = axes[1].twiny()
        ax2_2 = axes[2].twiny()

        # plot the same data in main axis, and twinx
        # interface and interface_score
        sns.heatmap(dfh_to_plot[0:2], ax=axes[0], xticklabels=False, cbar=False,
                    cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        sns.heatmap(dfh_to_plot[0:2], ax=ax2_0, cbar=False, cmap=cmap, annot=df_labels[0:2], fmt="s",
                    annot_kws={"color": "k", "fontsize" : fontsize, "verticalalignment" : "top"})

        # PREDDIMER, TMDOCK, THOIPA, LIPS
        sns.heatmap(dfh_to_plot[2:6], ax=axes[1], xticklabels=False, cbar=False,
                    cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        sns.heatmap(dfh_to_plot[2:6], ax=ax2_1, cbar=False, xticklabels=[''] * dfh_to_plot.shape[1],
                    cmap=cmap, annot=df_labels[2:6], fmt="s", annot_kws={"color": "k", "fontsize" : fontsize})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,

        # conservation, relative polarity, coevolution
        sns.heatmap(dfh_to_plot[6:9], ax=axes[2], cbar=False,
                    cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        sns.heatmap(dfh_to_plot[6:9], ax=ax2_2, cbar=False, xticklabels=[''] * dfh_to_plot.shape[1],
                    cmap=cmap, annot=df_labels[6:9], fmt="s", annot_kws={"color": "k", "fontsize" : fontsize})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,

        # set fontsize and rotation of y-labels
        axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=fontsize, rotation=0)
        axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=fontsize, rotation=0)
        axes[2].set_yticklabels(axes[2].get_yticklabels(), fontsize=fontsize, rotation=0)

        # bottom residue numbers
        axes[2].set_xticklabels(dfh.index, fontsize=fontsize, rotation=0)
        axes[2].tick_params(axis="x", direction='out', pad=1.5, tick2On=False)

        # figure title
        # residue letters at top
        ax2_0.set_xlabel(fig_label, fontsize=fontsize)
        ax2_0.set_xticks(axes[2].get_xticks())
        ax2_0.xaxis.tick_top()
        ax2_0.set_xticklabels(dfh.residue_name, fontsize=fontsize)
        ax2_0.tick_params(axis="x", direction='out', pad=-0.1, tick2On=False)

        plt.tight_layout()
        fig.savefig(heatmap_path, dpi=240)
        fig.savefig(heatmap_pdf_path)

        with pd.ExcelWriter(heatmap_data_xlsx_path) as writer:
            dfm.to_excel(writer, sheet_name="dfm")
            dfh.to_excel(writer, sheet_name="dfh")
            dfh_to_plot.to_excel(writer, sheet_name="dfh_to_plot")

        sys.stdout.write("\n{} heatmap finished. ({})".format(savename, heatmap_path))
        sys.stdout.flush()




        #
        # plt.close("all")
        # # sns.set_context("paper", rc={"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10})
        # """
        # IMPORTANT!!
        # The default fontsize controls the spacing between the subplots, EVEN IF THERE ARE NO TITLES or XLABELS!
        # """
        # plt.rcParams['font.size'] = fontsize / 2
        # # relative height of each subplot heatmap
        # gridspec_kw = {"height_ratios": [4, 8, 6]}
        #
        # # create 3 subplots
        # fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(10, 3), gridspec_kw=gridspec_kw)
        #
        # # fig.subplots_adjust(hspace=-1)
        # # fig.subplots_adjust(bottom=0.1, top=0.1)
        #
        # # SOMEWHAT INELEGANT: In order to create a second xticklabels, a second axis was created, and the heatmap rendered twice
        # # in ax1, the aa numbers are used as the xticklabels at the bottom
        # # in ax2, the amino acid letters are used as the xticklabels on the top
        # # ax2 = ax.twiny()
        # ax2_0 = axes[0].twiny()
        # ax2_1 = axes[1].twiny()
        # ax2_2 = axes[2].twiny()
        # # ax2_2 = axes[2].twiny()
        # # ax = axes[2]
        # # plot in ax and ax2
        # # sns.heatmap(dfh_to_plot, ax=ax, cbar=False, cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot, ax=ax2, cbar=False, cmap=cmap, annot=df_labels, fmt="s", annot_kws={"color": "k"})
        # sns.heatmap(dfh_to_plot[0:2], ax=axes[0], xticklabels=False, cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[0:2], ax=ax2_0, cbar=False, cmap=cmap, annot=df_labels[0:2], fmt="s",
        #             annot_kws={"color": "k", "fontsize": fontsize, "verticalalignment": "top"})
        # sns.heatmap(dfh_to_plot[2:6], ax=axes[1], xticklabels=False, cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[2:6], ax=ax2_1, cbar=False,
        #             cmap=cmap, annot=df_labels[2:6], fmt="s", annot_kws={"color": "k", "fontsize": fontsize})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[6:9], ax=axes[2], cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[6:9], ax=ax2_2, cbar=False, xticklabels=[''] * dfh_to_plot.shape[1],
        #             cmap=cmap, annot=df_labels[6:9], fmt="s", annot_kws={"color": "k", "fontsize": fontsize})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[0:2,:], ax=axes[0], xticklabels=False, cbar=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[0:2,:], ax=ax2 , cbar=False, cmap=cmap, annot=df_labels[0:2], fmt="s",
        # #             annot_kws={"color": "k"})
        # # sns.heatmap(dfh_to_plot.iloc[2:6,:], ax=axes[1], cbar=False,xticklabels=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # # sns.heatmap(dfh_to_plot.iloc[2:6,:], ax=ax2_1 ,cbar=False, cmap=cmap, annot=df_labels.iloc[2:6,:], fmt="s",
        # # #             annot_kws={"color": "k"})
        # # sns.heatmap(dfh_to_plot.iloc[6:9,:], ax=axes[2], cbar=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[6:9, :], ax=ax2_2,  cbar=False, cmap=cmap, annot=df_labels[6:9], fmt="s",
        # #             annot_kws={"color": "k"})
        # # set aa position and letter labels
        #
        # axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=fontsize, rotation=0)
        # axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=fontsize, rotation=0)
        # axes[2].set_yticklabels(axes[2].get_yticklabels(), fontsize=fontsize, rotation=0)
        # # ax2_2.set_yticklabels(df_labels[6:9].index, fontsize=5)
        #
        # axes[2].set_xticklabels(dfh.index, fontsize=fontsize, rotation=0)
        # axes[2].tick_params(axis="x", direction='out', pad=1.5, tick2On=False)
        # # ax2_1.set_yticklabels(axes[1].get_yticklabels(), fontsize=fontsize / 10, rotation=0)
        #
        #
        # # axes[2].set_xlabel("")
        # # ax2.set_xlim(ax.get_xlim())
        # # ax2_0.set_ylim(ax2_0.get_ylim())
        # # axes[2].set_ylim(axes[2].get_ylim())
        #
        # ax2_0.set_xlabel(fig_label, fontsize=fontsize)
        # ax2_0.set_xticks(axes[2].get_xticks())
        # ax2_0.xaxis.tick_top()
        # ax2_0.set_xticklabels(dfh.residue_name, fontsize=fontsize)
        # ax2_0.tick_params(direction='out', pad=-0.1, tick2On=False)
        #
        # # axes[2].tick_params(direction='out', pad=1.8, tick1On=False)
        # # ax2_0.tick_params(direction='out', pad=0.1, tick2On=False)
        #
        # # ax2_2.tick_params("x", direction='out', pad=-0.1, tick2On=False)
        #
        #
        #
        #
        # # attempt to remove any labels causing spaces
        # axes[1].set_xticklabels("", fontsize=fontsize / 10, rotation=0)
        # ax2_1.set_xticklabels("", fontsize=fontsize / 10, rotation=0)
        # axes[1].set_title("", fontsize=fontsize / 10)
        # ax2_1.set_title("", fontsize=fontsize / 10)
        # axes[1].set_xticklabels("", fontsize=fontsize / 10)
        # ax2_1.set_xticklabels("", fontsize=fontsize / 10)
        # axes[1].set_xlabel("", fontsize=fontsize / 10)
        # ax2_1.set_xlabel("", fontsize=fontsize / 10)
        #
        # plt.tight_layout()
        # fig.savefig(heatmap_path, dpi=240)
        # fig.savefig(heatmap_pdf_path)
        #
        # with pd.ExcelWriter(heatmap_data_xlsx_path) as writer:
        #     dfm.to_excel(writer, sheet_name="dfm")
        #     dfh.to_excel(writer, sheet_name="dfh")
        #     dfh_to_plot.to_excel(writer, sheet_name="dfh_to_plot")
        #
        # sys.stdout.write("\n{} heatmap finished. ({})".format(savename, heatmap_path))
        # sys.stdout.flush()


        #
        #
        # plt.close("all")
        # # sns.set_context("paper", rc={"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10})
        # plt.rcParams['font.size'] = fontsize
        # gridspec_kw = {"height_ratios": [4, 8, 6]}
        # # gridspec_kw = {"height_ratios": [2, 4, 3]}
        #
        # # fig, ax = plt.subplots(figsize=(3.42, 1))
        # fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(10, 3), gridspec_kw=gridspec_kw)
        # fig.subplots_adjust(hspace=-1)
        # # fig.subplots_adjust(bottom=0.1, top=0.1)
        #
        # # SOMEWHAT INELEGANT: In order to create a second xticklabels, a second axis was created, and the heatmap rendered twice
        # # in ax1, the aa numbers are used as the xticklabels at the bottom
        # # in ax2, the amino acid letters are used as the xticklabels on the top
        # # ax2 = ax.twiny()
        # ax2 = axes[0].twiny()
        # ax3 = axes[1].twiny()
        # ax4 = axes[2].twiny()
        # # ax4 = axes[2].twiny()
        # ax = axes[2]
        # # plot in ax and ax2
        # # sns.heatmap(dfh_to_plot, ax=ax, cbar=False, cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot, ax=ax2, cbar=False, cmap=cmap, annot=df_labels, fmt="s", annot_kws={"color": "k"})
        # sns.heatmap(dfh_to_plot[0:2], ax=axes[0], xticklabels=False, cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[0:2], ax=ax2, cbar=False, cmap=cmap, annot=df_labels[0:2], fmt="s",
        #             annot_kws={"color": "k"})
        # sns.heatmap(dfh_to_plot[2:6], ax=axes[1], xticklabels=False, cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[2:6], ax=ax3, cbar=False, xticklabels=[''] * dfh_to_plot.shape[1],
        #             cmap=cmap, annot=df_labels[2:6], fmt="s", annot_kws={"color": "k"})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[6:9], ax=axes[2], cbar=False,
        #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # sns.heatmap(dfh_to_plot[6:9], ax=ax4, cbar=False, xticklabels=[''] * dfh_to_plot.shape[1],
        #             cmap=cmap, annot=df_labels[6:9], fmt="s", annot_kws={"color": "k"})  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[0:2,:], ax=axes[0], xticklabels=False, cbar=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[0:2,:], ax=ax2 , cbar=False, cmap=cmap, annot=df_labels[0:2], fmt="s",
        # #             annot_kws={"color": "k"})
        # # sns.heatmap(dfh_to_plot.iloc[2:6,:], ax=axes[1], cbar=False,xticklabels=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # # sns.heatmap(dfh_to_plot.iloc[2:6,:], ax=ax3 ,cbar=False, cmap=cmap, annot=df_labels.iloc[2:6,:], fmt="s",
        # # #             annot_kws={"color": "k"})
        # # sns.heatmap(dfh_to_plot.iloc[6:9,:], ax=axes[2], cbar=False,
        # #             cmap=cmap)  # fmt = "s", annot_kws={"Axes.set_facecolor", 0.5} ,
        # # sns.heatmap(dfh_to_plot.iloc[6:9, :], ax=ax4,  cbar=False, cmap=cmap, annot=df_labels[6:9], fmt="s",
        # #             annot_kws={"color": "k"})
        # # set aa position and letter labels
        # axes[2].set_xticklabels(dfh.index, fontsize=fontsize, rotation=0)
        # axes[2].set_yticklabels(axes[2].get_yticklabels(), fontsize=fontsize, rotation=0)
        # axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=fontsize, rotation=0)
        # axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=fontsize, rotation=0)
        #
        # ax.set_xlabel("")
        # # ax2.set_xlim(ax.get_xlim())
        # ax2.set_ylim(ax2.get_ylim())
        # ax.set_ylim(ax.get_ylim())
        #
        # ax2.set_xlabel(fig_label, fontsize=fontsize)
        # ax2.set_xticks(ax.get_xticks())
        # ax2.xaxis.tick_top()
        # ax2.set_xticklabels(dfh.residue_name, fontsize=fontsize)
        #
        # ax.tick_params(direction='out', pad=1.8, tick1On=False)
        # # ax2.tick_params(direction='out', pad=0.1, tick2On=False)
        # ax2.tick_params(direction='out', pad=-0.1, tick2On=False)
        # plt.tight_layout()
        # fig.savefig(heatmap_path, dpi=240)
        # fig.savefig(heatmap_pdf_path)
        #
        # with pd.ExcelWriter(heatmap_data_xlsx_path) as writer:
        #     dfm.to_excel(writer, sheet_name="dfm")
        #     dfh.to_excel(writer, sheet_name="dfh")
        #     dfh_to_plot.to_excel(writer, sheet_name="dfh_to_plot")
        #
        # sys.stdout.write("\n{} heatmap finished. ({})".format(savename, heatmap_path))
        # sys.stdout.flush()
