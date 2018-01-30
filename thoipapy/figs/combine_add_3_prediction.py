import os
import sys
import pandas as pd
import numpy as np
from Bio import pairwise2
import thoipapy
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
import pickle
import re


def combine_file_add_PREDDIMER_TMDOCK_THOIPA_prediction(s, df_set, logging):
    """Combines all available predictions for a particular testset.

    The testset is determined by the original set_number, not the "test_datasets" list.

    The combined predictions file is saved in the thoipapy/Merged folder, so as not to
    be confused with the "combined" folder that holds only features and interface_score.

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.


    """
    # columns_kept_in_combined_file = ['residue_num', 'residue_name', 'conservation', 'lipo_Hessa', 'CoevDImax_norm',
    #                                  'CoevDI4_norm',
    #                                  'CoevDI8_norm', 'CoevMImax_norm', 'CoevMI4_norm', 'CoevMI8_norm',
    #                                  'CumDI4_norm', 'CumDI8_norm', 'CumMI4_norm', 'CumMI8_norm', 'RelPos_TMD',
    #                                  'RelPos_fullseq', 'LIPS_L*E', 'LIPS_surface_ranked', 'Hydrophobicity_sAA',
    #                                  'TMDOCK', 'PREDDIMER']
    # add the THOIPA prediction name to the list of columns to keep
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])
    # for simplicity, keep only the predictions. Since the index is unique, it can be added later to the combined file.
    columns_kept_in_combined_file = ['residue_num', 'residue_name', pred_colname, 'TMDOCK', 'PREDDIMER','interface','interface_score',"LIPS_surface_ranked", 'LIPS_L*E',"polarity","conservation","coev_i4_DI"]

    #set_list = thoipapy.figs.fig_utils.get_set_lists(s)
    PREDDIMER_TMDOCK_folder = os.path.join(s["base_dir"], "figs", "FigBZ18-PreddimerTmdockComparison")
    #for set_number in set_list:
    #setname = "set{:02d}".format(int(set_number))
    #set_path = thoipapy.common.get_path_of_protein_set(setname, s["sets_folder"])
    #df_set = pd.read_excel(set_path, sheetname="proteins")
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        #if acc =="2axtM1":
        full_seq = df_set.loc[i, "full_seq"]
        database = df_set.loc[i, "database"]
        train_data_file = os.path.join(s["features_folder"], "combined", database,"{}.surr20.gaps5.combined_features.csv".format(acc))
        combined_data_file = os.path.join(s["dropbox_dir"], "THOIPA_data","Features","combined",database,
                                       "{}.surr20.gaps5.combined_features.csv".format(acc))
        thoipapy.utils.make_sure_path_exists(combined_data_file,isfile=True)
        #THOIPA_prediction_file = os.path.join(s["thoipapy_data_folder"], "Predictions", "testset_trainset",database, "{}.THOIPA.trainset04.csv".format(acc))
        #THOIPA_prediction_file = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.LOO.prediction.csv".format(acc, s["setname"]))
        THOIPA_prediction_csv = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.{}.LOO.prediction.csv".format(acc, database, s["setname"]))
        PREDDIMER_prediction_file = os.path.join(PREDDIMER_TMDOCK_folder, database, "{}.preddimer.closedist.csv".format(acc))
        TMDOCK_prediction_file = os.path.join(PREDDIMER_TMDOCK_folder, database, "{}.tmdock.closedist.csv".format(acc))
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Merged", database, "{}.merged.csv".format(acc))
        merged_data_xlsx_path = os.path.join(s["thoipapy_data_folder"], "Merged", database, "{}.merged.xlsx".format(acc))
        thoipapy.utils.make_sure_path_exists(merged_data_xlsx_path, isfile=True)
        #merge_4_files_ALIGNMENT_METHOD(acc, full_seq, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file)

        # load the full feature file as the start of dfm
        dfm = pd.read_csv(train_data_file)
        # set the unique index, based on the residue number in the full sequence
        dfm.set_index("res_num_full_seq", inplace=True)
        #dfm["conservation"] = -1 * dfm["Entropy"]
        file_list = [THOIPA_prediction_csv, PREDDIMER_prediction_file, TMDOCK_prediction_file]
        prediction_name_list = [pred_colname, "PREDDIMER", "TMDOCK"]
        n_files_merged = 0
        for n, file in enumerate(file_list):
            prediction_name = [prediction_name_list[n]]
            if os.path.isfile(file):
                df = pd.read_csv(file, index_col=0)
                seq = df["residue_name"].str.cat()
                if seq not in full_seq:
                    logging.warning(prediction_name)
                    logging.warning("Sequence in residue_name column of dataframe is not found in the original df_set sequence."
                                     "\nacc : {}\nfile number : {}\nTMD_seq : {}\nfull seq in df_set : {}\nTHOIPA_prediction_csv:{}".format(acc, n, seq, full_seq, THOIPA_prediction_csv))
                    if prediction_name == [pred_colname]:
                        df = thoipapy.utils.add_mutation_missed_residues_with_na(s,acc,database,df)
                        seq = df["residue_name"].str.cat()
                    # skip protein
                    #continue

                # add the residue number in the full sequence
                df = thoipapy.utils.add_res_num_full_seq_to_df(acc, df, seq, full_seq)

                if n == 0:
                    #logging.info(df.columns)
                    # the thoipa prediction file has the residue_num as the index, similar to the features
                    df.drop(["residue_name"], axis=1, inplace=True)
                else:
                    # the preddimer and TMDOCK files have a range as the index
                    df.drop(["residue_num", "residue_name"], axis=1, inplace=True)
                    df.columns = pd.Series(df.columns).replace({"closedist" : prediction_name})

                # set the unique index, based on the residue number in the full sequence
                df.set_index("res_num_full_seq", inplace=True)

                # merge the growing dfm file. All rows are included
                dfm = pd.concat([dfm, df], axis=1, join="outer")

                n_files_merged += 1
        # keep the desired columns
        #dfm.to_csv(combined_data_file)
        new_columns_kept_in_combined_file = list(set(columns_kept_in_combined_file).intersection(set(dfm.columns)))
        dfm = dfm[new_columns_kept_in_combined_file]
        # save to "Merged" folder, so as not to get confused with the "combined" files
        dfm.to_csv(merged_data_csv_path)
        logging.info("{} predictions combined. n_files_merged : {}. ({})".format(acc, n_files_merged, merged_data_csv_path))

def create_AUC_BoAUC_figs_THOIPA_PREDDIMER_TMDOCK(s,df_set,logging):

    logging.info("start create_AUC_BoAUC_figs_THOIPA_PREDDIMER_TMDOCK")

    names_excel_path = os.path.join(os.path.dirname(s["sets_folder"]), "ETRA_NMR_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)
    predictor_name_list = ["THOIPA_{}_LOO".format(s["set_number"]),"PREDDIMER", "TMDOCK", "LIPS_surface_ranked"] #"LIPS_L*E",
    AUC_AUBOC_df = pd.DataFrame()
    AUC_AUBOC_name_list = []
    linechar_name_list = []
    AUBOC10_list = []
    df_o_minus_r_mean_df = pd.DataFrame()
    auc_mean_list=[]
    AUC_AUBOC_file = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_AUC_AUBOC.csv".format(s["setname"]))
    mean_AUBOC_file = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_mean_AUBOC.csv".format(s["setname"]))
    predictors_AUC_barchart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_AUC_barchart.png".format(s["setname"]))
    predictors_BOAUC10_barchart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_BOAUC_barchart.png".format(s["setname"]))
    predictors_BOCURVE_linechart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_BOCURVE_linechart.png".format(s["setname"]))
    predictors_mean_auc_barchart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.4predictors_mean_auc_barchart.png".format(s["setname"]))
    for predictor_name in predictor_name_list:
        BO_data_df = pd.DataFrame()
        xv_dict = {}
        auc_dict = {}
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        auc_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.{}_AUC_data.pkl".format(s["setname"],predictor_name.replace('*',"")))
        BO_curve_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                "{}.{}_BO_Curve_data.csv".format(s["setname"],predictor_name.replace('*',"")))
        thoipapy.utils.make_sure_path_exists(BO_curve_data_csv, isfile=True)
        BO_data_excel = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                            "{}.{}_BO_curve_data.xlsx".format(s['setname'],predictor_name.replace('*',"")))
        BO_linechart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                               "{}.{}_BO_linechart.png".format(s['setname'],predictor_name.replace('*',"")))
        BO_barchart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                              "{}.{}_AUBOC10_barchart.png".format(s['setname'],predictor_name.replace('*',"")))
        df_o_minus_r_mean_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                              "{}.df_o_minus_r_mean.csv".format(s['setname']))

        for i in df_set.index:
            sys.stdout.write(".")
            sys.stdout.flush()
            acc = df_set.loc[i, "acc"]

            database = df_set.loc[i, "database"]
            acc_db = acc + "-" + database
            merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Merged", database, "{}.merged.csv".format(acc))
            merged_data_df = pd.read_csv(merged_data_csv_path,engine="python")
            merged_data_df["LIPS_L*E"] = -1 * merged_data_df["LIPS_L*E"]
            merged_data_df["PREDDIMER"] = -1 * merged_data_df["PREDDIMER"]
            merged_data_df["TMDOCK"] = -1 * merged_data_df["TMDOCK"]

            # Mark_is_testing_something = True
            # if Mark_is_testing_something:
            #     merged_data_df = merged_data_df.dropna(subset=["interface_score", "PREDDIMER"])

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

            fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name])
            auc_value = auc(fpr, tpr)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0
            auc_dict[acc_db] = auc_value
            xv_dict[acc_db] = {"fpr": fpr, "tpr": tpr, "auc": auc_value}

        # save dict as pickle
        with open(auc_pkl, "wb") as f:
            pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        BO_data_df.to_csv(BO_curve_data_csv)
        # THOIPA_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(THOIPA_BO_curve_data_csv, BO_curve_folder, names_excel_path)
        thoipapy.figs.Create_Bo_Curve_files.parse_BO_data_csv_to_excel(BO_curve_data_csv, BO_data_excel, logging, predictor_name)
        AUC_ser = pd.Series(auc_dict)
        AUC_ser.sort_values(inplace=True, ascending=False)
        auc_mean_list.append(AUC_ser.mean())
        AUBOC10_ser = pd.read_excel(BO_data_excel, sheetname="AUBOC10", index_col=0)["AUBOC10"].copy()

        df_o_minus_r = pd.read_excel(BO_data_excel, sheetname="df_o_minus_r", index_col=0)
        df_o_minus_r.columns = pd.Series(df_o_minus_r.columns).replace(namedict)
        df_o_minus_r_mean = df_o_minus_r.T.mean()
        df_o_minus_r_mean_df= pd.concat([df_o_minus_r_mean_df,df_o_minus_r_mean],axis=1, join="outer")
        AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)
        AUBOC10_list.append(AUBOC10)
        linechar_name_list.append(predictor_name)
        AUC_AUBOC_name_list.append("{}-AUC".format(predictor_name))
        AUC_AUBOC_name_list.append("{}-AUBOC10".format(predictor_name))
        thoipapy.figs.Create_Bo_Curve_files.save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict,
                                                 logging, AUC_ser)

        AUC_AUBOC_df = pd.concat([AUC_AUBOC_df,AUC_ser, AUBOC10_ser], axis=1, join="outer")
    #sys.stdout.write(auc_mean_list)
    auc_mean_df = pd.DataFrame.from_records([AUBOC10_list], columns=linechar_name_list)
    auc_mean_df.to_csv(mean_AUBOC_file)
    df_o_minus_r_mean_df.columns = linechar_name_list
    AUC_AUBOC_df.columns = AUC_AUBOC_name_list
    AUC_AUBOC_df.index.name = "acc_db"
    #sys.stdout.write(df_o_minus_r_mean_df, AUBOC10_list)
    #sys.stdout.write(AUC_AUBOC_df)
    AUC_AUBOC_df.to_csv(AUC_AUBOC_file)
    THOIPA_best_set = s["THOIPA_best_set"]
    create_4predictors_AUC_AUBOC10_barchart(AUC_AUBOC_df, predictors_AUC_barchart_png, predictors_BOAUC10_barchart_png, namedict, THOIPA_best_set)
    create_4predictors_bocurve_linechart(df_o_minus_r_mean_df, AUBOC10_list, linechar_name_list, predictors_BOCURVE_linechart_png)
    create_mean_AUC_barchart_comp(auc_mean_list, linechar_name_list, predictors_mean_auc_barchart_png)
    df_o_minus_r_mean_df.to_csv(df_o_minus_r_mean_csv)
    logging.info("finished create_AUC_BoAUC_figs_THOIPA_PREDDIMER_TMDOCK")


def create_mean_AUC_barchart_comp(auc_mean_list,linechar_name_list,predictors_mean_auc_barchart_png):
    plt.close("all")
    # plt.rcParams.update({'font.size': 2})
    mean_auc_name = [linechar_name_list[0],'\n{}\n'.format(linechar_name_list[1]),linechar_name_list[2],'\n{}\n'.format(linechar_name_list[3])]
    figsize = np.array([3.42, 3.42])   # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    # replace the protein names
    x = y_pos = np.arange(len(linechar_name_list))
    plt.bar(x, auc_mean_list, width=0.6, color = 'rgbk', alpha=0.5)
    plt.xticks(y_pos, mean_auc_name,fontsize=6)
    plt.ylabel("performance value\n(mean auc)")

    #ax.set_ylabel("performance value\n(auc)")
    ax.set_ylim(0, 0.70)
    ax.legend()  # (["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(predictors_mean_auc_barchart_png, dpi=240)


def create_4predictors_bocurve_linechart(df_o_minus_r_mean_df,AUBOC10_list,linechar_name_list, predictors_BOCURVE_linechart_png):
    # BO_linechart_png
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    color_list = 'rgbk'
    fig, ax = plt.subplots(figsize=figsize)
    for i,column in enumerate(df_o_minus_r_mean_df.columns):
        # df_o_minus_r_mean_df.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC10 : {:0.2f}".format(AUBOC10))
        label_name = "{}(AUBOC:{:.2f})".format(linechar_name_list[i] ,AUBOC10_list[i])
        df_o_minus_r_mean_df[column].plot(ax=ax,  linestyle="-",label=label_name, color = color_list[i])
    ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="random", alpha=0.5)
    ax.grid(False)
    ax.set_ylabel("performance value\n(observed - random)", color="#0f7d9b")
    ax.tick_params('y', colors="#0f7d9b")

    ax.spines['left'].set_color("#0f7d9b")
    ax.legend()
    fig.tight_layout()
    fig.savefig(predictors_BOCURVE_linechart_png, dpi=140)

def create_4predictors_AUC_AUBOC10_barchart(AUC_AUBOC_df, predictors_AUC_barchart_png, predictors_BOAUC10_barchart_png,namedict, THOIPA_best_set):
    THOIPA_best_setnumber = int(THOIPA_best_set[3:])
    #colname = "THOIPA_5_LOOAUC"
    THOIPA_x_LOOAUC = "THOIPA_{}_LOO-AUC".format(THOIPA_best_setnumber)
    THOIPA_x_LOOAUBOC10 = "THOIPA_{}_LOO-AUBOC10".format(THOIPA_best_setnumber)
    # auc_list = AUC_AUBOC_df.columns[[0, 2, 4, 6]]
    # bo_auc_list = AUC_AUBOC_df.columns[[1, 3, 5, 7]]

    auc_list = AUC_AUBOC_df.columns[::2]
    bo_auc_list = AUC_AUBOC_df.columns[1::2]
    AUC_AUBOC_df = AUC_AUBOC_df.sort_values(by=[THOIPA_x_LOOAUC], ascending=False)
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
    fig.savefig(predictors_AUC_barchart_png, dpi=240)

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
    fig.savefig(predictors_BOAUC10_barchart_png, dpi=240)

def create_AUBOC10_4predictors_3databases_figs(s,df_set,logging):

    logging.info("start create_AUBOC10_4predictors_3databases_figs")
    predictor_name_list = ["THOIPA_{}_LOO".format(s["set_number"]),"PREDDIMER", "TMDOCK", "LIPS_surface_ranked"] #"LIPS_L*E",
    databases = ["crystal", "NMR", "ETRA"]
    for database in databases:
        df_o_minus_r_mean_df = pd.DataFrame()
        AUBOC10_list = []
        mean_AUBOC_file = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_mean_AUBOC10.csv".format(s["setname"],
                                                                                                         database))
        mean_AUBOC_barplot_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                       "{}.{}.4predictors_mean_AUBOC10.png".format(s["setname"],
                                                                                   database))
        predictors_BOCURVE_linechart_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_BOCURVE_linechart.csv".format(s["setname"],
                                                                                                         database))
        predictors_BOCURVE_linechart_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_BOCURVE_linechart.png".format(s["setname"],database))
        for predictor_name in predictor_name_list:
            BO_data_excel = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                         "{}.{}_BO_curve_data.xlsx".format(s['setname'], predictor_name.replace('*', "")))
            df_o_minus_r = pd.read_excel(BO_data_excel,sheetname="df_o_minus_r",index_col=0)
            df_o_minus_r = df_o_minus_r.filter(regex=database, axis=1)
            df_o_minus_r_mean = df_o_minus_r.T.mean()
            AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)
            AUBOC10_list.append(AUBOC10)
            df_o_minus_r_mean_df = pd.concat([df_o_minus_r_mean_df, df_o_minus_r_mean], axis=1, join="outer")

        auboc_mean_df = pd.DataFrame.from_records([AUBOC10_list], columns=predictor_name_list)
        auboc_mean_df.to_csv(mean_AUBOC_file)
        df_o_minus_r_mean_df.columns = predictor_name_list
        df_o_minus_r_mean_df.to_csv(predictors_BOCURVE_linechart_csv)
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
        fig.savefig(predictors_BOCURVE_linechart_png, dpi=140)

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
        mean_auc_list = []
        mean_tpr_list= []
        mean_AUC_file = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_mean_AUC.csv".format(s["setname"],
                                                                                                         database))
        mean_AUC_barplot_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                       "{}.{}.4predictors_mean_AUC.png".format(s["setname"],
                                                                                   database))
        predictors_ROC_curve_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_AUC_ROC.csv".format(s["setname"],
                                                                                                         database))
        predictors_AUC_ROC_png = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                                        "{}.{}.4predictors_AUC_ROC.png".format(s["setname"],database))
        plt.close("all")
        figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
        fig, ax = plt.subplots(figsize=figsize)
        for predictor_name in predictor_name_list:
            mean_auc = []
            mean_tpr = 0.0
            mean_fpr = np.linspace(0, 1, 100)
            n=0
            AUC_data_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_predictors",
                                         "{}.{}_AUC_data.pkl".format(s["setname"],predictor_name.replace('*',"")))
            # open pickle file
            with open(AUC_data_pkl, "rb") as f:
                xv_dict = pickle.load(f)
                for k,v in xv_dict.items():
                    if re.search(database,k):
                        mean_auc.append(v['auc'])
                        fpr = v['fpr']
                        tpr = v['tpr']
                        mean_tpr += interp(mean_fpr, fpr, tpr)
                        mean_tpr[0] = 0.0
                        n = n+1
            mean_tpr /= n
            mean_tpr[-1] = 1.0
            mean_auc = np.mean(mean_auc)
            mean_auc_list.append(mean_auc)
            mean_tpr_list.append(mean_tpr)
            ax.plot(mean_fpr, mean_tpr, lw=1, label="{} (area = {:.2f})".format(predictor_name, mean_auc), alpha=0.8)
        ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
        ax.set_xlim([-0.05, 1.05])
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(predictors_AUC_ROC_png, dpi=240)
        # fig.savefig(crossvalidation_png[:-4] + ".pdf")
        fig.savefig(thoipapy.utils.pdf_subpath(predictors_AUC_ROC_png))

        df_tpr = pd.DataFrame.from_records(list(map(list, zip(*mean_tpr_list))),
                                                    columns=predictor_name_list)
        df_tpr.to_csv(predictors_ROC_curve_csv)

        auc_mean_df = pd.DataFrame.from_records([mean_auc_list], columns=predictor_name_list)
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


def merge_4_files_ALIGNMENT_METHOD(acc, full_seq, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file):
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

    df_train = thoipapy.utils.add_res_num_full_seq_to_df(acc,df_train, df_train_seq, full_seq)
    df_thoipa = thoipapy.utils.add_res_num_full_seq_to_df(acc,df_thoipa, df_thoipa_seq, full_seq)
    df_preddimer = thoipapy.utils.add_res_num_full_seq_to_df(acc,df_preddimer, df_preddimer_seq, full_seq)
    df_tmdock = thoipapy.utils.add_res_num_full_seq_to_df(acc,df_tmdock, df_tmdock_seq, full_seq)

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

def create_ROC_Curve_comp_4predictors(s,df_set,logging):

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
            merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Merged", database,
                                                "{}.merged.csv".format(acc))
            dfm = pd.read_csv(merged_data_csv_path, engine="python", index_col=0)
            dfm.dropna(inplace=True)
            interface = dfm["interface"].values
            if n ==0 or n == 3:
                predict = dfm[predictor].values
            else:
                predict = -1 * dfm[predictor].values
            fpr, tpr, thresholds = roc_curve(interface, predict)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0
        mean_tpr /= len(df_set.index)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        mean_tpr_list.append(mean_tpr)
        ax.plot(mean_fpr, mean_tpr, lw=1,label="{} (area = {:.2f})".format(predictor, mean_auc), alpha=0.8)
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