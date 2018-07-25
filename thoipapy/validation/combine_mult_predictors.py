import os
import pandas as pd
import thoipapy


def merge_predictions(s, df_set, logging):
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
    # add the THOIPA prediction name to the list of columns to keep
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])
    # for simplicity, keep only the predictions. Since the index is unique, it can be added later to the combined file.
    columns_kept_in_combined_file = ['residue_num', 'residue_name', pred_colname, 'TMDOCK', 'PREDDIMER','interface','interface_score',"LIPS_surface","LIPS_surface_ranked", 'LIPS_L*E',"relative_polarity","conservation","DI4mean"]

    #set_list = thoipapy.figs.fig_utils.get_set_lists(s)
    PREDDIMER_TMDOCK_folder = os.path.join(s["base_dir"], "figs", "FigBZ18-PreddimerTmdockComparison")

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        #if acc =="2axtM1":
        full_seq = df_set.loc[i, "full_seq"]
        database = df_set.loc[i, "database"]
        train_data_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,"{}.surr20.gaps5.combined_features.csv".format(acc))
        combined_data_file = os.path.join(s["dropbox_dir"], "THOIPA_data","Features","combined",database,
                                       "{}.surr20.gaps5.combined_features.csv".format(acc))
        thoipapy.utils.make_sure_path_exists(combined_data_file, isfile=True)
        THOIPA_prediction_csv = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.{}.LOO.prediction.csv".format(acc, database, s["setname"]))
        PREDDIMER_prediction_file = os.path.join(PREDDIMER_TMDOCK_folder, database, "{}.preddimer.closedist.csv".format(acc))
        TMDOCK_prediction_file = os.path.join(PREDDIMER_TMDOCK_folder, database, "{}.tmdock.closedist.csv".format(acc))
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))
        #merged_data_xlsx_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.xlsx".format(acc))
        thoipapy.utils.make_sure_path_exists(merged_data_csv_path, isfile=True)
        #merge_4_files_alignment_method_deprecated(acc, full_seq, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file)

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
                        df = thoipapy.utils.add_mutation_missed_residues_with_na(s, acc, database, df)
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
        new_columns_kept_in_combined_file = list(set(columns_kept_in_combined_file).intersection(set(dfm.columns)))
        dfm = dfm[new_columns_kept_in_combined_file]
        # save to "Merged" folder, so as not to get confused with the "combined" files
        dfm.to_csv(merged_data_csv_path)
        logging.info("{} predictions combined. n_files_merged : {}. ({})".format(acc, n_files_merged, merged_data_csv_path))