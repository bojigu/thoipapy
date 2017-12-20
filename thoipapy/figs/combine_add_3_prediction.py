import os
import sys
import pandas as pd
from Bio import pairwise2
import thoipapy


def combine_file_add_PREDDIMER_TMDOCK_THOIPA_prediction(s, df_set, logging):
    # columns_kept_in_combined_file = ['residue_num', 'residue_name', 'conservation', 'lipo_Hessa', 'CoevDImax_norm',
    #                                  'CoevDI4_norm',
    #                                  'CoevDI8_norm', 'CoevMImax_norm', 'CoevMI4_norm', 'CoevMI8_norm',
    #                                  'CumDI4_norm', 'CumDI8_norm', 'CumMI4_norm', 'CumMI8_norm', 'RelPos_TMD',
    #                                  'RelPos_fullseq', 'LIPS_L*E', 'LIPS_surface_ranked', 'Hydrophobicity_sAA',
    #                                  'TMDOCK', 'PREDDIMER']
    # add the THOIPA prediction name to the list of columns to keep
    pred_colname = "THOIPA_{}_LOO".format(s["set_number"])
    # for simplicity, keep only the predictions. Since the index is unique, it can be added later to the combined file.
    columns_kept_in_combined_file = ['residue_num', 'residue_name', pred_colname, 'TMDOCK', 'PREDDIMER']

    #set_list = thoipapy.figs.fig_utils.get_set_lists(s)
    PREDDIMER_TMDOCK_folder = os.path.join(s["base_dir"], "figs", "FigBZ18-PreddimerTmdockComparison")
    #for set_number in set_list:
    #setname = "set{:02d}".format(int(set_number))
    #set_path = thoipapy.common.get_path_of_protein_set(setname, s["sets_folder"])
    #df_set = pd.read_excel(set_path, sheetname="proteins")
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        full_seq = df_set.loc[i, "full_seq"]
        #if acc == "O75460":
        database = df_set.loc[i, "database"]
        train_data_file = os.path.join(s["features_folder"], "combined", database,"{}.surr20.gaps5.combined_features.csv".format(acc))
        #THOIPA_prediction_file = os.path.join(s["thoipapy_data_folder"], "Predictions", "testset_trainset",database, "{}.THOIPA.trainset04.csv".format(acc))
        THOIPA_prediction_file = os.path.join(s["thoipapy_data_folder"], "Predictions", "leave_one_out", database, "{}.{}.LOO.prediction.csv".format(acc, s["setname"]))
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
        file_list = [THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file]
        prediction_name_list = [pred_colname, "PREDDIMER", "TMDOCK"]
        n_files_merged = 0
        for n, file in enumerate(file_list):
            prediction_name = [prediction_name_list[n]]
            if os.path.isfile(file):
                df = pd.read_csv(file, index_col=0)
                seq = df["residue_name"].str.cat()
                if seq not in full_seq:
                    logging.warning("Sequence in residue_name column of dataframe is not found in the original df_set sequence."
                                     "\nacc : {}\nfile number : {}\nTMD_seq : {}\nfull seq in df_set : {}".format(acc, n, seq, full_seq))
                    # skip protein
                    continue

                # add the residue number in the full sequence
                df = thoipapy.utils.add_res_num_full_seq_to_df(acc, df, seq, full_seq)

                if n == 0:
                    print(df.columns)
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
        columns_kept_in_combined_file = list(set(columns_kept_in_combined_file).intersection(set(dfm.columns)))
        dfm = dfm[columns_kept_in_combined_file]
        # save to "Merged" folder, so as not to get confused with the "combined" files
        dfm.to_csv(merged_data_csv_path)
        logging.info("{} predictions combined. n_files_merged : {}. ({})".format(acc, n_files_merged, merged_data_csv_path))

def merge_4_files_ALIGNMENT_METHOD(acc, full_seq, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file, TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file):
    all_files_exist = True
    for path in [train_data_file, THOIPA_prediction_file,PREDDIMER_prediction_file,TMDOCK_prediction_file]:
        if not os.path.isfile(path):
            all_files_exist = False
            print("{} does not exist".format(path))
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

    df_train = thoipapy.utils.add_res_num_full_seq_to_df(df_train, df_train_seq, full_seq)
    df_thoipa = thoipapy.utils.add_res_num_full_seq_to_df(df_thoipa, df_thoipa_seq, full_seq)
    df_preddimer = thoipapy.utils.add_res_num_full_seq_to_df(df_preddimer, df_preddimer_seq, full_seq)
    df_tmdock = thoipapy.utils.add_res_num_full_seq_to_df(df_tmdock, df_tmdock_seq, full_seq)

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
        sys.stdout.write("4 sequences has more than 2 unique sequences, alignment not possible. protein skipped.")
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
