import os
import sys
import pandas as pd
from Bio import pairwise2
import thoipapy


def combine_file_add_PREDDIMER_TMDOCK_THOIPA_prediction(s,columns_kept_in_combined_file):
    set_list = thoipapy.figs.fig_utils.get_set_lists(s)
    for set_number in set_list:
        setname = "set{:02d}".format(int(set_number))
        set_path = thoipapy.common.get_path_of_protein_set(setname, s["set_path"])
        df_set = pd.read_excel(set_path,sheetname="proteins")
        for i in df_set.index:
            acc = df_set.loc[i, "acc"]
            #if acc == "O75460":
            database = df_set.loc[i, "database"]
            train_data_file = os.path.join(s["thoipapy_feature_folder"],"combined",database,"{}.surr20.gaps5.combined_features.csv".format(acc))
            THOIPA_prediction_file = os.path.join(s["thoipapy_folder"],"Predictions","testset_trainset",database,"{}.THOIPA.trainset04.csv".format(acc))
            PREDDIMER_prediction_file = os.path.join(s["PREDDIMER_TMDOCK_folder"],database,"{}.preddimer.closedist.csv".format(acc))
            TMDOCK_prediction_file = os.path.join(s["PREDDIMER_TMDOCK_folder"], database,"{}.tmdock.closedist.csv".format(acc))
            merged_data_xlsx_path = os.path.join(s["combined_floder"],database,"{}.combined.xlsx".format(acc))
            thoipapy.utils.make_sure_path_exists(os.path.dirname(merged_data_xlsx_path))
            merge_4_files(acc,train_data_file,THOIPA_prediction_file,PREDDIMER_prediction_file,TMDOCK_prediction_file,merged_data_xlsx_path,columns_kept_in_combined_file)

def merge_4_files(acc, train_data_file, THOIPA_prediction_file, PREDDIMER_prediction_file,TMDOCK_prediction_file, merged_data_xlsx_path, columns_kept_in_combined_file):
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
    df_thoipa = pd.read_csv(THOIPA_prediction_file,index_col = 0)
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
        sys.stdout.write("\nstarting reindexing of different TM lengths.")
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

        sys.stdout.write("Check the 4 aligned sequences."
                         "\n{}\n".format(dfs.aligned_seq.tolist()))
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


    dfm = pd.concat([df_train, df_thoipa, df_preddimer,df_tmdock], axis=1, join="outer")
    dfm["aa_pos_in_dfm"] = dfm.index.str[1:].astype(int)
    # start numbering at 1 instead of 0
    dfm["aa_pos_in_dfm"] = dfm["aa_pos_in_dfm"] + 1
    dfm.sort_values("aa_pos_in_dfm", inplace=True)
    dfm["Polarity"] = dfm["lipo_Hessa"]
    # use -entropy, named as conservation
    dfm["Conservation"] = -dfm["Entropy"]
    dfm = dfm[columns_kept_in_combined_file]
    with pd.ExcelWriter(merged_data_xlsx_path) as writer:
        df_train.to_excel(writer, sheet_name="df_train")
        df_thoipa.to_excel(writer, sheet_name="df_thoipa")
        df_preddimer.to_excel(writer, sheet_name="df_preddimer")
        df_tmdock.to_excel(writer, sheet_name="df_tmdock")
        dfm.to_excel(writer, sheet_name="dfm")

    sys.stdout.write("\n{} finished. Merged data saved to {}\n\n".format(acc, merged_data_xlsx_path))
    sys.stdout.flush()
