import os

import pandas as pd

import thoipapy


def add_random_interface_to_combined_features_mult_prot(s, df_set, logging):
    """Creates combined data with a randomly chosen interface residues.

    Sorts TMDs according to length.
    Adds the shortest one to the bottom of the list again.
    Iterates through each TMD. The real interface is added as "real_interface"
    The random interface is created by adding the interface for the previous one.
    Saves combined data in a "rand_int" subfolder.

    Only run if generate_randomised_interfaces is in settings excel file.

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
    logging.info("starting add_random_interface_to_combined_features")
    df_set.sort_values("TMD_len", inplace=True)
    df_set.index = range(df_set.shape[0])

    logging.info("df_set.tail(2) before adding protein01 to end of protein list\n{}".format(df_set.tail(2)))

    # add the first TMD again at the end of the list
    df_set.loc[df_set.index.max() + 1, :] = df_set.loc[0, :]

    logging.info("df_set.tail(2) after adding protein01 to end of protein list\n{}".format(df_set.tail(2)))

    df_experiment_data_prev_TMD = pd.DataFrame()
    interface_list_prev_TMD = []

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        feature_combined_file_rand_int = os.path.join(s["thoipapy_data_folder"], "features", "combined", "rand_int", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        thoipapy.utils.make_sure_path_exists(feature_combined_file_rand_int, isfile=True)

        if database == "ETRA":
            experimental_data_file = os.path.join(s["dropbox_dir"], "ETRA_data", "Average_with_interface", "{}_mul_scan_average_data.xlsx".format(acc))
        else:
            #experimental_data_file = os.path.join(s["thoipapy_data_folder"], "features", 'structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))
            experimental_data_file = os.path.join(s["thoipapy_data_folder"], "features", 'structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))

        df_combined = pd.read_csv(feature_combined_file, index_col=0)

        if database == "ETRA":
            df_experiment_data = pd.read_excel(experimental_data_file)
            df_experiment_data = df_experiment_data.rename(columns={"aa_position" : "residue_num", "orig_aa" : "residue_name", "Interface" : "interface", "Disruption" : "interface_score"})
        else:
            df_experiment_data = pd.read_csv(experimental_data_file)
            df_experiment_data = df_experiment_data.rename(columns={"bind": "interface", "closedist": "interface_score"})


        if i > 0:
            # make sure that a range index is used
            df_experiment_data_prev_TMD.index = range(df_experiment_data_prev_TMD.shape[0])
            # reindex so that the length of the interface residue list from the previous TMD now matches that of the current TMD
            df_experiment_data_prev_TMD = df_experiment_data_prev_TMD.reindex(df_combined.index)
            # convert to a list that can be added to the combined file.
            # Replace nan with 0 (non-interface), where the prev TMD is shorter than the current TMD (all except protein 1, when sorted by TMD_len)
            interface_list_prev_TMD = df_experiment_data_prev_TMD.interface.fillna(0).astype(int).tolist()  # [:df_combined.shape[0]]
            """ Instead of adding the full DF of experimental data, simply add a RANDOM list of interface residues(from the previous TMD).
            
            current list of interface residues that belongs to the current TMD
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            list of interface residues extracted from the previous TMD
            [0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            """
            df_combined["interface"] = interface_list_prev_TMD

            # actually add the real interface
            df_combined["real_interface"] = df_experiment_data.interface.tolist()

            # overwrite existing combined features file
            df_combined.to_csv(feature_combined_file_rand_int)
            #logging.info("{} add_random_interface_to_combined_features finished ({})".format(acc, feature_combined_file_rand_int))

        logging.info("{}    real interface for this TMD {}".format(acc, df_experiment_data.interface.tolist()))
        logging.info("{} RANDOM interface from prev TMD {}".format(acc, interface_list_prev_TMD))

        df_experiment_data_prev_TMD = df_experiment_data

    logging.info("add_random_interface_to_combined_features finished")