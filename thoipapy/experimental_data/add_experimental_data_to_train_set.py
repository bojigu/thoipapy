import os
import sys

import pandas as pd

from thoipapy.utils import normalise_between_2_values


def add_experimental_data_to_combined_features_mult_prot(s, df_set, logging):
    """Run add_experimental_data_to_combined_features for a list of proteins.

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
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        if database == "ETRA":
            #experimental_data_file = os.path.join(s["base_dir"], "data_xy", "Figure", "Show_interface", "Interface_xlsx", "{}.xlsx".format(acc))
            experimental_data_file = os.path.join(s["dropbox_dir"], "ETRA_data", "Average_with_interface", "{}_mul_scan_average_data.xlsx".format(acc))
        else:
            #experimental_data_file = os.path.join(s["thoipapy_data_folder"], "features", 'structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))
            experimental_data_file = os.path.join(s["thoipapy_data_folder"], "features", 'structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))

        add_experimental_data_to_combined_features(acc, database, TMD_seq, feature_combined_file, experimental_data_file, logging)


def add_experimental_data_to_combined_features(acc, database, TMD_seq, feature_combined_file, experimental_data_file, logging):
    """Add the "bind" experimental data "interface_score" to the csv file with features.

    Parameters
    ----------
	acc : str
        Protein accession (e.g. UniProt, PDB)
    database : str
        Database name, e.g. "crystal", "NMR" or "ETRA".
    TMD_seq : str
        TMD sequence
    feature_combined_file : str
        Path to csv with all features combined
    experimental_data_file : str
        Path to csv file with the interface_score from experimental data.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    df_combined = pd.read_csv(feature_combined_file, index_col=0)

    if not os.path.isfile(experimental_data_file):
        # try searching for the data files with uppercase accession
        experimental_data_file = experimental_data_file.replace(acc, acc.upper())
        if os.path.isfile(experimental_data_file):
            logging.warning("experimental_data_file IS IN UPPERCASE ({})".format(experimental_data_file))

    if os.path.isfile(experimental_data_file):
        if database == "ETRA":
            sys.stdout.write(experimental_data_file)
            df_experiment_data = pd.read_excel(experimental_data_file, index_col=0)
            # confirm that correct index_col is chosen
            assert list(df_experiment_data.index) == list(range(1, df_experiment_data.shape[0] + 1))
            df_experiment_data = df_experiment_data.rename(columns={"aa_position" : "residue_num", "orig_aa" : "residue_name", "Interface" : "interface", "Disruption" : "interface_score"})
        else:
            df_experiment_data = pd.read_csv(experimental_data_file)
            df_experiment_data = df_experiment_data.rename(columns={"bind": "interface", "closedist": "interface_score"})
            # confirm that correct index_col is chosen
            assert list(df_experiment_data.index) == list(range(df_experiment_data.shape[0]))
            closedist_notes = False
            if closedist_notes:
                min_closedist = df_experiment_data.interface_score.min()
                if min_closedist < 3:
                    sys.stdout.write("\n{} {} min_closedist {:0.2f}".format(acc, database, min_closedist))

        # join the two dataframes together
        # if either the residue_num or residue_name don't match, the rows will be dropped
        df_combined_plus_exp_data = df_experiment_data.merge(df_combined, on=["residue_num", "residue_name"])

        TMD_seq_in_merged_file = df_combined_plus_exp_data.residue_name.str.cat()
        if TMD_seq != TMD_seq_in_merged_file:
            TMD_seq_in_combined_file = df_combined.residue_name.str.cat()
            TMD_seq_in_bind_file = df_experiment_data.residue_name.str.cat()
            sys.stdout.write("\n{}, TMD_seq in protein set   = {}".format(acc, TMD_seq))
            sys.stdout.write("\n{}, TMD_seq_in_combined_file = {}".format(acc, TMD_seq_in_combined_file))
            sys.stdout.write("\n{}, TMD_seq_in_bind_file     = {}".format(acc, TMD_seq_in_bind_file))
            sys.stdout.write("\n{}, TMD_seq_in_merged_file   = {}\n".format(acc, TMD_seq_in_merged_file))
            #sys.stdout.write("TMD_seq in original settings file and final merged features dataframe does not match.")
            raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

        # create normalised interface scores from both ETRA and closedist(NMR/crystal) data
        # 0 = non-interface
        # 0.5 = intermediate
        # 1 = definitely an interface
        if database == "crystal" or database == "NMR":
            # normalize crystal and NMR closedistance to between 0 and 1 with invert, min and max values were set as 2 and 10 angstrom
            df_combined_plus_exp_data["interface_score_norm"] = normalise_between_2_values(df_combined_plus_exp_data["interface_score"], 2, 10, invert=True)
        elif database == "ETRA":
            ###normalize ETRA experimental disruption value to the range of 0 to 1 without invert, the min and max values were set as -0.4 and 0.4
            df_combined_plus_exp_data["interface_score_norm"] = normalise_between_2_values(df_combined_plus_exp_data["interface_score"], -0.4, 0.4)

        new_index: pd.Series = acc + "-" + database + "_" + df_combined_plus_exp_data["residue_num"].apply(lambda x: f"{x:02d}") + df_combined_plus_exp_data["residue_name"]
        df_combined_plus_exp_data.index = new_index

        # overwrite existing combined features file
        df_combined_plus_exp_data.to_csv(feature_combined_file)
        logging.info("{} add_experimental_data_to_combined_features_mult_prot finished ({})".format(acc, experimental_data_file))

    else:
        logging.warning("{} add_experimental_data_to_combined_features failed, {} not found".format(acc, experimental_data_file))