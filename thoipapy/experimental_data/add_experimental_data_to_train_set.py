import os
import sys

import pandas as pd

from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import normalise_between_2_values


def add_experimental_data_to_combined_features_mult_prot(paths: ArtefactPaths, df_set, inter_pair_max: int, logging):
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
        feature_combined_file = paths.combined_features_csv(database, acc)
        if database == "ETRA":
            experimental_data_file = paths.etra_experimental_csv(acc)
        else:
            experimental_data_file = paths.experimental_interface_csv(database, acc, inter_pair_max)

        add_experimental_data_to_combined_features(
            acc, database, TMD_seq, feature_combined_file, experimental_data_file, logging
        )


def add_experimental_data_to_combined_features(
    acc, database, TMD_seq, feature_combined_file, experimental_data_file, logging
):
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
            logging.warning(f"experimental_data_file IS IN UPPERCASE ({experimental_data_file})")

    if os.path.isfile(experimental_data_file):
        if database == "ETRA":
            sys.stdout.write(experimental_data_file)
            df_experiment_data = pd.read_csv(experimental_data_file)
            # The file previously carried an unnamed index column duplicating aa_position, which was
            # read as the index and checked here. The column is gone; aa_position carries the same
            # 1-based residue numbering, so assert on it directly.
            assert list(df_experiment_data.aa_position) == list(range(1, df_experiment_data.shape[0] + 1))
            df_experiment_data = df_experiment_data.rename(
                columns={
                    "aa_position": "residue_num",
                    "orig_aa": "residue_name",
                    "Interface": "interface",
                    "Disruption": "interface_score",
                }
            )
        else:
            df_experiment_data = pd.read_csv(experimental_data_file)
            df_experiment_data = df_experiment_data.rename(
                columns={"bind": "interface", "closedist": "interface_score"}
            )
            # confirm that correct index_col is chosen
            assert list(df_experiment_data.index) == list(range(df_experiment_data.shape[0]))
            closedist_notes = False
            if closedist_notes:
                min_closedist = df_experiment_data.interface_score.min()
                if min_closedist < 3:
                    sys.stdout.write(f"\n{acc} {database} min_closedist {min_closedist:0.2f}")

        # join the two dataframes together
        # if either the residue_num or residue_name don't match, the rows will be dropped
        df_combined_plus_exp_data = df_experiment_data.merge(df_combined, on=["residue_num", "residue_name"])

        TMD_seq_in_merged_file = df_combined_plus_exp_data.residue_name.str.cat()
        if TMD_seq != TMD_seq_in_merged_file:
            TMD_seq_in_combined_file = df_combined.residue_name.str.cat()
            TMD_seq_in_bind_file = df_experiment_data.residue_name.str.cat()
            sys.stdout.write(f"\n{acc}, TMD_seq in protein set   = {TMD_seq}")
            sys.stdout.write(f"\n{acc}, TMD_seq_in_combined_file = {TMD_seq_in_combined_file}")
            sys.stdout.write(f"\n{acc}, TMD_seq_in_bind_file     = {TMD_seq_in_bind_file}")
            sys.stdout.write(f"\n{acc}, TMD_seq_in_merged_file   = {TMD_seq_in_merged_file}\n")
            # sys.stdout.write("TMD_seq in original settings file and final merged features dataframe does not match.")
            raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

        # create normalised interface scores from both ETRA and closedist(NMR/crystal) data
        # 0 = non-interface
        # 0.5 = intermediate
        # 1 = definitely an interface
        if database == "crystal" or database == "NMR":
            # normalize crystal and NMR closedistance to between 0 and 1 with invert, min and max values were set as 2 and 10 angstrom
            df_combined_plus_exp_data["interface_score_norm"] = normalise_between_2_values(
                df_combined_plus_exp_data["interface_score"], 2, 10, invert=True
            )
        elif database == "ETRA":
            ###normalize ETRA experimental disruption value to the range of 0 to 1 without invert, the min and max values were set as -0.4 and 0.4
            df_combined_plus_exp_data["interface_score_norm"] = normalise_between_2_values(
                df_combined_plus_exp_data["interface_score"], -0.4, 0.4
            )

        new_index: pd.Series = (
            acc
            + "-"
            + database
            + "_"
            + df_combined_plus_exp_data["residue_num"].apply(lambda x: f"{x:02d}")
            + df_combined_plus_exp_data["residue_name"]
        )
        df_combined_plus_exp_data.index = new_index

        # overwrite existing combined features file
        df_combined_plus_exp_data.to_csv(feature_combined_file)
        logging.info(f"{acc} add_experimental_data_to_combined_features_mult_prot finished ({experimental_data_file})")

    else:
        logging.warning(f"{acc} add_experimental_data_to_combined_features failed, {experimental_data_file} not found")
