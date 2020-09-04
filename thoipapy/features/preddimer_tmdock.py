import os
from pathlib import Path
from typing import Union

import pandas as pd

from thoipapy.utils import normalise_between_2_values


def add_PREDDIMER_TMDOCK_to_combined_features_mult_prot(s, df_set, logging):
    """Run add_PREDDIMER_TMDOCK_to_combined_features for a list of proteins.

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
        merged_data_csv_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/predictions/merged/{database}.{acc}.merged.csv"

        add_PREDDIMER_TMDOCK_to_combined_features(acc, feature_combined_file, merged_data_csv_path, logging)


def add_PREDDIMER_TMDOCK_to_combined_features(acc, feature_combined_file, merged_data_csv_path, logging):
    """Add the "closedist" predictions from PREDDIMER and TMDOCK to the csv file with features.

    The "MERGED" csv file should contain the output from all available predictors.

    Parameters
    ----------
	acc : str
        Protein accession (e.g. UniProt, PDB)
    feature_combined_file : str
        Path to csv with all features combined
    merged_data_csv_path : str
        Path to csv file with the prediction results from THOIPA, PREDDIMER, TMDOCK, etc
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    df_combined = pd.read_csv(feature_combined_file, index_col=0)
    orig_df_combined_index = df_combined.index
    df_combined.set_index("res_num_full_seq", inplace=True, drop=False)

    if not os.path.isfile(merged_data_csv_path):
        logging.warning("merged_data_csv_path NOT FOUND. TMDOCK and PREDDIMER not added to combined file. ({})".format(merged_data_csv_path))

    if os.path.isfile(merged_data_csv_path):
        df_predictions = pd.read_csv(merged_data_csv_path, index_col=0)
        cols_to_keep = ["PREDDIMER", "TMDOCK"]
        df_predictions = df_predictions.reindex(columns=cols_to_keep, index=df_combined.res_num_full_seq)
        df_predictions.columns = ["PREDDIMER_feature", "TMDOCK_feature"]
        # fill any missing positions with a very high closedist value
        df_predictions.fillna(10, inplace=True)
        for col in ["PREDDIMER_feature", "TMDOCK_feature"]:
            df_predictions[col] = normalise_between_2_values(df_predictions[col], 0.5, 10, invert=True)

        # # fill any missing positions with a very high closedist value
        # df_predictions.fillna(20, inplace=True)
        #
        # for col in cols_to_keep:
        #     if col not in df_predictions.columns:
        #         cols_to_keep.remove(col)
        #     else:
        #         # normalise closedist "10 to 0.5" to "0 to 1"
        #         df_predictions[col] = normalise_between_2_values(df_predictions[col], 0.5, 10, invert=True)
        # df_predictions = df_predictions.reindex(columns=cols_to_keep, index=df_combined.res_num_full_seq)

        # join the two dataframes together
        # if either the residue_num or residue_name don't match, the rows will be dropped
        #df_combined_new = df_predictions.merge(df_combined, on=["residue_num", "residue_name"])
        df_combined_new = pd.concat([df_combined, df_predictions], axis=1, join="outer")

        # TMD_seq_in_new_combined_file = df_combined_new.residue_name.str.cat()
        # if TMD_seq != TMD_seq_in_new_combined_file:
        #     TMD_seq_in_combined_file = df_combined.residue_name.str.cat()
        #     TMD_seq_in_predictions_file = df_predictions.residue_name.str.cat()
        #     sys.stdout.write("\n{}, TMD_seq in protein set   = {}".format(acc, TMD_seq))
        #     sys.stdout.write("\n{}, TMD_seq_in_combined_file = {}".format(acc, TMD_seq_in_combined_file))
        #     sys.stdout.write("\n{}, TMD_seq_in_predictions_file     = {}".format(acc, TMD_seq_in_predictions_file))
        #     sys.stdout.write("\n{}, TMD_seq_in_new_combined_file   = {}\n".format(acc, TMD_seq_in_new_combined_file))
        #     #sys.stdout.write("TMD_seq in original settings file and final merged features dataframe does not match.")
        #     raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

        # revert back to original index. Maybe not necessary?
        df_combined_new.index = orig_df_combined_index
        # overwrite existing combined features file
        df_combined_new.to_csv(feature_combined_file)
        logging.info("{} add_PREDDIMER_TMDOCK_to_combined_features finished ({})".format(acc, merged_data_csv_path))

    else:
        logging.warning("{} add_PREDDIMER_TMDOCK_to_combined_features failed, {} not found".format(acc, merged_data_csv_path))