"""Merge predictions from several predictors into one table.

NOT MAINTAINED. Reachable only from the run_validation pipeline stage, which is switched off in both
shipped settings files and disabled in the functional test. Last exercised for the 2020
publication. Not covered by any test, and excluded from the refactoring applied to the maintained
part of the package.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

import thoipapy
from thoipapy.artefacts import ArtefactPaths


def merge_predictions(paths: ArtefactPaths, df_set, logging, testsetname: str, trainsetname: str):
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
    logging.info("\n--------------- starting merge_predictions ---------------\n")
    # add the THOIPA prediction name to the list of columns to keep
    THOIPA_pred_colname = f"THOIPA_{paths.set_number}_LOO"

    other_predictors_dir = paths.data_dir / "Predictions/other_predictors"

    thoipa_trainsetname = f"thoipa.train{trainsetname}"

    # for simplicity, keep only the predictions. Since the index is unique, it can be added later to the combined file.
    columns_kept_in_combined_file = [
        "residue_num",
        "residue_name",
        THOIPA_pred_colname,
        thoipa_trainsetname,
        "TMDOCK",
        "PREDDIMER",
        "interface",
        "interface_score",
        "LIPS_surface",
        "LIPS_surface_ranked",
        "LIPS_L*E",
        "relative_polarity",
        "conservation",
        "DI4mean",
    ]

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        full_seq = df_set.loc[i, "full_seq"]
        database = df_set.loc[i, "database"]
        # inputs
        train_data_file = paths.combined_features_csv(database, acc)
        THOIPA_LOO_prediction_csv = paths.results_dir / f"predictions/THOIPA_LOO/{database}.{acc}.LOO.prediction.csv"
        PREDDIMER_prediction_file = os.path.join(other_predictors_dir, database, f"{acc}.preddimer.closedist.csv")
        TMDOCK_prediction_file = os.path.join(other_predictors_dir, database, f"{acc}.tmdock.closedist.csv")
        # output
        merged_data_csv_path: Path | str = paths.results_dir / f"predictions/merged/{database}.{acc}.merged.csv"
        thoipapy.utils.make_sure_path_exists(merged_data_csv_path, isfile=True)

        # load the full feature file as the start of dfm
        dfm = pd.read_csv(train_data_file, index_col=0)
        dfm["acc_db_resnum_resname"] = dfm.index
        # set the unique index, based on the residue number in the full sequence
        dfm.set_index("res_num_full_seq", inplace=True)
        # dfm["entropy"] = -1 * dfm["entropy"]
        file_list: list[Path | str] = [THOIPA_LOO_prediction_csv, PREDDIMER_prediction_file, TMDOCK_prediction_file]
        prediction_name_list = [THOIPA_pred_colname, "PREDDIMER", "TMDOCK"]
        if paths.setname == testsetname:
            THOIPA_testset_trainset_csv = (
                paths.results_dir
                / f"predictions/thoipa.train{trainsetname}/{database}.{acc}.thoipa.train{trainsetname}.csv"
            )
            file_list.append(THOIPA_testset_trainset_csv)
            prediction_name_list.append(thoipa_trainsetname)
        n_files_merged = 0
        for n, file in enumerate(file_list):
            prediction_name = prediction_name_list[n]
            if os.path.isfile(file):
                df = pd.read_csv(file, index_col=None)
                assert prediction_name in df.columns.to_list() or "closedist" in df.columns.to_list()
                TMD_seq = df["residue_name"].str.cat()
                if TMD_seq not in full_seq:
                    logging.warning(prediction_name)
                    logging.warning(
                        "Sequence in residue_name column of dataframe is not found in the original df_set sequence."
                        f"\nacc : {acc}\nfile number : {n}\nTMD_seq : {TMD_seq}\nfull_seq in df_set : {full_seq}\n"
                        f"THOIPA_LOO_prediction_csv:{THOIPA_LOO_prediction_csv}\ncsv file:{file}"
                    )
                    if prediction_name in [THOIPA_pred_colname, thoipa_trainsetname]:
                        df = thoipapy.utils.add_mutation_missed_residues_with_na(
                            paths.combined_features_csv(database, acc), acc, database, df
                        )
                        TMD_seq = df["residue_name"].str.cat()
                    # skip protein
                    # continue

                # add the residue number in the full sequence
                df = thoipapy.utils.add_res_num_full_seq_to_df(acc, df, TMD_seq, full_seq, prediction_name, file)

                # set the unique index, based on the residue number in the full sequence
                df.set_index("res_num_full_seq", inplace=True)

                # for structural predictors, rename "closedist" column to the predictor name
                df.columns = pd.Series(df.columns).replace({"closedist": prediction_name})

                # drop all columns but the prediction
                df = df.reindex(index=df.index, columns=[prediction_name])

                assert not df[prediction_name].isnull().values.any()

                # merge the growing dfm file. All rows are included
                dfm = pd.concat([dfm, df], axis=1, join="outer")

                n_files_merged += 1
            else:
                logging.warning(f"Input file not found: {file}")

        # keep the desired columns
        present = set(dfm.columns)
        new_columns_kept_in_combined_file = [c for c in columns_kept_in_combined_file if c in present]
        dfm = dfm[new_columns_kept_in_combined_file]

        # add a completely random "prediction"
        interf_score = dfm.interface_score.dropna().copy()
        dfm.loc[interf_score.index, "random"] = np.random.rand(len(interf_score))
        # NOT RECOMMENDED, AS HEAVY ATOM DIST IS OPPOSITE OF DISRUPTION
        # dfm.loc[interf_score.index, "random"] = shuffle(interf_score.to_numpy())

        # save to "Merged" folder, so as not to get confused with the "combined" files
        dfm.to_csv(merged_data_csv_path)

        logging.info(f"{acc} predictions combined. n_files_merged : {n_files_merged}. ({merged_data_csv_path})")
    logging.info("\n--------------- finished merge_predictions ---------------\n")
