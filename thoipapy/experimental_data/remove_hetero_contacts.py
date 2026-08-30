import os

import pandas as pd

from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import make_sure_path_exists


def remove_crystal_hetero_contact_residues_mult_prot(paths: ArtefactPaths, df_set, logging):
    """Run remove_crystal_hetero_contact_residues_mult_prot for a list of proteins.

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
    hetero_contact_residues_csv = paths.hetero_contact_residues_csv()
    make_sure_path_exists(hetero_contact_residues_csv, isfile=True)

    n_hetero_contact_residues = 0
    hetero_res_summary_dict = {}
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        if database == "crystal":
            feature_combined_file = paths.combined_features_csv(database, acc)
            no_hetero_feature_combined_file = paths.combined_features_csv(database, acc, nohetero=True)
            homo_hetero_contact_file = paths.homohetero_contact_csv(database, acc)

            hetero_contact_num = remove_crystal_hetero_contact_residues(
                acc, feature_combined_file, homo_hetero_contact_file, no_hetero_feature_combined_file, logging
            )
            n_hetero_contact_residues += hetero_contact_num
            hetero_res_summary_dict[acc] = hetero_contact_num

    hetero_ser = pd.Series(hetero_res_summary_dict).T
    hetero_ser.to_csv(hetero_contact_residues_csv)

    logging.info(f"n_hetero_contact_residues: {n_hetero_contact_residues}  ")


def remove_crystal_hetero_contact_residues(
    acc, feature_combined_file, homo_hetero_contact_file, no_hetero_feature_combined_file, logging
):
    """remove the hetero contact residues from the combined csv file with features.

    The "homo_hetero" csv file should contain and mark both the homo and hetero contact residues

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    feature_combined_file : str
        Path to csv with all features combined
    homo_hetero_contact_file : str
        Path to csv file with the homo and hetero interfaces.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    df_combined = pd.read_csv(feature_combined_file, index_col=0)
    # df_addhetreo_combined = pd.DataFrame()
    # df_addhetreo_combined = df_combined
    if not os.path.isfile(homo_hetero_contact_file):
        raise FileNotFoundError(
            f"homo_hetero_contact_file NOT FOUND. hetero contact residues not calculated and added to combined file.\n({homo_hetero_contact_file})"
        )

    hetero_contact_num = 0
    if os.path.isfile(homo_hetero_contact_file):
        df_contacts = pd.read_csv(homo_hetero_contact_file, index_col=0)
        combined_interface = df_combined.interface
        homo_hetero_interface = df_contacts.bind
        hetero_inter = [
            1 if homo_hetero_interface.iloc[i] == 1 and combined_interface.iloc[i] == 0 else 0
            for i in homo_hetero_interface.index
        ]
        hetero_contact_num = hetero_inter.count(1)
        # df_addhetreo_combined["hetero_interface"] = hetero_inter
        # df_addhetreo_combined.to_csv(add_hetero_feature_combined_file)
        hetero_inter_index = []
        for i in range(len(hetero_inter)):
            if hetero_inter[i] == 1:
                hetero_inter_index.append(i)
        df_combined = df_combined.drop(df_combined.index[hetero_inter_index])
        # df_combined["interface"] = hetero_inter
        df_combined.to_csv(no_hetero_feature_combined_file)

        logging.info(f"{acc} add_hetero_contact_to_crystal_combined_files finished ({no_hetero_feature_combined_file})")

    else:
        logging.warning(f"{acc} add_hetero_contact_to_crystal_combined_file failed, {feature_combined_file} not found")
    return hetero_contact_num

    ###################################################################################################
    #                                                                                                 #
    #            combining test data and add physical parameters                                      #
    #                                                                                                 #
    ###################################################################################################
