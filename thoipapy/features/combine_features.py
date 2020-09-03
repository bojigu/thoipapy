import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import thoipapy
from thoipapy import utils as utils
from thoipapy.features.phobius import return_num_tmd
from thoipapy.features.normalise_features import normalise_features

import matplotlib as mpl
# set matplotlib backend to Agg when run on a server
from thoipapy.utils import make_sure_path_exists

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')

def combine_all_features(s, full_seq, acc, database, TMD_seq, TMD_start, feature_combined_file, entropy_file, rate4site_csv, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv, motifs_file, alignment_summary_csv,full_seq_fasta_file,phobius_outfile, logging):
    """Combine all the training features for a particular protein.

    Parameters
    ----------
    s : dict
        Settings dictionary
    full_seq : str
        Full protein sequence
    acc : str
        Protein accession (e.g. UniProt, PDB)
    database : str
        Database name, e.g. "crystal", "NMR" or "ETRA".
    TMD_seq : str
        TMD sequence
    TMD_start : int
        Start of TMD in full sequence
    feature_combined_file : str
        Path to csv with all features combined
    entropy_file : str
        Path to csv file with entropy (conservation) data
    pssm_csv : str
        Path to csv file with the PSSM for the TMD region.
    lipo_csv : str
        Path to csv with the lipophilicity features
    freecontact_parsed_csv : str
        Path to csv with coevolution features
    relative_position_file : str
        Path to csv file with the relative position of each residue in the TMD and full protein sequence
    LIPS_parsed_csv : str
        Path to csv with LIPS output organised into features for machine learning.
    motifs_file : str
        Path to csv containing the features related to sequence motifs
    alignment_summary_csv : str
		Path to csv file containing the summary of the alignments and homologues
		(e.g. how many homologues before and after filtering)
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    thoipapy.utils.make_sure_path_exists(feature_combined_file, isfile=True)

    for n, filepath in enumerate([entropy_file, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv]):
        if not os.path.isfile(filepath):
            raise FileNotFoundError("combine_all_features_mult_prot failed. {} not found".format(filepath))

    entropy_file_df = pd.read_csv(entropy_file)
    pssm_csv_df = pd.read_csv(pssm_csv)
    lipophilicity_file_df = pd.read_csv(lipo_csv)
    # residue number according to full sequence is used as the index here
    freecontact_parsed_csv_df = pd.read_csv(freecontact_parsed_csv, index_col=0)
    relative_position_file_df = pd.read_csv(relative_position_file)
    LIPS_parsed_csv_df = pd.read_csv(LIPS_parsed_csv)
    motifs_df = pd.read_csv(motifs_file)
    rate4site_df = pd.read_csv(rate4site_csv)

    list_of_dfs = [entropy_file_df, pssm_csv_df, lipophilicity_file_df, freecontact_parsed_csv_df,  relative_position_file_df, LIPS_parsed_csv_df, motifs_df]
    for n, df in enumerate(list_of_dfs):
        if True in df.columns.str.contains("Unnamed").tolist():
            raise ValueError("unnamed column found in dataframe number {}".format(n))

    merge1 = entropy_file_df.merge(pssm_csv_df, on=['residue_num','residue_name'])
    merge2 = merge1.merge(lipophilicity_file_df, on=['residue_num','residue_name'])
    merge3 = merge2.merge(freecontact_parsed_csv_df, on=["residue_num","residue_name"])
    merge4 = merge3.merge(relative_position_file_df, on=["residue_num","residue_name"])
    merge5 = merge4.merge(LIPS_parsed_csv_df, on=["residue_num", "residue_name"])
    merge6 = merge5.merge(rate4site_df, on=["residue_num", "residue_name"])
    df_features_single_protein = merge6.merge(motifs_df, on=["residue_num","residue_name"])

    test_indexing = False
    if test_indexing:
        file_list = ["entropy_file", "pssm_csv", "lipo_csv", "freecontact_parsed_csv", "relative_position_file", "LIPS_parsed_csv", "rate4site_csv"]
        df_list = ["entropy_file_df", "pssm_csv_df", "lipophilicity_file_df", "freecontact_parsed_csv_df", "relative_position_file_df", "LIPS_parsed_csv_df", "rate4site_df"]
        sys.stdout.write("{}".format(entropy_file_df))
        #entropy_file_df.loc[0, "residue_name"] = "X"
        #entropy_file_df.loc[0, "Entropy"] = 9999
        #entropy_file_df.loc[10, "Entropy"] = 7777
        #entropy_file_df["residue_num"] = range(3, entropy_file_df.shape[0] + 3)
        for df in [merge1, merge2, merge3, merge4, merge5, merge6, df_features_single_protein]:
            sys.stdout.write(df.shape)
        for n, df in enumerate([entropy_file_df, pssm_csv_df, lipophilicity_file_df, freecontact_parsed_csv_df, relative_position_file_df, LIPS_parsed_csv_df, rate4site_df]):
            dfname = df_list[n]
            TMD_seq_this_df = df.residue_name.str.cat()
            sys.stdout.write("\n{} ({})".format(TMD_seq_this_df, dfname))

    # Raise an error if the TMD sequence does not match original seq in settings file
    TMD_seq_in_merged_file = df_features_single_protein.residue_name.str.cat()
    if TMD_seq != TMD_seq_in_merged_file:
        sys.stdout.write("acc = {}\nTMD_seq in protein set = {}\nmerged                 = {}\n".format(acc, TMD_seq, TMD_seq_in_merged_file))

        sys.stdout.write("\n{}, TMD_seq in protein set   = {}".format(acc, TMD_seq))
        sys.stdout.write("\n{}, TMD_seq_in_merged_file   = {}".format(acc, TMD_seq_in_merged_file))
        sys.stdout.write("\n{}, entropy_file_df          = {}".format(acc, entropy_file_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, rate4site_df             = {}".format(acc, rate4site_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, pssm_csv_df              = {}".format(acc, pssm_csv_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, lipophilicity_file_df    = {}".format(acc, lipophilicity_file_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, freecontact_parsed_csv_df= {}".format(acc, freecontact_parsed_csv_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, relative_position_file_df= {}".format(acc, relative_position_file_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, motifs_df                = {}".format(acc, motifs_df.residue_name.str.cat()))
        raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

    single_prot_aln_result_ser = utils.open_csv_as_series(alignment_summary_csv)

    # n_homologues and n_TMDs are TMD-properties, not a residue properties, and are therefore the same for all residues)

    # add the cube root of the number of homologues
    n_homologues_orig = single_prot_aln_result_ser["n_uniq_TMD_seqs_for_PSSM_FREECONTACT"]
    n_homologues = np.cbrt(n_homologues_orig)
    df_features_single_protein["n_homologues"] = n_homologues

    # add number of TMDs in protein predicted by phobius
    df_features_single_protein["n_TMDs"] = return_num_tmd(s, acc, full_seq, full_seq_fasta_file, phobius_outfile, logging)

    # add the residue number in the full protein sequence
    # this assumes that the index is a range, starting from 0 to x,
    # and that the residue_name exactly matches the original TMD_seq
    df_features_single_protein["res_num_full_seq"] = df_features_single_protein.index + TMD_start

    df_features_single_protein = normalise_features(df_features_single_protein)

    df_features_single_protein.to_csv(feature_combined_file)
    logging.info("{} combine_all_features_mult_prot finished ({})".format(acc, feature_combined_file))


def combine_all_features_mult_prot(s, df_set, logging):
    """Run combine_all_features for all proteins in a list.

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
    logging.info('Combining features into traindata')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        TMD_start = df_set.loc[i, "TMD_start"]
        full_seq = df_set.loc[i, "full_seq"]
        scalename = s["lipophilicity_scale"]
        lipo_csv = os.path.join(s["thoipapy_data_folder"], "Features", "Lipophilicity", database, "{}_{}_lipo.csv".format(acc, scalename))
        relative_position_file = os.path.join(s["thoipapy_data_folder"], "Features", "relative_position", database, "%s.relative_position%s.csv") % (acc, s["surres"])
        LIPS_parsed_csv = os.path.join(s["thoipapy_data_folder"], "Features", "lips_score", database, "{}.surr{}.LIPS_score_parsed.csv".format(acc, s["num_of_sur_residues"]))
        pssm_csv = os.path.join(s["thoipapy_data_folder"], "Features", "pssm", database, "{}.surr{}.gaps{}.pssm.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        entropy_file = os.path.join(s["thoipapy_data_folder"], "Features", "Entropy", database, "{}.surr{}.gaps{}.uniq.entropy.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        rate4site_csv: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}_rate4site.csv")
        freecontact_parsed_csv = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        motifs_file = os.path.join(s["thoipapy_data_folder"], "Features", "motifs", database, "{}.motifs.csv".format(acc))
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        alignment_summary_csv = os.path.join(alignments_dir, "{}.surr{}.gaps{}.alignment_summary.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        full_seq_fasta_file = os.path.join(s["thoipapy_data_folder"], "Proteins", database, "{}.fasta".format(acc))
        phobius_outfile = os.path.join(s["thoipapy_data_folder"], "Proteins", database, "{}.phobius".format(acc))
        combine_all_features(s, full_seq, acc, database, TMD_seq, TMD_start, feature_combined_file, entropy_file, rate4site_csv, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv, motifs_file, alignment_summary_csv, full_seq_fasta_file, phobius_outfile, logging)


def combine_all_train_data_for_machine_learning(s, df_set, logging):
    """Combine training (or test) data for multiple proteins

    Effectively stacks the CSVs on top of each other.

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

    Saved Files
    -----------
    train_data_csv : csv
        csv file with stacked feature data for multiple proteins
        index = range(0, ..)
        columns =
    """
    logging.info('creating train or test data for machine learning')

    train_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/01_train_data_orig.csv"
    make_sure_path_exists(train_data_csv, isfile=True)

    df_set_nonred = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    # set up a dataframe to hold the features for all proteins
    df_all = pd.DataFrame()
    for i in df_set_nonred.index:
        acc = df_set_nonred.loc[i, "acc"]
        database = df_set_nonred.loc[i, "database"]
        if s["remove_crystal_hetero"]:
            if database == "crystal":
                feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                     "{}.nohetero.surr{}.gaps{}.combined_features.csv".format(acc, s[
                                                         "num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
            else:
                feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        else:
            feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                 "{}.surr{}.gaps{}.combined_features.csv".format(acc, s[
                                                     "num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        df_features_new_protein = pd.read_csv(feature_combined_file, index_col=0)
        df_features_new_protein["acc_db"] = "{}-{}".format(acc, database)
#
        # for the first protein, replace the empty dataframe
        if df_all.empty:
            df_all = df_features_new_protein
        else:
            # concatenate the growing dataframe of combined proteins and new dataframe
            df_all = pd.concat([df_all, df_features_new_protein])

    # drop any positions where there is no interface_score (e.g. no mutations, or hetero contacts?)
    if "interface_score" in df_all.columns:
        df_all.dropna(subset=["interface_score"], inplace=True)
    else:
        logging.warning("No experimental data has been added to this dataset!!!")
        raise KeyError("Dataframe is missing 'interface_score' column containing experimentally determined interface residues.")

    # reset the index to be a range (0,...).
    df_all.index = range(df_all.shape[0])

    # reorder the columns
    column_list = ['acc_db', 'interface', 'interface_score', 'residue_num', 'residue_name', 'n_homologues']
    df_all = thoipapy.utils.reorder_dataframe_columns(df_all, column_list)
    new_index: pd.Series = df_all["acc_db"] + "_" + df_all["residue_num"].apply(lambda x: f"{x:02d}") + df_all["residue_name"]
    df_all.index = new_index
    # # remove crystal hetero_interface residues and drop "hetero_interface" column
    # hetero_inter_index = []
    # for i in range(df_all.shape[0]):0
    #     if df_all.loc[i,"hetero_interface"] == 1:
    #         hetero_inter_index.append(i)
    # df_all = df_all.drop(df_all.index[hetero_inter_index])
    # df_all.drop(["hetero_interface"],axis=1, inplace=True)
    df_all.to_csv(train_data_csv)
    logging.info('Finished creating train or test data for machine learning.')