import os

import pandas as pd

import thoipapy


def motifs_from_seq_mult_protein(s, df_set, logging):
    """Runs motifs_from_seq for multiple proteins

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
    logging.info('start parsing lips output to cons and lips scores')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        TMD_seq_pl_surr = df_set.loc[i, "TMD_seq_pl_surr"]
        motifs_file = os.path.join(s["thoipapy_data_folder"], "Features", "motifs", database, "{}.motifs.csv".format(acc))
        thoipapy.utils.make_sure_path_exists(motifs_file, isfile=True)
        tm_surr_left = int(df_set.loc[i, "tm_surr_left"])
        tm_surr_right = int(df_set.loc[i, "tm_surr_right"])
        motifs_from_seq(TMD_seq, TMD_seq_pl_surr, tm_surr_left, tm_surr_right, motifs_file, logging)


def motifs_from_seq(TMD_seq, TMD_seq_pl_surr, tm_surr_left, tm_surr_right, motifs_file, logging):
    """Generate features related to the presence or absence of simple sequence motifs, such as GxxxG.

    Parameters
    ----------
    TMD_seq : str
        TMD sequence
    TMD_seq_pl_surr : str
        TMD sequence plus surrounding residues (e.g. 20 each side)
    tm_surr_left : int
        Number of surrounding residues to the N-terminus (left)
    tm_surr_right : int
        Number of surrounding residues to the C-terminus (right)
    motifs_file : str
        Path to csv containing the features related to sequence motifs
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """

    motif_dict = {"GxxxG" : {"motif_ss" : r"[G].{3}[G]", "motif_len" : 4},
                  "SmxxxSm": {"motif_ss": r"[GASC].{3}[GASC]", "motif_len": 4},
                  "PolarxxxPolar": {"motif_ss": r"[DKERQPHNSGYTWAMC].{3}[DKERQPHNSGYTWAMC]", "motif_len": 4}}
                  # This simple method doesn't work for "LxxLLxL", as the internal residues are not labelled
                  #"LxxLLxL": {"motif_ss": r"([LVI].{2}[LVI][LVI].{1}[LVI])", "motif_len": 6}}

    df_motifs = pd.DataFrame()
    df_motifs["residue_num"] = range(1, len(TMD_seq) + 1)
    df_motifs["residue_name"] = list(TMD_seq)

    for motif_name in motif_dict:
        # motif search string
        motif_ss = motif_dict[motif_name]["motif_ss"]
        # length of the searched segment - 1 (e.g. i, i+4 for GxxxG).
        motif_len = motif_dict[motif_name]["motif_len"]
        list_residues_in_motif = thoipapy.utils.get_list_residues_in_motif(TMD_seq_pl_surr, motif_ss, motif_len)
        #sys.stdout.write(TMD_seq)
        #sys.stdout.write("".join([str(x) for x in list_residues_in_motif])[tm_surr_left:len(TMD_seq_pl_surr) - tm_surr_right])
        # slice out the TMD region
        list_residues_in_motif_TMD_only = list_residues_in_motif[tm_surr_left: len(TMD_seq_pl_surr) - tm_surr_right]
        df_motifs[motif_name] = list_residues_in_motif_TMD_only
    df_motifs.to_csv(motifs_file, index=False)
    logging.info("motifs_from_seq finished ({})".format(motifs_file))