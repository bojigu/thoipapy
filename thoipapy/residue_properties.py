import csv
import math
import os
import platform
import re
import sys
from difflib import SequenceMatcher
from math import isclose
from pathlib import Path
from shutil import copyfile
import numpy as np
import pandas as pd
import scipy as sc
import scipy.stats
import thoipapy
import thoipapy.utils as utils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pandas import Series
from thoipapy.utils import normalise_0_1, normalise_between_2_values
from weighslide import calculate_weighted_windows

import matplotlib as mpl
# set matplotlib backend to Agg when run on a server
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

def calc_lipophilicity(seq, method = "mean"):
    """ Calculates the average hydrophobicity of a sequence according to the Hessa biological scale.

    Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81

    The Hessa scale has been calculated empirically, using the glycosylation assay of TMD insertion.
    Negative values indicate hydrophobic amino acids with favourable membrane insertion.

    Other hydrophobicity scales are in the settings folder. They can be generated as follows.
    hydrophob_scale_path = r"D:\korbinian\korbinian\settings\hydrophobicity_scales.xlsx"
    df_hs = pd.read_excel(hydrophob_scale_path, skiprows=2)
    df_hs.set_index("1aa", inplace=True)
    dict_hs = df_hs.Hessa.to_dict()
    hessa_scale = np.array([value for (key, value) in sorted(dict_hs.items())])
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
     'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
     'W', 'Y']

    Parameters:
    -----------
    seq : string
        Sequence to be analysed. Gaps (-) and unknown amino acids (x) should be ignored.
    method : string
        Method to be used to average the hydrophobicity values over the whole sequence.
        The hydrophobicity score is positive for polar/charged aa, negative for hydrophobic aa.
            "sum" will return the sum of the hydrophobicity scores over the sequence
            "mean" will return the mean of the hydrophobicity scores over the sequence

    Returns:
    --------
    mean hydrophobicity value for the sequence entered

    Usage:
    ------
    from korbinian.utils import calc_lipophilicity
    # for a single sequence
    s = "SAESVGEVYIKSTETGQYLAG"
    calc_lipophilicity(s)
    # for a series of sequences
    TMD_ser = df2.TM01_SW_match_seq.dropna()
    hydro = TMD_ser.apply(lambda x : calc_lipophilicity(x))

    Notes:
    ------
    %timeit results:
    for a 20aa seq: 136 µs per loop
    for a pandas series with 852 tmds: 118 ms per loop
    """
    # hydrophobicity scale
    hessa_scale = np.array([0.11, -0.13, 3.49, 2.68, -0.32, 0.74, 2.06, -0.6, 2.71,
                            -0.55, -0.1, 2.05, 2.23, 2.36, 2.58, 0.84, 0.52, -0.31,
                            0.3, 0.68])
    # convert to biopython analysis object
    analysed_seq = ProteinAnalysis(seq)
    # biopython count_amino_acids returns a dictionary.
    aa_counts_dict = analysed_seq.count_amino_acids()
    # get the number of AA residues used to calculated the hydrophobicity
    # this is not simply the sequence length, as the sequence could include gaps or non-natural AA
    aa_counts_excluding_gaps = np.array(list(aa_counts_dict.values()))
    number_of_residues = aa_counts_excluding_gaps.sum()
    # if there are no residues, don't attempt to calculate a mean. Return np.nan.
    if number_of_residues == 0:
        return np.nan
    # convert dictionary to array, sorted by aa
    aa_counts_arr = np.array([value for (key, value) in sorted(aa_counts_dict.items())])
    multiplied = aa_counts_arr * hessa_scale
    sum_of_multiplied = multiplied.sum()
    if method == "mean":
        return sum_of_multiplied / number_of_residues
    if method == "sum":
        return sum_of_multiplied

def mem_a3m_homologues_filter(s,logging):
    p_r_i_n_t = print
    tmp_list_loc = s["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0][0:5]
        TMD_start = int(row.strip().split(",")[4]) + 1
        TMD_end = int(row.strip().split(",")[4]) + (int(row.strip().split(",")[3]) - int(row.strip().split(",")[2]) + 1)
        #return TMD_end
        homo_a3m_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.parse.a3m") % acc
        if os.path.isfile(homo_a3m_file):

            homo_a3m_file_handle=open(homo_a3m_file,"r")
            homo_filter_file=os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.a3m.mem.uniq.2gaps") %acc
            homo_filter_file_handle = open(homo_filter_file,"w")
            homo_mem_lips_input_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.mem.lips.input") %acc
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            logging.info("starting parsing a3m file: %s\n" %homo_filter_file)
            i = 0
            tm_query=""
            for line in homo_a3m_file_handle:
                tm_str = line[(TMD_start - 1):TMD_end]
                if i == 0:
                    tm_query = tm_str
                    i = i + 1
                    p_r_i_n_t("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                    p_r_i_n_t("{}".format(tm_str), file=homo_filter_file_handle)
                    continue
                mean_hydrophobicity = calc_lipophilicity(tm_str)
                ratio = SequenceMatcher(None, tm_query, tm_str).ratio()
                if not re.search("-", tm_str) and not re.search("X",
                                                                tm_str) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"]and mean_hydrophobicity < s["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    p_r_i_n_t("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                gap_num = tm_str.count("-")
                if (gap_num <= 3 and not re.search("X",
                                                   tm_str) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"] and mean_hydrophobicity < s["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    p_r_i_n_t("{}".format(tm_str), file=homo_filter_file_handle)
                    # homo_filter_file_handle.write(line)
                    # homo_mem_lips_input_file_handle.write(tm_str
            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()
        homo_a3m_file_handle.close()
    tmp_file_handle.close()



def create_PSSM_from_MSA_mult_prot(s, df_set, logging):
    """ Runs create_PSSM_from_MSA for each protein in a list.

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

    Returns
    -------

    """
    logging.info('start pssm calculation')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        TMD_seq_pl_surr5 = df_set.loc[i, "TMD_seq_pl_surr5"]
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir,"{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"],s["max_n_gaps_in_TMD_subject_seq"]))
        pssm_csv = os.path.join(s["thoipapy_data_folder"], "Features", "pssm", database, "{}.surr{}.gaps{}.pssm.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, TMD_seq, logging)

        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.txt".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))

        pssm_csv_surr5 = os.path.join(s["thoipapy_data_folder"], "Features", "pssm", database, "{}.surr5.gaps{}.pssm.csv".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        create_PSSM_from_MSA(path_uniq_TMD_seqs_surr5_for_LIPO, pssm_csv_surr5, acc, TMD_seq_pl_surr5, logging)

def create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, TMD_seq, logging):
    """Creates a PSSM from a multiple sequence alignment in FASTA format.

    Parameters
    ----------
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    pssm_csv : str
        Path to csv file with the PSSM for the TMD region.
    acc : str
        Protein accession (e.g. UniProt, PDB)
    logging : logging.Logger
        Python object with settings for logging to console and file.

    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        thoipapy.utils.make_sure_path_exists(pssm_csv, isfile=True)
        with open(pssm_csv, 'w') as pssm_file_handle:
            mat = []
            writer = csv.writer(pssm_file_handle)
            writer.writerow(['residue_num', 'residue_name', 'A', 'I', 'L', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'M', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P'])
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, "r") as f:
                for line in f.readlines():
                    # strip removes \n. Input currently is not FASTA, no need to check for >
                    mat.append(list(line.strip()))
                    # if not re.search("^>", line):
                    #     mat.append(list(line))

            """ mat is a nested list of all sequences
            [['L', 'G', 'C', 'S', 'A', 'V', 'G', 'G'], ['L', 'A', 'A', 'S', 'A', 'V', 'G', 'P', 'G', 'I', 'G', 'E', 'G'], [.. and so on
            """
            # number of residues in each sequence
            n_residues = len(mat[0])
            if n_residues != len(TMD_seq):
                raise ValueError("Alignment length does not match TMD length. Check that appropriate TMD_seq (or TMD_seq_surr5) has been inserted."
                                 "TMD_seq : {}, first line of alignment = {}".format(TMD_seq, mat[0]))
            # number of sequences in alignment
            n_seqs = len(mat)
            column = []
            # write 20 amino acids as the header of pssm output file
            # pssm_file_handle.write(
            # 'residue'+' '+'A' + ' ' + 'I' + ' ' + 'L' + ' ' + 'V' + ' ' + 'F' + ' ' + 'W' + ' ' + 'Y' + ' ' + 'N' + ' ' + 'C' + ' ' + 'Q' + ' ' + 'M' + ' ' + 'S' + ' ' + 'T' + ' ' + 'D' + ' ' + 'E' + ' ' + 'R' + ' ' + 'H' + ' ' + 'K' + ' ' + 'G' + ' ' + 'P' + '\n')

            #for j in range(0, rowlen - 1):

            # no need to exclude final \n anymore. iterate through number of residues
            for j in range(0, n_residues):
                # iterate through sequences
                for i in range(0, n_seqs):
                    # append the residue in that column
                    # this results in a list of residues at that column position, for all sequences
                    column.append(mat[i][j])
                aa_num = [column.count('A') / n_seqs, column.count('I') / n_seqs, column.count('L') / n_seqs, column.count('V') / n_seqs, column.count('F') / n_seqs,
                          column.count('W') / n_seqs, column.count('Y') / n_seqs, column.count('N') / n_seqs, column.count('C') / n_seqs, column.count('Q') / n_seqs,
                          column.count('M') / n_seqs, column.count('S') / n_seqs, column.count('T') / n_seqs, column.count('D') / n_seqs, column.count('E') / n_seqs,
                          column.count('R') / n_seqs, column.count('H') / n_seqs, column.count('K') / n_seqs, column.count('G') / n_seqs, column.count('P') / n_seqs]
                # add the residue name to the second column
                aa_num.insert(0, TMD_seq[j])
                #aa_num.insert(0, mat[0][j])  #DEPRECATED: Assumes that first sequence is always the original sequence (will break as soon as this is not the case...)
                # add the residue number to the first column
                aa_num.insert(0, j + 1)
                """ the aa_num now looks like this:
                [4, 'S', 0.0, 0.014, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0074, 0.0, 0.0, 0.933, 0.037, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007]
                """
                writer.writerow(aa_num)
                column = []
            logging.info('{} pssm calculation finished ({})'.format(acc, pssm_csv))

    else:
        logging.warning("{} homo_filter_fasta_file does not exist({})".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))


def lipo_from_pssm_mult_prot(s, df_set, logging):
    """Calculates lipophilicity from a PSSM for a list of proteins.

    This function executes lipo_from_pssm for an input list of proteins.

    The calculation is based on the desired hydrophobicity scale. This is currently a fixed object, e.g. "Hessa" in lipo_from_pssm.

    Parameters
    ----------
    s : dict
        settings dictionary
    acc_list : list
        list of accessions (uniprot or pdb)
    database : str
        Database name, e.g. "crystal", "NMR" or "ETRA".
    """
    logging.info('~~~~~~~~~~~~                 starting lipo_from_pssm_mult_prot              ~~~~~~~~~~~~')

    # set name of hydrophobicity scale
    # current options are KyteDoolittle, Wimley, Hessa, Elazar, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, Engelman(GES)
    scalename = s["lipophilicity_scale"]
    failed_acc_list = []
    plot_linechart = True

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]

        tm_surr_left = int(df_set.loc[i, "tm_surr_left"])
        tm_surr_right = int(df_set.loc[i, "tm_surr_right"])
        # reset to 5. AM NOT SURE WHY.
        if tm_surr_left >= 5:
            tm_surr_left = 5
        if tm_surr_right >= 5:
            tm_surr_right = 5

        # pssm_csv = os.path.join(s["thoipapy_data_folder"], "Features", "pssm", database, "{}.mem.2gap.pssm_surr5.csv".format(acc))
        pssm_csv_surr5 = os.path.join(s["thoipapy_data_folder"], "Features", "pssm", database, "{}.surr5.gaps{}.pssm.csv".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        lipo_csv = os.path.join(s["thoipapy_data_folder"], "Features", "Lipophilicity", database, "{}_{}_lipo.csv".format(acc, scalename))

        result_tuple = lipo_from_pssm(acc, pssm_csv_surr5, lipo_csv, tm_surr_left, tm_surr_right, scalename, logging, plot_linechart)

        if result_tuple[1] is False:
            failed_acc_list.append(acc)

    if len(failed_acc_list) > 0:
        sys.stdout.write("\n\nlipo_from_pssm_mult_prot failed for the following : {}".format(failed_acc_list))
    logging.info('~~~~~~~~~~~~                 finished lipo_from_pssm_mult_prot              ~~~~~~~~~~~~')


def lipo_from_pssm(acc, pssm_csv_surr5, lipo_csv, tm_surr_left, tm_surr_right, scalename, logging, plot_linechart=False):
    """Calculates polarity from for a single protein. Residue propensities are taken from a PSSM, generated from a MSA.

    Takes a PSSM as an input. The PSSM should have the fractions of each residue at each position.
    1) Normalises the PSSM so that it adds up to 1.0 (excluding gaps)
    2) Calculates the polarity for each position by multiplying the residue propensities by a hydrophobicity scale
    3) Uses the weighslide module to get the average lipophilicity of the 3 N-terminal aa to each position (e.g. polarity3Nmean)
       and the 3 C-terminal aa to each position (e.g. polarity3Cmean), and the 3 residues centred on the residue of interest (polarity1mean)

    Utilises the excel file containing hydrophobicity scales, within the THOIPApy module.

    Parameters
    ----------
    acc : str
        UniProt or PDB accession.
    pssm_csv_surr5 : str
        Path to the csv containing the PSSM using the alignment of the TMD plus and minus 5 residues
    lipo_csv : str
        Path to csv with the lipophilicity features
    tm_surr_left : int
        Number of surrounding residues to the N-terminus (left)
    tm_surr_right : int
        Number of surrounding residues to the C-terminus (right)
    scalename : str
        Name of hydrophobicity scale. Should match the excel file with the scales.
        Current options are KyteDoolittle, Wimley, Hessa, Elazar, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, Engelman(GES)
    logging : logging.Logger
        Python object with settings for logging to console and file.
    plot_linechart : bool
        Whether to plot the linecharts showing the hydrophobicity.
        Plotting the figure slows down the processing time.

    Saved files
    -----------
    lipo_csv - csv
        csv output file with the polarity values for each position
        index = ["A1", "I2", ....]
        columns = polarity  polarity3Nmean  polarity3Cmean  polarity1mean, etc

    Returns
    -------
    result_tuple : tuple
        (string, bool, string)
        E.g. acc, True, "" if run successfully.
        or if an error occurred,
        acc, False, "pssm_csv not found"
    """
    lipo_excel = lipo_csv[:-4] + ".xlsx"
    lipo_linechart = lipo_csv[:-4] + "_linechart.png"

    if not os.path.isfile(pssm_csv_surr5):
        sys.stdout.write("\n{} skipped for lipo_from_pssm, pssm_csv not found. ".format(acc))
        sys.stdout.write("pssm_csv : {}".format(pssm_csv_surr5))
        sys.stdout.flush()
        return acc, False, "pssm_csv not found"

    if not os.path.isdir(os.path.dirname(lipo_csv)):
        os.makedirs(os.path.dirname(lipo_csv))

    df = pd.read_csv(pssm_csv_surr5, index_col=0)

    """
                residue_name     A     I     L     V    F    W    Y    N    C  ...      M     S     T     D    E    R     H    K     G     P
    residue_num                                                                ...                                                          
    1                      V  0.00  0.06  0.12  0.74  0.0  0.0  0.0  0.0  0.0  ...   0.08  0.00  0.00  0.00  0.0  0.0  0.00  0.0  0.00  0.00
    2                      S  0.00  0.00  0.00  0.00  0.0  0.0  0.0  0.0  0.0  ...   0.00  0.78  0.04  0.02  0.0  0.0  0.02  0.0  0.14  0.00
    3                      P  0.06  0.00  0.00  0.00  0.0  0.0  0.0  0.0  0.0  ...   0.00  0.00  0.00  0.00  0.0  0.0  0.00  0.0  0.02  0.92

    [3 rows x 21 columns]
    """

    df["aa_pos"] = df.index.astype(str)
    df["newindex"] = df.residue_name + df.aa_pos
    df.set_index("newindex", inplace=True, drop=False)

    """
             residue_name     A     I     L     V    F    W    Y    N    C    ...        T     D    E    R     H    K     G     P  aa_pos  newindex
    newindex                                                                  ...                                                                  
    V1                  V  0.00  0.06  0.12  0.74  0.0  0.0  0.0  0.0  0.0    ...     0.00  0.00  0.0  0.0  0.00  0.0  0.00  0.00       1        V1
    S2                  S  0.00  0.00  0.00  0.00  0.0  0.0  0.0  0.0  0.0    ...     0.04  0.02  0.0  0.0  0.02  0.0  0.14  0.00       2        S2
    P3                  P  0.06  0.00  0.00  0.00  0.0  0.0  0.0  0.0  0.0    ...     0.00  0.00  0.0  0.0  0.00  0.0  0.02  0.92       3        P3

    [3 rows x 23 columns]"""

    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    dfa = df[aa_list].copy()
    df["summed"] = dfa.sum(axis=1)
    df["adds_up_to_1.0"] = df["summed"].apply(lambda a: isclose(a, 1.0))

    """For positions that contained gaps, the summed aa propensity does not add to 1.0

             residue_name     A     I     L     V     F    W     Y    N     C       ...           E    R     H    K     G     P  aa_pos  newindex  summed  adds_up_to_1.0
    newindex                                                                        ...                                                                                  
    L30                 L  0.00  0.02  0.84  0.02  0.00  0.0  0.00  0.0  0.04       ...        0.02  0.0  0.00  0.0  0.00  0.00      30       L30    1.00            True
    V31                 V  0.02  0.10  0.00  0.70  0.12  0.0  0.00  0.0  0.00       ...        0.00  0.0  0.00  0.0  0.00  0.00      31       V31    0.98           False
    P32                 P  0.02  0.00  0.00  0.00  0.00  0.0  0.04  0.0  0.00       ...        0.00  0.0  0.02  0.0  0.02  0.76      32       P32    0.98           False

    [3 rows x 25 columns]"""

    # Find all the positions that don't add up to 1.0
    # Normalise all these positions and replace in the df
    # COMMENT OUT THIS SECTION IF THE PSSM IS ALREADY NORMALISED
    positions_that_dont_add_to_1 = df.loc[df["adds_up_to_1.0"] == False].index
    for pos in positions_that_dont_add_to_1:
        dfa.loc[pos, :] = dfa.loc[pos, :] / df.loc[pos, "summed"]
    df["summed_after_norm"] = dfa.sum(axis=1)

    """Now the "summed_after_norm" should add up to 1.0 for each column.

             residue_name     A     I     L     V     F    W     Y    N     C        ...            R     H    K     G     P  aa_pos  newindex  summed  adds_up_to_1.0  summed_after_norm
    newindex                                                                         ...                                                                                                 
    L30                 L  0.00  0.02  0.84  0.02  0.00  0.0  0.00  0.0  0.04        ...          0.0  0.00  0.0  0.00  0.00      30       L30    1.00            True                1.0
    V31                 V  0.02  0.10  0.00  0.70  0.12  0.0  0.00  0.0  0.00        ...          0.0  0.00  0.0  0.00  0.00      31       V31    0.98           False                1.0
    P32                 P  0.02  0.00  0.00  0.00  0.00  0.0  0.04  0.0  0.00        ...          0.0  0.02  0.0  0.02  0.76      32       P32    0.98           False                1.0

    [3 rows x 26 columns]"""

    # DEPRECATED DUE TO PROBLEMS WITH LOCATING SETTINGS FILE IN DOCKER
    # thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    # hydrophob_scale_path = os.path.join(thoipapy_module_path, "setting", "hydrophobicity_scales.xlsx")
    # df_hs = pd.read_excel(hydrophob_scale_path, skiprows=2)
    # df_hs.set_index("1aa", inplace=True)
    # df_hs.sort_index(inplace=True)
    # hs_arr = df_hs[scalename].to_numpy()

    # hard-coded Engelman (GES) hydrophobicity scale
    # if re-implementing flexible scale, use the csv instead
    # "FAULTS: error: can't copy 'setting\hydrophobicity_scales.xlsx': doesn't exist or not a regular file"
    hs_arr = np.array([1.6,2.,-9.2,-8.2,3.7,1.,-3.,3.1,-8.8,2.8,3.4,-4.8,-0.2,-4.1,-12.3,0.6,1.2,2.6,1.9,-0.7])

    """The series should look like this for the Hessa scale. For speed, this is typically converted to a numpy array, sorted alphabetically according to the residue.

    1aa
    A    0.11
    C   -0.13
    D    3.49
    E    2.68
    F   -0.32
    G    0.74
    H    2.06
    I   -0.60
    K    2.71
    L   -0.55
    M   -0.10
    N    2.05
    P    2.23
    Q    2.36
    R    2.58
    S    0.84
    T    0.52
    V   -0.31
    W    0.30
    Y    0.68
    Name: Hessa, dtype: float64
    Out[5]:
    array([ 0.11, -0.13,  3.49,  2.68, -0.32,  0.74,  2.06, -0.6 ,  2.71,
           -0.55, -0.1 ,  2.05,  2.23,  2.36,  2.58,  0.84,  0.52, -0.31,
            0.3 ,  0.68])
    """

    dfh = dfa * hs_arr

    """ dfh is the aa propensity multiplied by the hydrophob scale. It's sum (not mean) is the average hydrophobicity at that position.

                   A    C       D    E    F       G       H      I    K      L      M    N       P    Q    R       S       T       V    W    Y
    newindex                                                                                                                                  
    V1        0.0000 -0.0  0.0000  0.0 -0.0  0.0000  0.0000 -0.036  0.0 -0.066 -0.008  0.0  0.0000  0.0  0.0  0.0000  0.0000 -0.2294  0.0  0.0
    S2        0.0000 -0.0  0.0698  0.0 -0.0  0.1036  0.0412 -0.000  0.0 -0.000 -0.000  0.0  0.0000  0.0  0.0  0.6552  0.0208 -0.0000  0.0  0.0
    P3        0.0066 -0.0  0.0000  0.0 -0.0  0.0148  0.0000 -0.000  0.0 -0.000 -0.000  0.0  2.0516  0.0  0.0  0.0000  0.0000 -0.0000  0.0  0.0
    """

    df_lipo = pd.DataFrame()
    # calculate the mean polarity at each position
    df_lipo["polarity".format(scalename)] = dfh.sum(axis=1)
    # add a range index
    df_lipo.index = range(df_lipo.shape[0])

    # make sure a high polarity value means the residues are polar
    list_scales_polar_high_values = ["Hessa", "Elazar", "Hopp - Woods"]
    list_scales_polar_low_values = ["KyteDoolittle", "Wimley", "Cornette", "Eisenberg", "Rose", "Janin", "Engelman(GES)"]

    if scalename in list_scales_polar_high_values:
        # everything is good in the world. Polar stuff is high.
        pass
    elif scalename in list_scales_polar_low_values :
        # reverse all values so that polar = high
        df_lipo = - df_lipo
    else:
        raise ValueError("Panic. The scalename is wrong somehow.")

    # add the lowest value, so polarity starts at zero and gets larger
    # e.g. lowest value in Hessa scale is -0.6 for Ile
    min_value_dict = {"Hessa" : 0.6, "Elazar" : 1.92, "Hopp - Woods" : 3.4, "KyteDoolittle" : 4.5, "Wimley" : 1.85,
                      "Cornette" : 5.7, "Eisenberg" : 1.38, "Rose" : 0.91, "Janin" : 0.9, "Engelman(GES)" : 3.7}
    lowest_polarity_value = min_value_dict[scalename]
    df_lipo = df_lipo + lowest_polarity_value

    # lowest_polarity_value = 0.6
    # columns = ["polarity", "polarity3Nmean", "polarity3Cmean", "polarity1mean"]
    # for col in columns:
    #     df_lipo[col] = df_lipo[col] + lowest_polarity_value

    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [1, 1, 1, "x", 0, 0, 0]
    polarity_i1_i3_N = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [0, 0, 0, "x", 1, 1, 1]
    polarity_i1_i3_C = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    window = [1, 1, 1]
    polarity1mean = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    # calculate polarity of central position relative to 6 surrounding residues
    window = [1, 1, 1, "x", 1, 1, 1]
    mean_polarity_surr_6_res = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    relative_polarity = df_lipo["polarity"] / mean_polarity_surr_6_res

    # replace positions with nan that could not be properly calculated
    # this inlcludes the first 3 positions of polarity_i1_i3_N, and last 3 positions of polarity_i1_i3_C
    # these nan positions should be OUTSIDE THE TMD
    polarity_i1_i3_N.iloc[0:3] = np.nan
    polarity_i1_i3_C.iloc[-3:] = np.nan

    df_lipo["polarity3Nmean"] = polarity_i1_i3_N
    df_lipo["polarity3Cmean"] = polarity_i1_i3_C
    df_lipo["polarity1mean"] = polarity1mean
    df_lipo["relative_polarity"] = relative_polarity

    """df_lipo now looks like this.

              polarity  polarity3Nmean  polarity3Cmean  polarity1mean
    position                                                             
    0         0.587961               NaN          0.752395       1.051956
    1         1.515952               NaN          0.216839       1.055965
    2         1.063982               NaN          0.039074       0.752395
    3        -0.322749          0.527983          0.214767       0.289119
    4         0.126125          0.376197          0.171906       0.065124
    """

    df_lipo["residue_num"] = df["aa_pos"].values
    df_lipo["residue_name"] = df["residue_name"].values
    #df_lipo.index = range(df_lipo.shape[0])
    #df_lipo.index=[int(i) -tm_surr_left for i in df_lipo.index]
    df_lipo["residue_num"] = [int(i) - tm_surr_left for i in df_lipo["residue_num"]]
    df_lipo.index = df_lipo["residue_num"]

    df_lipo = df_lipo[["residue_name", "polarity", "polarity3Nmean", "polarity3Cmean", "polarity1mean", "relative_polarity"]]
    #df_lipo.set_index("IND", inplace=True)
    if tm_surr_right ==0 :
        df_lipo = df_lipo[tm_surr_left:]
    else:
        df_lipo=df_lipo[tm_surr_left:-tm_surr_right]

    df_lipo.to_csv(lipo_csv)

    # if plot_linechart:
    #     # plot a linechart with lipo, polarity_i1_i3_N, polarity_i1_i3_C
    #     fig, ax = plt.subplots()
    #     try:
    #         df_lipo.plot(ax=ax)
    #         ax.set_ylabel("lipophilicity, {} scale".format(scalename))
    #         ax.set_xticks(range(len(df_lipo)))
    #         ax.set_xticklabels(df_lipo.index, rotation=90)
    #         ax.set_xlabel("")
    #         ax.grid(False)
    #         fig.tight_layout()
    #         fig.savefig(lipo_linechart, dpi=200)
    #         plt.close("all")
    #     except:
    #         logging.info("the acc :{} has no effective alignment".format(acc))

    logging.info("{} lipo_from_pssm_mult_prot finished using {} scale ({})".format(acc, scalename, lipo_csv))
    return acc, True, ""


def entropy_calculation_mult_prot(s, df_set, logging):
    """ Runs entropy_calculation for a set of proteins

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
    logging.info('start entropy calculation')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        #homo_filter_fasta_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m",database,"%s.a3m.mem.uniq.2gaps%s") % (acc,s["surres"])
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        entropy_file = os.path.join(s["thoipapy_data_folder"], "Features", "Entropy", database, "{}.surr{}.gaps{}.uniq.entropy.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, TMD_seq, entropy_file, logging)


def entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, TMD_seq, entropy_file, logging):
    """Calculates conservation of positions using entropy formula.

    S = - sum (Pi + log(Pi))
    See the scipy documentation for more information regarding the entropy formula.
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html

    Note that this code currently considers gaps to be the 21st amino acid.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    entropy_file : str
        Path to csv file with entropy (conservation) data
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):

        if not os.path.isdir(os.path.dirname(entropy_file)):
            os.makedirs(os.path.dirname(entropy_file))

        with open(entropy_file, 'w') as entropy_file_handle:
            mat = []
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
                # iterate through each sequence
                for line in f.readlines():
                    # append each seq as a list to nested list, mat
                    mat.append(list(line.strip()))
                    # if not re.search("^>", line):
                    #     mat.append(list(line))
                n_residues = len(mat[0])
                n_seqs = len(mat)
                column = []
            writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "Entropy"])
            # iterate through each residue position
            for j in range(0, n_residues):
                # iterate through each sequence
                for i in range(0, n_seqs):
                    # add the residue at that column position to list
                    column.append(mat[i][j])
                """
                UP UNTIL NOW, THE CODE MIRRORS EXACTLY WHAT IS DONE FOR THE PSSM CALCULATION
                """

                # convert to series
                column_serie = Series(column)
                # calculates the probabilities
                p_data = column_serie.value_counts() / len(column_serie)
                # input probabilities to get the entropy
                # note that this currently includes gaps "-" as the 21st amino acid
                entropy = sc.stats.entropy(p_data)
                csv_header_for_ncbi_homologues_file = [j + 1, TMD_seq[j], entropy]
                writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_ncbi_homologues_file)
                # entropy_file_handle.write(mat[0][j]+' '+ str(entropy)+'\n')
                column = []
            # entropy_file_handle.close()
            logging.info('{} entropy_calculation finished ({})'.format(acc, entropy_file))
    else:
        logging.warning("{} entropy_calculation failed. {} input file not found".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))


def rate4site_calculation(s, df_set, logging):
    """Calculates conservation of positions using rate4site.

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
    if not "Linux" in platform.system():
        logging.warning("Aborting rate4site calculation, ao is only implemented on linux.")
        return False

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.fas".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        cons_cdhit_input_fasta: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.lipo_seqs_cdhit_input.fas")
        cons_cdhit_output_fasta: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.lipo_seqs_cdhit_output.fas")

        rate4site_input: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.rate4site_input.fas")
        rate4site_orig_output: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.rate4site_orig_output.txt")

        with open(str(cons_cdhit_input_fasta), "w") as f_out:
            with open(path_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
                for line in f_in:
                    f_out.write(line.replace("-",""))

        # delete output file if it exists
        if cons_cdhit_output_fasta.is_file():
            cons_cdhit_output_fasta.unlink()

        len_cdhit_cluster_reps = 1000

        cutoff = 1.0
        cutoff_decrease_per_round = 0.01
        rerun = False

        print(f"decreasing cdhit cutoff for {acc}: ")

        while len_cdhit_cluster_reps > 200:
            if rerun:
                temp = str(cons_cdhit_output_fasta)[:-4] + "temp.fas"
                if Path(temp).is_file():
                    Path(temp).unlink()
                os.rename(cons_cdhit_output_fasta, temp)
                cdhit_cluster_reps = run_cdhit(temp, cons_cdhit_output_fasta, cutoff)
            else:
                cdhit_cluster_reps = run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)

            len_cdhit_cluster_reps = len(cdhit_cluster_reps)
            sys.stdout.write(f"{cutoff:0.2f}({len_cdhit_cluster_reps}), ")
            sys.stdout.flush()
            cutoff -= cutoff_decrease_per_round
            if cutoff <= 0.20:
                to_be_truncated_fasta: str = str(cons_cdhit_output_fasta)[:-4] + "to_be_truncated.fas"
                os.rename(cons_cdhit_output_fasta, to_be_truncated_fasta)
                with open(to_be_truncated_fasta, "r") as f_in:
                    with open(cons_cdhit_output_fasta, "w") as f_out:
                        for n, line in enumerate(f_in):
                            f_out.write(line)
                            if n >= 400:
                                break
                break
            rerun = True

        final_cutoff_used = cutoff + cutoff_decrease_per_round
        sys.stdout.write("\n")
        logging.info(f"cd-hit for rate4site finished. Final cutoff = {final_cutoff_used:0.2f}. Clusters = {len_cdhit_cluster_reps}. Output = {cons_cdhit_output_fasta}")

        #if len(cdhit_cluster_reps) > 200:
        #    cutoff = 0.7
        #    run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)
        #    logging.info(f"repeating cd-hit with a stricter cutoff\nafter cd-hit analysis, {cons_cdhit_output_fasta} has {len(cdhit_cluster_reps)} clusters")


        copy_sequence = False
        with open(path_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
            with open(str(rate4site_input), "w") as f_out:
                for line in f_in:
                    if line[0] == ">":
                        if line[1:] in cdhit_cluster_reps:
                            f_out.write(line)
                            copy_sequence = True
                    else:
                        if copy_sequence:
                            f_out.write(line)
                            copy_sequence = False


        if not rate4site_orig_output.is_file() or (s["rerun_rate4site"] in [True, 1]):

            if os.path.isfile(rate4site_input):

                if not rate4site_orig_output.parent.is_dir():
                    rate4site_orig_output.parent.mkdir(parents=True)

                exect_str = f"rate4site -s {rate4site_input} -o {rate4site_orig_output}"
                print(exect_str)
                command = utils.Command(exect_str)
                command.run(timeout=1200, log_stderr=False)

                if not Path(rate4site_orig_output).is_file():
                    raise FileNotFoundError("rate4site output file is not found")

                # cleanup temp files
                temp_output_files = ["r4s.res", "r4sOrig.res", "TheTree.txt"]
                for temp_output_file in temp_output_files:
                    if Path(temp_output_file).is_file():
                        Path(temp_output_file).unlink()

                logging.info('{} rate4site finished ({})'.format(acc, rate4site_orig_output))

            else:
                logging.warning("{} rate4site failed. {} input file not found".format(acc, path_uniq_TMD_seqs_surr5_for_LIPO))
        else:
            logging.info(f"skipping rate4site algo for existing file {rate4site_orig_output}. Set 'rerun_rate4site' to True to rerun calculation.")

        # convert text output to standard csv
        df = pd.read_csv(rate4site_orig_output, skiprows=range(13), index_col=0, header=None, delim_whitespace=True, error_bad_lines=False, comment="#")
        df.columns = ["seq", "score", "qq-interval", "std", "msa-data"]
        df.to_csv(str(rate4site_orig_output)[:-4] + ".orig.csv")

        # convert standard csv to csv for thoipa features
        df_cons = df.reindex(["seq", "score"], axis=1)
        df_cons.columns = ["residue_name", "rate4site"]
        df_cons.index.name = "residue_num"

        rate4site_csv = rate4site_orig_output.parent / (rate4site_orig_output.name.split(".")[0] + "_rate4site.csv")
        df_cons.to_csv(rate4site_csv)

        logging.info(f"{rate4site_csv} saved")

    logging.info("rate4site_calculation finished")


def get_word_size(cutoff):
    word_size_dict = {
        0.70: 5,
        0.60: 4,
        0.50: 3,
        0.40: 2,
        0.00: 1
    }
    for threshold, word_size in word_size_dict.items():
        if cutoff >= threshold:
            return word_size

    return "error"

def run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff):
    word_size = get_word_size(cutoff)
    word_size_command = "" if cutoff == 1.0 else f"-n {word_size}"
    exect_str = f"cdhit -i {cons_cdhit_input_fasta} -o {cons_cdhit_output_fasta} -c {cutoff:0.2f} {word_size_command}"
    command = utils.Command(exect_str)
    command.run(timeout=120, log_stderr=False)
    assert cons_cdhit_output_fasta.is_file()
    cdhit_cluster_reps = []
    with open(str(cons_cdhit_output_fasta), "r") as f:
        for line in f:
            if line[0] == ">":
                cdhit_cluster_reps.append(line[1:])
    return cdhit_cluster_reps


def coevolution_calculation_with_freecontact_mult_prot(s, df_set, logging):
    """Runs coevolution_calculation_with_freecontact for a set of protein sequences.

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
    logging.info('start coevolution calculation using freecontact')
    freecontact_loc=s["freecontact_dir"]

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_file = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, freecontact_loc, logging)

def coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, freecontact_loc, logging):
    """Runs freecontact from command line on a multiple sequence alignment in FASTA format.

    Parameters
    ----------
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    freecontact_file : str
        Path to freecontact output file.
    freecontact_loc : str
        Path to the executable freecontact file (e.g. "freecontact" or "D:/programs/freecontact/bin/freecontact")
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        try:
            thoipapy.utils.make_sure_path_exists(freecontact_file, isfile=True)
            exect_str = "grep -v '^>' {aln_file} |sed 's/[a-z]//g'|{freecontact} >{freecontact_output_file}".format(
                aln_file=path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact=freecontact_loc, freecontact_output_file=freecontact_file)

            command = utils.Command(exect_str)
            command.run(timeout=400, log_stderr=False)

            logging.info("Output file: %s\n" % freecontact_file)
        except:
            logging.warning("freecontact gives an error")
    else:
        logging.warning("{} does not exist".format(path_uniq_TMD_seqs_for_PSSM_FREECONTACT))

def parse_freecontact_coevolution_mult_prot(s, df_set, logging):
    """Runs parse_freecontact_coevolution on a set of proteins.

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
    logging.info('cumulative co-evolutionary strength parsing')

    for i in df_set.index:
        sys.stdout.write(".")
        sys.stdout.flush()
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_start = int(df_set.loc[i, "TMD_start"])
        TMD_end = int(df_set.loc[i, "TMD_end"])
        freecontact_file = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_parsed_csv = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, TMD_start, TMD_end, logging)
    sys.stdout.write("\n")

def parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, TMD_start, TMD_end, logging):
    """Parses the freecontact output file to create a number of predictive coevolution features

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    freecontact_file : str
        Path to freecontact output file.
    freecontact_parsed_csv : str
        Path to csv with coevolution features
    TMD_start : int
        TMD start in full sequence
    TMD_end : int
        TMD end in full sequence
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    #######################################################################################################
    #                                                                                                     #
    #                     Dictionary method to add initial coevolutionary scores                          #
    #                                                                                                     #
    #######################################################################################################

    if not os.path.isfile(freecontact_file):
        logging.warning("{} parse_freecontact_coevolution failed, {} not found.".format(acc, freecontact_file))
        raise FileNotFoundError("{} parse_freecontact_coevolution failed, {} not found.".format(acc, freecontact_file))
        #return

    dict_di = {}
    dict_mi = {}
    dict_di_list = {}
    dict_mi_list = {}
    s = "-"
    freecontact_parsed_csv_handle = open(freecontact_parsed_csv, 'w')
    freecontact_file_handle = open(freecontact_file, 'r')
    dict_residuenum_residuename = {}
    for row in freecontact_file_handle:  # get the last line of the *.freecontact file which contains the tmd length infor
        residue_pairs = row.strip().split()
        """ residue_pairs
        ['1', 'I', '2', 'T', '0.243618', '0.454792']
        ['1', 'I', '3', 'L', '0.40476', '0.44558']
        NOTE! For the CumDI4 etc, originally the dictionary for mi and di was switched.
        """
        dict_mi[s.join([residue_pairs[0], residue_pairs[2]])] = residue_pairs[4]
        dict_di[s.join([residue_pairs[0], residue_pairs[2]])] = residue_pairs[5]
        dict_residuenum_residuename[int(residue_pairs[0])] = residue_pairs[1]
        dict_residuenum_residuename[int(residue_pairs[2])] = residue_pairs[3]
        if not residue_pairs[0] in dict_di_list:
            dict_di_list[residue_pairs[0]] = []
            dict_di_list[residue_pairs[0]].append(residue_pairs[5])
        else:
            dict_di_list[residue_pairs[0]].append(residue_pairs[5])
        if not residue_pairs[2] in dict_di_list:
            dict_di_list[residue_pairs[2]] = []
            dict_di_list[residue_pairs[2]].append(residue_pairs[5])
        else:
            dict_di_list[residue_pairs[2]].append(residue_pairs[5])
        if not residue_pairs[0] in dict_mi_list:
            dict_mi_list[residue_pairs[0]] = []
            dict_mi_list[residue_pairs[0]].append(residue_pairs[4])
        else:
            dict_mi_list[residue_pairs[0]].append(residue_pairs[4])
        if not residue_pairs[2] in dict_mi_list:
            dict_mi_list[residue_pairs[2]] = []
            dict_mi_list[residue_pairs[2]].append(residue_pairs[4])
        else:
            dict_mi_list[residue_pairs[2]].append(residue_pairs[4])

    """
    dict_di = {'1-2': '0.243618', '1-3': '0.40476', '1-4': '0.0177035', '1-5': '0.106223', '1-6': '0.244482',
    dict_di_list:
    {'1': ['0.454792', '0.44558', '-1.06626', '-0.731704', '-0.252246', '-0.0125942', '-0.0222937', '1.59152', '-0.129083', '1.1289', '-0.242027', '-0.853413', '-1.64731', '-0.745289', '-1.35698', '3.11027', '-1.72332', '-1.01224', '0.598161', '0.603701', '2.23205', '-0.340545'], 
    '2': ['0.454792', '-1.34186', '2.22772', '-0.764966', '0.648499', '-0.844334', '1.09294', '1.74287', '2.06549', '-1.04338', '1.05392', '-2.15485', '-0.468028', '-0.97496', '-0.24502', '-1.48226', '0.550665', '0.913346', '-0.651264', '-2.15379', '0.665787', '0.163698'], 
    '3': ['0.44558', '-1.34186', '-0.972858', '-0.0938075', '-0.319951', '1.78156', '-1.04316', '0.015566', '0.0821186', '0.460809', '-1.37461', '-0.981004', '0.268589', '0.650184', '0.13531', '0.240688', '-0.39947', '2.78247', '1.4023', '-0.697562', '1.62713', '-0.590952'],
    """

    tmd_length = int(row.strip().split()[2])
    CumDI4 = [0] * tmd_length
    CumMI4 = [0] * tmd_length
    CumDI8 = [0] * tmd_length
    CumMI8 = [0] * tmd_length
    coev_all_top4_DI = [0] * tmd_length
    coev_all_top4_MI = [0] * tmd_length
    coev_all_top8_DI = [0] * tmd_length
    coev_all_top8_MI = [0] * tmd_length
    coev_all_max_DI_1 = [0] * tmd_length  # the sum of the top 8
    coev_all_max_MI_1 = [0] * tmd_length
    for key in dict_di_list:
        coev_all_top4_DI[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:4])) / 4
        coev_all_top4_MI[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:4])) / 4
        coev_all_top8_DI[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:8])) / 8
        coev_all_top8_MI[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:8])) / 8
        # MT note: as far as I can see, the coev_all_max_DI is not limited to [0:8], and is therefore the max DI value between
        # the residue of interest, and any of the residues in the TMD
        coev_all_max_DI_1[int(key) - 1] = sorted(dict_di_list[key], reverse=True)[0]
        coev_all_max_MI_1[int(key) - 1] = sorted(dict_mi_list[key], reverse=True)[0]
        # sys.stdout.write(str(key)+"corresponding to"+str(dict_di_list[key]))
    dict_di_value_sort = sorted(dict_di.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]
    dict_mi_value_sort = sorted(dict_mi.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]

    for i in range(0, 4):
        res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
        res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
        CumDI4[res0] = CumDI4[res0] + float(dict_di_value_sort[i][1])
        CumDI4[res1] = CumDI4[res1] + float(dict_di_value_sort[i][1])
        CumMI4[res0] = CumMI4[res0] + float(dict_mi_value_sort[i][1])
        CumMI4[res1] = CumMI4[res1] + float(dict_mi_value_sort[i][1])

    for i in range(0, 8):
        res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
        res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
        CumDI8[res0] = CumDI8[res0] + float(dict_di_value_sort[i][1])
        CumDI8[res1] = CumDI8[res1] + float(dict_di_value_sort[i][1])
        CumMI8[res0] = CumMI8[res0] + float(dict_mi_value_sort[i][1])
        CumMI8[res1] = CumMI8[res1] + float(dict_mi_value_sort[i][1])

    writer = csv.writer(freecontact_parsed_csv_handle, delimiter=',', quotechar='"',
                        lineterminator='\n',
                        quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
    writer.writerow(["residue_num", "residue_name", "coev_all_max_DI", "coev_all_top4_DI", "coev_all_top8_DI", "CumDI4", "CumDI8", "coev_all_max_MI", "coev_all_top4_MI", "coev_all_top8_MI", "CumMI4", "CumMI8"])
    for index in range(len(CumDI8)):
        csv_header_for_cumulative_strength_file = [(index + 1), dict_residuenum_residuename[(index + 1)],
                                                   coev_all_max_DI_1[index], coev_all_top4_DI[index], coev_all_top8_DI[index], CumDI4[index], CumDI8[index], coev_all_max_MI_1[index], coev_all_top4_MI[index],
                                                   coev_all_top8_MI[index], CumMI4[index], CumMI8[index]]
        # writer = csv.writer(freecontact_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
        # quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
        writer.writerow(csv_header_for_cumulative_strength_file)
    freecontact_file_handle.close()
    freecontact_parsed_csv_handle.close()

    # freecontact_parsed_csv_handle.write(str(residue_di[index])+"\t"+str(residue_mi[index])+"\n")

    #######################################################################################################
    #                                                                                                     #
    #                     Pandas method to add further coevolutionary scores                              #
    #                                                                                                     #
    #######################################################################################################

    df = pd.read_csv(freecontact_file, sep=" ", header=None)
    df.columns = ["n1", "res1", "n2", "res2", "MI", "DI"]
    # due to uniprot indexing, residue_num should equal res + TMD_start - 1
    df["n1"] = df.n1 + TMD_start - 1
    df["n2"] = df.n2 + TMD_start - 1
    """df is a csv like this:
           n1 res1  n2 res2        MI        DI
    0  92    I  93    T  0.243618  0.454792
    1  92    I  94    L  0.404760  0.445580
    2  92    I  95    I  0.017704 -1.066260
    3  92    I  96    I  0.106223 -0.731704
    4  92    I  97    F  0.244482 -0.252246
    """

    # open up the existing parsed freecontact file, with CumDI4, CumDI8, etc
    df_out = pd.read_csv(freecontact_parsed_csv)

    # reindex based on the start of the TMD according to the original settings file
    index_range = range(TMD_start, TMD_start + df_out.shape[0])

    df_out.index = index_range

    for XI in ["MI", "DI"]:
        dfp = df.pivot_table(index="n1", columns="n2", values=XI)

        """ asymmetrical pivoted data
        
        Note the padding of 4 residues at all sides of the data, to allow easy indexing.
        
            n2   231  232  233  234  235       236       237       238       239       240 ...        252       253       254       255       256       257  258  259  260  261
        n1                                                                             ...                                                                                 
        231  NaN  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        232  NaN  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        233  NaN  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        234  NaN  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        235  NaN  NaN  NaN  NaN  NaN  0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523  NaN  NaN  NaN  NaN
        236  NaN  NaN  NaN  NaN  NaN       NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577  NaN  NaN  NaN  NaN
        237  NaN  NaN  NaN  NaN  NaN       NaN       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414  NaN  NaN  NaN  NaN
        238  NaN  NaN  NaN  NaN  NaN       NaN       NaN       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822  NaN  NaN  NaN  NaN
        """

        # get full list of residues
        position_list_unique = np.array(list(set(dfp.index.tolist() + dfp.columns.tolist())))
        # padding allows -4 and +4 indexing at ends
        padding = 4
        # DEPRECATED min_ method. Get start and end from df_set
        # get min and max of the actual TMD position (should be equivalent to TMD_start and
        # min_ = int(position_list_unique.min())
        # max_ = int(position_list_unique.max())
        #position_list = range(min_ - padding, max_ + padding + 1)
        position_list = range(TMD_start - padding, TMD_end + padding + 1)
        dfp = dfp.reindex(index = position_list, columns=position_list)
        #put data on both sides of the table for easy indexing
        for col in dfp.columns:
            start = col + 1
            dfp.loc[start:, col] = dfp.loc[col, start:]
        # drop rows with only nan
        #dfp.dropna(how="all", inplace=True)
        """ now is symmetrical, with nans in the central positions
    
            n2   231  232  233  234       235       236       237       238       239       240 ...        252       253       254       255       256       257  258  259  260  261
        n1                                                                                  ...                                                                                 
        231  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        232  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        233  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        234  NaN  NaN  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN ...        NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN  NaN  NaN
        235  NaN  NaN  NaN  NaN       NaN  0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523  NaN  NaN  NaN  NaN
        236  NaN  NaN  NaN  NaN  0.243618       NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577  NaN  NaN  NaN  NaN
        237  NaN  NaN  NaN  NaN  0.404760  0.332451       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414  NaN  NaN  NaN  NaN
        238  NaN  NaN  NaN  NaN  0.017704  0.140595  0.062405       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822  NaN  NaN  NaN  NaN
        """

        # calculate max and mean over all connections (within TMD) for that residue position
        # DEPRECATED. SAME AS ABOVE.
        #df_out["coev_all_max_XI"] = dfp.max(axis=1)
        df_out["coev_all_mean_XI"] = dfp.mean(axis=1)

        for pos in range(TMD_start, TMD_end + 1):
            #iterate from i-1 and i+1 to i-5 and i+5
            for n in range(1, 6):
                i_minus = pos - n if pos - n in dfp.columns else dfp.columns.min()
                i_plus = pos + n if pos + n in dfp.columns else dfp.columns.max()
                # select the two datapoints (e.g. i-4 and i+4 relative to i)
                sel_XI_ser = dfp.loc[pos, [i_minus, i_plus]]
                df_out.loc[pos, "coev_i{}_XI".format(n)] = sel_XI_ser.mean()

                if n == 4:
                    # add the direct connection between i-4 and i+4 (excluding i)
                    sel_XI_ser["m4_to_p4_value"] = dfp.loc[i_plus, i_minus]
                    df_out.loc[pos, "coev_i4_cx_XI"] = sel_XI_ser.mean()

                    # calculate mean and max of all pairwise values between i and from i-4 to i+4 (total of 8 values)
                    m4_to_p4_ser = dfp.loc[pos, i_minus:i_plus]
                    df_out.loc[pos, "coev_i1-i4_XI"] = m4_to_p4_ser.mean()
                    df_out.loc[pos, "coev_i1-i4_max_XI"] = m4_to_p4_ser.max()

        # HEPTAD MOTIF
        a, b, c, d, e, f, g = 1, np.nan, np.nan, 1, 1, np.nan, 1
        # extend the list longer than any TMD
        # e.g. [1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, ......
        hep_list = [a, b, c, d, e, f, g]*10

        highest_XI_face_value = 0

        # iterate through heptad faces, assuming perfect helix
        for face in range(7):
            # truncate heptad [1,nan...] list to match the length of the TMD
            end = dfp.shape[0] + face
            hep_list_trunc = hep_list[face: end]
            hep_list_trunc = np.array(hep_list_trunc)
            # get indices of the residues corresponding to that heptad motif
            #e.g. [ 88  91  92  94  95  98  99 101 102 105 106 108 109 112 113 115 116]
            hep_cols = hep_list_trunc * np.array(dfp.index)
            hep_cols = hep_cols[~np.isnan(hep_cols)].astype(int)

            XI_of_face = dfp.loc[hep_cols, hep_cols]

            """
            XI_of_face is a selection of the full dataframe
            It's still symmetrical. Each pairwise value is there twice.
            The mean shows the average connection between all residues on that face
            
                 88   91        92        94        95        98        99        101       102       105       106       108       109       112       113  115  116                                                                                                                                                    
            88   NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN
            91   NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN
            92   NaN  NaN       NaN  0.404760  0.017704  0.159686  0.042351  0.240405  0.173562  0.031921  0.185876  0.297854  0.132235  0.360217  0.320984  NaN  NaN
            94   NaN  NaN  0.404760       NaN  0.062405  0.451966  0.112373  0.318711  0.266932  0.266313  0.413066  0.347342  0.336857  0.521322  0.538269  NaN  NaN
            95   NaN  NaN  0.017704  0.062405       NaN  0.045973  0.236728  0.014745  0.020983  0.007767  0.093697  0.014614  0.119658  0.114663  0.064796  NaN  NaN
            98   NaN  NaN  0.159686  0.451966  0.045973       NaN  0.096514  0.155630  0.145974  0.083631  0.224182  0.065785  0.314746  0.414121  0.273921  NaN  NaN
            99   NaN  NaN  0.042351  0.112373  0.236728  0.096514       NaN  0.035880  0.032697  0.018841  0.172811  0.035557  0.169880  0.257919  0.108648  NaN  NaN
            101  NaN  NaN  0.240405  0.318711  0.014745  0.155630  0.035880       NaN  0.072537  0.203192  0.172131  0.050251  0.261064  0.371482  0.318819  NaN  NaN
            102  NaN  NaN  0.173562  0.266932  0.020983  0.145974  0.032697  0.072537       NaN  0.037867  0.161862  0.159869  0.227201  0.407755  0.192598  NaN  NaN
            105  NaN  NaN  0.031921  0.266313  0.007767  0.083631  0.018841  0.203192  0.037867       NaN  0.145918  0.156221  0.141419  0.147921  0.111840  NaN  NaN
            106  NaN  NaN  0.185876  0.413066  0.093697  0.224182  0.172811  0.172131  0.161862  0.145918       NaN  0.163811  0.231322  0.573612  0.476921  NaN  NaN
            108  NaN  NaN  0.297854  0.347342  0.014614  0.065785  0.035557  0.050251  0.159869  0.156221  0.163811       NaN  0.203200  0.344927  0.153646  NaN  NaN
            109  NaN  NaN  0.132235  0.336857  0.119658  0.314746  0.169880  0.261064  0.227201  0.141419  0.231322  0.203200       NaN  0.842323  0.438536  NaN  NaN
            112  NaN  NaN  0.360217  0.521322  0.114663  0.414121  0.257919  0.371482  0.407755  0.147921  0.573612  0.344927  0.842323       NaN  0.461507  NaN  NaN
            113  NaN  NaN  0.320984  0.538269  0.064796  0.273921  0.108648  0.318819  0.192598  0.111840  0.476921  0.153646  0.438536  0.461507       NaN  NaN  NaN
            115  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN
            116  NaN  NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN  NaN  NaN
            
            """
            XI_of_face_mean = XI_of_face.mean().mean()
            # keep a record of the "best face"
            if XI_of_face_mean > highest_XI_face_value:
                index_positions_highest_face = hep_cols

        # label the best face as 1, and all other residues as 0
        index_positions_highest_face = set(index_positions_highest_face).intersection(set(df_out.index))
        df_out.loc[index_positions_highest_face, "highest_face_XI"] = 1
        df_out["highest_face_XI"] = df_out["highest_face_XI"].fillna(0).astype(int)

        # replace the XI in the column names with either MI or DI
        new_column_names = pd.Series(df_out.columns).str.replace("XI", XI)
        df_out.columns = new_column_names

    #"""

    # normalise all columns except for residue_num and residue_name
    column_list = ["residue_num", "residue_name"]
    coev_colname_list = df_out.columns.tolist()[len(column_list):]

    # Specifically overwrite normalised values for Cum metrics. Convert to 0 or 1.
    coev_cum_colname_list = ["CumDI4", "CumDI8", "CumMI4", "CumMI8"]
    for col in coev_cum_colname_list:
        df_out[col] = df_out[col].apply(lambda x : 1 if x > 0 else 0)

    # Normalise each. Add as new ColumnName_norm
    for col in coev_colname_list:
        # if "DI" in col:
        #     df_out["{}_norm".format(col)] = normalise_0_1(df_out[col])[0]
        # elif "MI" in col:
        #     df_out["{}_nonnorm".format(col)] = df_out[col]
        #     df_out[col] = normalise_0_1(df_out[col])[0]

        # ASSUME NORMALISED VALUES ARE TO BE USED IN ALL CASES
        df_out["{}_nonnorm".format(col)] = df_out[col]
        df_out[col] = normalise_0_1(df_out[col])[0]
    #"""

    df_out.to_csv(freecontact_parsed_csv)

def calc_relative_position_mult_prot(s, df_set, logging):
    """calculate the residue relative position on the TMD

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
    logging.info('start to calculate the relative positions')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #TMD_start = int(row.strip().split(",")[2])
        #seqlen = int(row.strip().split(",")[1])
        TMD_start = df_set.loc[i, "TMD_start"]
        tm_seq = df_set.loc[i, "full_seq"]
        seqlen = df_set.loc[i, "seqlen"]

        relative_position_file = os.path.join(s["thoipapy_data_folder"], "Features", "relative_position", database, "%s.relative_position%s.csv") % (acc, s["surres"])
        thoipapy.utils.make_sure_path_exists(relative_position_file, isfile=True)
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        calc_relative_position(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, relative_position_file, TMD_start, seqlen, logging)


def calc_relative_position(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, relative_position_file, TMD_start, seqlen, logging):
    """Calculate the residue relative position on the TMD

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    relative_position_file : str
        Path to csv file with the relative position of each residue in the TMD and full protein sequence
    TMD_start : int
        Start of TMD in full sequence
    seqlen : int
        Length of full sequence
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        relative_position_file_handle = open(relative_position_file, 'w')
        mat = []
        writer = csv.writer(relative_position_file_handle, delimiter=',', quotechar='"',
                            lineterminator='\n',
                            quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
        writer.writerow(['residue_num', 'residue_name', 'RelPos_TMD', 'RelPos_fullseq'])

        with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
            mat = []
            for line in f.readlines():
                if not re.search("^>", line):
                    mat.append(list(line))
            tm_seq = mat[0]
            tm_len = len(tm_seq)
            for i in range(1, tm_len):
                rp1 = i / tm_len
                rp2 = (i + TMD_start - 1) / seqlen
                writer.writerow([i, tm_seq[i - 1], rp1, rp2])
        relative_position_file_handle.close()
        logging.info('{} relative position calculation finished ({})'.format(acc, relative_position_file))
        dfrp = pd.read_csv(relative_position_file, index_col=0)
        dfrp["residue_depth"] = (normalise_0_1(1 - abs(dfrp.RelPos_TMD - 0.5))[0]).round(1)
        dfrp.to_csv(relative_position_file)

    else:
        logging.warning("{} calc_relative_position failed, file not found ({})".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))


def LIPS_score_calculation_mult_prot(s, df_set, logging):
    """Run LIPS_score_calculation for a list of proteins.

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

    Returns
    -------

    """
    logging.info('start lips score calculation')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #LIPS_input_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m",database, "%s.mem.lips.input%s") % (acc,s["surres"])
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, s["num_of_sur_residues"]))

        if os.path.isfile(path_uniq_TMD_seqs_no_gaps_for_LIPS):
            # LIPS_output_file = os.path.join(s["thoipapy_data_folder"], "Features", "lips_score", "zpro/NoRedundPro/%s.mem.lips.output") % acc
            #path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, s["num_of_sur_residues"]))

            #LIPS_output_file = os.path.join(s["thoipapy_data_folder"], "Features", "lips_score", database, "%s.mem.lips.output%s") % (acc, s["surres"])

            LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, s["num_of_sur_residues"]))

            LIPS_score_calculation(path_uniq_TMD_seqs_no_gaps_for_LIPS, LIPS_output_file)
        else:
            logging.warning("{} path_uniq_TMD_seqs_no_gaps_for_LIPS not found".format(acc))


def LIPS_score_calculation(input_seq_file, LIPS_output_file):
    """Python version of the LIPS algorithm by Adamian and Liang (2016) Prediction of transmembrane helix orientation in polytopic membrane protenis.

    This script should give exactly the same output as the original perl algorithm.

    Parameters
    ----------
    input_seq_file : str
        Path to text file with a list of sequences
    LIPS_output_file : str
        Path to file with LIPS output result.
    """
    # use different reference to the print function, to aid in finding areas that are being debugged
    p_r_i_n_t = print

    thoipapy.utils.make_sure_path_exists(LIPS_output_file, isfile=True)

    with open(input_seq_file, "r") as file:
        sequence = ' '.join(line.strip() for line in file)

    with open(LIPS_output_file, "w") as LIPS_output_file_handle:
        n = 0
        sump = 0
        sumlip = 0
        sume = {}  # the sum of entropy for each one of the seven surfaces
        sumf = 0
        sumim = {}  # the sum of lipophilicity for each surfaces
        aanum = {}  # the number of residues for each one of seven surfaces
        resnum = 1
        amino = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        ###propi and propm are both from TMLIP scale, propi means the residue lipophilicity propensity in membrane headgroup region
        ##and propm is the residue lipophilicity propensity in hydrophobic core region
        ##in TMLIP scale paper, the membrane headgroup regions is defined as the first (tmlen/5) and (tmlen-tmlen/5) residues which
        ##more likely the membrane bilayer
        ##while the other residues are defined as hydrophobic core region


        propi = {
            'A': 0.71,
            'R': 1.47,
            'N': 0.96,
            'D': 1.20,
            'C': 1.16,
            'Q': 0.61,
            'E': 0.90,
            'G': 0.48,
            'H': 0.82,
            'I': 1.11,
            'L': 1.18,
            'K': 2.38,
            'M': 1.38,
            'F': 1.57,
            'P': 0.99,
            'S': 0.69,
            'T': 0.72,
            'W': 2.45,
            'Y': 1.23,
            'V': 0.98
        }

        propm = dict(
            A=0.82,
            R=0.18,
            N=0.19,
            D=0.29,
            C=1.01,
            Q=0.26,
            E=0.19,
            G=0.35,
            H=0.12,
            I=1.88,
            L=1.71,
            K=0.42,
            M=1.02,
            F=1.97,
            P=0.65,
            S=0.55,
            T=0.66,
            W=1.65,
            Y=0.94,
            V=1.77
        )

        tmp = sequence.split(' ')
        nrow = len(tmp)
        ncol = len(tmp[0])
        bnum = ncol / 5
        oc = {}
        prob = {}
        entropy = {}  ##residue entropy
        exp_entropy = {}  # the final exponential entropy
        lips = {}  ##the lipophilicity score

        for i in range(nrow):
            for j in range(ncol):
                residue = tmp[i][j]
                res_j = ' '.join((residue, str(j)))
                if (res_j in oc.keys()):
                    oc[res_j] = oc[res_j] + 1
                else:
                    oc[res_j] = 1

        for j in range(ncol):
            for res in amino:
                if (' '.join((res, str(j))) in oc):
                    prob[res] = oc[' '.join((res, str(j)))] / nrow
                    if (j in entropy.keys()):
                        entropy[j] = entropy[j] + prob[res] * math.log(prob[res])  # the entropy calculation
                    else:
                        entropy[j] = prob[res] * math.log(prob[res])
                    if ((j <= bnum) or (j > ncol - bnum)):  ###here is the membrane headgroup residues
                        if (j in lips.keys()):
                            lips[j] = lips[j] + prob[res] * propi[res]
                        else:
                            lips[j] = prob[res] * propi[res]
                    else:  ###here is the hydrophobic region residues
                        if (j in lips.keys()):
                            lips[j] = lips[j] + prob[res] * propm[res]
                        else:
                            lips[j] = prob[res] * propm[res]
            exp_entropy[j] = 2.718 ** ((-1) * entropy[j])  # expontional entropy

        for j in sorted(exp_entropy):
            res = tmp[0][j]
            m = resnum + j
            sump = sump + exp_entropy[j]
            sumlip = sumlip + lips[j]

        for i in range(4):  # for the first 4 surfaces
            p_r_i_n_t("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
            # LIPS_output_file_handle.write("SURFACE", "%s" % i, ":")
            j = i
            while j < ncol:
                res = tmp[0][j]
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[j]  # the sum of lipophilicity for surface i
                else:
                    sumim[i] = lips[j]
                prop = lips[j]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[j]  # the sum of entropy for surface i
                else:
                    sume[i] = exp_entropy[j]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1  # the sum of the residue numbers for surface i
                else:
                    aanum[i] = 1
                rn = j + resnum
                # r3=residuename123(res)
                p_r_i_n_t("%3s" % rn, res, "%6.3f" % prop,
                      "%6.3f" % exp_entropy[j],
                      file=LIPS_output_file_handle)  # print residue information which is in surface i
                # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prop,"%6.3f" % exp_entropy[j])
                k = j + 3
                while (k <= j + 4):  # here add the the residues of i+3 and i+4 into surface i to form heptad repeat
                    if (k < ncol):
                        res = tmp[0][k]
                        # r3=residuename123(res)
                        if (i in sumim.keys()):
                            sumim[i] = sumim[i] + lips[k]
                        else:
                            sumim[i] = lips[k]
                        prob = lips[k]
                        if (i in sume.keys()):
                            sume[i] = sume[i] + exp_entropy[k]
                        else:
                            sume[i] = exp_entropy[k]
                        if (i in aanum.keys()):
                            aanum[i] = aanum[i] + 1
                        else:
                            aanum[i] = 1
                        rn = k + resnum
                        p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
                              file=LIPS_output_file_handle)
                        # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
                j = j + 7
        for i in range(4, 7):  # for surfaces from 4 to 6
            p_r_i_n_t("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
            # LIPS_output_file_handle.write("SURFACE", "%s" % i, ":")
            j = i
            while j < ncol:
                res = tmp[0][j]
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[j]
                else:
                    sumim[i] = lips[j]
                prob = lips[j]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[j]
                else:
                    sume[i] = exp_entropy[j]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1
                else:
                    aanum[i] = 1
                rn = j + resnum
                # r3=residuename123(res)
                p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j], file=LIPS_output_file_handle)
                # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j])
                k = j + 3
                while (k <= j + 4):
                    if (k < ncol):
                        res = tmp[0][k]
                        # r3=residuename123(res)
                        if (i in sumim.keys()):
                            sumim[i] = sumim[i] + lips[k]
                        else:
                            sumim[i] = lips[k]
                        prob = lips[k]
                        if (i in sume.keys()):
                            sume[i] = sume[i] + exp_entropy[k]
                        else:
                            sume[i] = exp_entropy[k]
                        if (i in aanum.keys()):
                            aanum[i] = aanum[i] + 1
                        else:
                            aanum[i] = 1
                        rn = k + resnum
                        p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
                              file=LIPS_output_file_handle)
                        # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
                j = j + 7
            k = i - 4
            while (k <= i - 3):  # here adding residues at the first 7 positions
                if (k < ncol):
                    res = tmp[0][k]
                    # r3=residuename123(res)
                    if (i in sumim.keys()):
                        sumim[i] = sumim[i] + lips[k]
                    else:
                        sumim[i] = lips[k]
                    prob = lips[k]
                    if (i in sume.keys()):
                        sume[i] = sume[i] + exp_entropy[k]
                    else:
                        sume[i] = exp_entropy[k]
                    if (i in aanum.keys()):
                        aanum[i] = aanum[i] + 1
                    else:
                        aanum[i] = 1
                    rn = k + resnum
                    p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k], file=LIPS_output_file_handle)
                    # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                k = k + 1
        p_r_i_n_t("SURFACE LIPOPHILICITY ENTROPY   LIPS", file=LIPS_output_file_handle)
        # LIPS_output_file_handle.write("SURFACE LIPOPHILICITY ENTROPY   LIPS")
        for i in sumim.keys():
            avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
            avpim = avpim * 2
            ave = sume[i] / aanum[i]  # average entropy for surface i
            peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
            p_r_i_n_t("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
                  "%8.3f" % peim,
                  file=LIPS_output_file_handle)  # print seven surfaces and see which surface with lowewst LIPS score
            # LIPS_output_file_handle.write("%s" % i, "%10.3f" % avpim, "%8.3f" % ave, "%8.3f" % peim)

def parse_LIPS_score_mult_prot(s, df_set, logging):
    """Runs parse_LIPS_score for a list of sequences.

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
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, s["num_of_sur_residues"]))
        LIPS_parsed_csv = os.path.join(s["thoipapy_data_folder"], "Features", "lips_score", database, "{}.surr{}.LIPS_score_parsed.csv".format(acc, s["num_of_sur_residues"]))
        parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging)

def parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging):
    """Parse the LIPS output to create a CSV with features for input in machine learning algorithm.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    LIPS_output_file : str
        Path to file with LIPS output result.
    LIPS_parsed_csv : str
        Path to csv with LIPS output organised into features for machine learning.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """

    thoipapy.utils.make_sure_path_exists(LIPS_parsed_csv, isfile=True)

    if os.path.isfile(LIPS_output_file):
        # try:
        surface_num = 0
        surface_lips = 100  ##100 is an initialized big number assuming lips score will not bigger than this number
        with open(LIPS_output_file, "r") as LIPS_output_handle:
            with open(LIPS_parsed_csv, "w") as LIPS_parsed_csv_handle:
                i = 0
                array = []
                dict = {}
                for row in LIPS_output_handle:
                    if re.search("^\s+\d+\s+[A-Z]", row):
                        array = row.split()
                        if not int(array[0]) in dict:
                            dict[int(array[0])] = " ".join([array[1], array[2], array[3]])

                    if re.search("^\d{1}\s+", row):
                        surface_num1 = row.split()[0]
                        surface_lips1 = row.split()[3]
                        if (float(surface_lips1) < float(surface_lips)):
                            surface_lips = surface_lips1
                            surface_num = surface_num1
                LIPS_output_handle.close()

                surface_find = 0
                dict1 = {}
                LIPS_output_handle = open(LIPS_output_file, "r")
                for row in LIPS_output_handle:
                    if re.search("^SURFACE\s" + surface_num, row):
                        surface_find = 1
                        continue
                    if surface_find == 1 and re.search("^\s+\d+\s+[A-Z]", row):
                        array = row.split()
                        if not int(array[0]) in dict1:
                            dict1[int(array[0])] = " ".join([array[1], array[2], array[3]])
                    else:
                        surface_find = 0
                LIPS_output_handle.close()

                writer = csv.writer(LIPS_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num", "residue_name", "LIPS_polarity", "LIPS_entropy", "LIPS_surface"])
                for k, v in sorted(dict.items()):
                    v1 = v.split()
                    v1.insert(0, k)
                    if k in dict1:
                        v1.insert(4, 1)
                    else:
                        v1.insert(4, 0)
                    csv_header_for_cons_lips_score_file = v1
                    writer.writerow(csv_header_for_cons_lips_score_file)
                LIPS_parsed_csv_handle.close()
                logging.info('{} lips score parse finished ({})'.format(acc, LIPS_parsed_csv))
    else:
        logging.warning("{} LIPS_output_file not found.")


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


def convert_bind_data_to_csv(s, df_set, logging):
    """Convert bind data (interface-residue or non-interface-residue) to a csv, for all proteins in a list.

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
        bind_file = os.path.join(s["thoipapy_data_folder"], "Features", "Structure", "%s.4.0closedist") % acc
        csv_output_file = os.path.join(s["thoipapy_data_folder"], "Features", "Structure", "%s.4.0closedist.csv") % acc
        if os.path.isfile(bind_file):
            try:
                with open(bind_file,"r") as bind_file_handle:
                    #csv_output_file=os.path.join(s["thoipapy_data_folder"], "Features", "Structure","NoRedundPro/%s.csv") %acc
                    with open(csv_output_file,"w") as csv_output_file_handle:
                        writer = csv.writer(csv_output_file_handle, delimiter=',', lineterminator='\n')
                        writer.writerow(["residue_num", "residue_name", "bind","closedist"])
                        i=1
                        for row in bind_file_handle:
                            array=row.split()
                            if len(array) ==6:
                             csv_header_for_bind_file=[i,array[3].strip('"'),array[5],array[4]]
                             writer.writerow(csv_header_for_bind_file)
                             i=i+1
            except:
                sys.stdout.write("bind file parsing occures errors")
        else:
            sys.stdout.write("{} convert_bind_data_to_csv failed. {} not found".format(acc, bind_file))


###################################################################################################
#                                                                                                 #
#            combining train data and add physical parameters                                     #
#                                                                                                 #
###################################################################################################

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
        freecontact_parsed_csv = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        motifs_file = os.path.join(s["thoipapy_data_folder"], "Features", "motifs", database, "{}.motifs.csv".format(acc))
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        alignment_summary_csv = os.path.join(alignments_dir, "{}.surr{}.gaps{}.alignment_summary.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        full_seq_fasta_file = os.path.join(s["thoipapy_data_folder"], "Proteins", database, "{}.fasta".format(acc))
        phobius_outfile = os.path.join(s["thoipapy_data_folder"], "Proteins", database, "{}.phobius".format(acc))
        combine_all_features(s, full_seq,acc, database, TMD_seq, TMD_start, feature_combined_file, entropy_file, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv, motifs_file, alignment_summary_csv,full_seq_fasta_file,phobius_outfile, logging)




def combine_all_features(s, full_seq, acc, database, TMD_seq, TMD_start, feature_combined_file, entropy_file, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv, motifs_file, alignment_summary_csv,full_seq_fasta_file,phobius_outfile, logging):
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

    list_of_dfs = [entropy_file_df, pssm_csv_df, lipophilicity_file_df, freecontact_parsed_csv_df,  relative_position_file_df, LIPS_parsed_csv_df, motifs_df]
    for n, df in enumerate(list_of_dfs):
        if True in df.columns.str.contains("Unnamed").tolist():
            raise ValueError("unnamed column found in dataframe number {}".format(n))

    merge1 = entropy_file_df.merge(pssm_csv_df, on=['residue_num','residue_name'])
    merge2 = merge1.merge(lipophilicity_file_df, on=['residue_num','residue_name'])
    merge3 = merge2.merge(freecontact_parsed_csv_df, on=["residue_num","residue_name"])
    merge4 = merge3.merge(relative_position_file_df, on=["residue_num","residue_name"])
    merge5 = merge4.merge(LIPS_parsed_csv_df, on=["residue_num", "residue_name"])
    df_features_single_protein = merge5.merge(motifs_df, on=["residue_num","residue_name"])

    test_indexing = False
    if test_indexing:
        file_list = ["entropy_file", "pssm_csv", "lipo_csv", "freecontact_parsed_csv", "relative_position_file", "LIPS_parsed_csv"]
        df_list = ["entropy_file_df", "pssm_csv_df", "lipophilicity_file_df", "freecontact_parsed_csv_df", "relative_position_file_df", "LIPS_parsed_csv_df"]
        sys.stdout.write("{}".format(entropy_file_df))
        #entropy_file_df.loc[0, "residue_name"] = "X"
        #entropy_file_df.loc[0, "Entropy"] = 9999
        #entropy_file_df.loc[10, "Entropy"] = 7777
        #entropy_file_df["residue_num"] = range(3, entropy_file_df.shape[0] + 3)
        for df in [merge1, merge2, merge3, merge4, df_features_single_protein]:
            sys.stdout.write(df.shape)
        for n, df in enumerate([entropy_file_df, pssm_csv_df, lipophilicity_file_df,freecontact_parsed_csv_df,  relative_position_file_df, LIPS_parsed_csv_df]):
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
        sys.stdout.write("\n{}, pssm_csv_df              = {}".format(acc, pssm_csv_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, lipophilicity_file_df    = {}".format(acc, lipophilicity_file_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, freecontact_parsed_csv_df= {}".format(acc, freecontact_parsed_csv_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, relative_position_file_df= {}".format(acc, relative_position_file_df.residue_name.str.cat()))
        sys.stdout.write("\n{}, motifs_df                = {}".format(acc, motifs_df.residue_name.str.cat()))
        raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

    single_prot_aln_result_ser = utils.open_csv_as_series(alignment_summary_csv)
    n_homologues = single_prot_aln_result_ser["n_uniq_TMD_seqs_for_PSSM_FREECONTACT"]
    # add the number of homologues (same number for all residues)
    df_features_single_protein["n_homologues"] = n_homologues
    # add if the protein is multipass or singlepass
    # TEMPORARY! NEEDS TO BE REPLACED WITH A TOPOLOGY PREDICTOR LATER
    df_features_single_protein["n_TMDs"] = return_num_tmd(s, acc,  full_seq, full_seq_fasta_file,phobius_outfile, logging)
    # if database == "crystal":
    #     df_features_single_protein["n_TMDs"] = 2
    # # elif database == "NMR":
    # #     df_features_single_protein["n_TMDs"] = 3
    # else:
    #     df_features_single_protein["n_TMDs"] = 1

    # add the residue number in the full protein sequence
    # this assumes that the index is a range, starting from 0 to x,
    # and that the residue_name exactly matches the original TMD_seq
    df_features_single_protein["res_num_full_seq"] = df_features_single_protein.index + TMD_start

    df_features_single_protein = normalise_features(df_features_single_protein)
    df_features_single_protein = thoipapy.utils.rename_features(df_features_single_protein)

    df_features_single_protein.to_csv(feature_combined_file)
    logging.info("{} combine_all_features_mult_prot finished ({})".format(acc, feature_combined_file))

def return_num_tmd(s, acc,  full_seq, full_seq_fasta_file,phobius_outfile,logging):
    """Calculate the number of TMDs for the protein using the phobius prediction algorithm.

    Important when mixing crystal dataset (multipass) with single-pass protein datasets.
    Gives the extra feature n_TMDs (number of TMDs) in the protein.
    Helps the machine learning technique.

    Requires Phobius (http://phobius.sbc.su.se/data.html). Follow instructions in readme,
    including link creation so "phobius" is recognised in the console.

    Parameters
    ----------
    s : dict
        Settings dictionary
    acc : str
        Protein accession number. (e.g. UniProt acc, PDB_acc + chain letter)
    database : str
        Database name, e.g. "crystal", "NMR" or "ETRA".
    full_seq : str
        Full protein sequence
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    utils.make_sure_path_exists(full_seq_fasta_file, isfile=True)
    with open(full_seq_fasta_file, 'w') as f:
        f.write(">{}\n{}".format(acc, full_seq))
    f.close()
    if "Windows" in platform.system():
        logging.warning("phobius currently not run for Windows! Skipping phobius prediction.")
    else:
        #perl_dir = s["perl_dir"]
        #phobius_dir = s["phobius_dir"]
        #exect_str = "{} {} {}> {}".format(perl_dir, phobius_dir, full_seq_fasta_file, phobius_outfile)
        # use sudo ln -s /path/to/phobius.pl /usr/local/bin/phobius to create a link,
        # so the perl and phobius directory are not necessary
        exect_str = "phobius {}> {}".format(full_seq_fasta_file, phobius_outfile)
        command = utils.Command(exect_str)
        command.run(timeout=400, log_stderr=False)

    if os.path.exists(phobius_outfile):
        tm_num = 1
        with open(phobius_outfile) as file:
            for line in file:
                if re.search('TRANSMEM', line):
                    tm_num = tm_num + 1

        # Options are only 0 (no TMD predicted by phobius), 1, 2, 3, or 4 (4 or more TMDs predicted by phobius)
        if tm_num > 4:
            tm_num = 4
        return tm_num
    else:
        #sys.stdout.write("no phobius output file found, try to check the reason")
        #return None
        raise FileNotFoundError("{} Phobius output not found ({})".format(acc, phobius_outfile))


def normalise_features(df_features_single_protein):
    """Normalise selected features.

    Number of homologues -> normalised between 1 and 8
    conservation -> as minus entropy
    LIPS L*E -> calculated here
    LIPS surface ranked -> calculated based on LIPS_surface and LIPS L*E
    CS, DE, etc -> %C+%S, %D+%E, etc

    Parameters
    ----------
    df_features_single_protein : pd.DataFrame
        Dataframe with all features for a protein
        Index : range index
        Columns : "residue_num", "residue_name", "Entropy", etc
    """

    # normalise number of homologues to 1,2 or 3
    df_features_single_protein["n_homol_norm"] = df_features_single_protein["n_homologues"].apply(normalise_number_of_homologues)
    # convert entropy to conservation by inverting, and adding 3 to give positive values
    # low values are poorly conserved. High values are highly conserved.
    df_features_single_protein["conservation"] = - df_features_single_protein["Entropy"] + 3
    # calculate LIPS L*E for later BO curve, etc
    df_features_single_protein["LIPS_L*E"] = df_features_single_protein.LIPS_polarity * df_features_single_protein.LIPS_entropy
    # rank the LIPS score by adding a fraction of the L*E to the predicted interface (0 or 1)
    df_features_single_protein["LIPS_surface_ranked"] = df_features_single_protein["LIPS_surface"] - (df_features_single_protein["LIPS_L*E"] / 20)
    df_features_single_protein["LIPS_surface_ranked_norm"] = normalise_0_1(df_features_single_protein["LIPS_surface_ranked"])[0]

    # create our own conservation + polarity and conservation*polarity
    df_features_single_protein["cons+polarity"] = df_features_single_protein["conservation"] + df_features_single_protein["polarity"]
    df_features_single_protein["cons*polarity"] = df_features_single_protein["conservation"] * df_features_single_protein["polarity"]

    # add the mean polarity or conservation of positions i, i+4 and i-4
    window = [1,0,0,0,1,0,0,0,1]
    df_features_single_protein["cons4mean"] = calculate_weighted_windows(df_features_single_protein["conservation"], window, statistic="mean", full_output=False)
    df_features_single_protein["polarity4mean"] = calculate_weighted_windows(df_features_single_protein["polarity"], window, statistic="mean", full_output=False)

    df_features_single_protein["CS"] = df_features_single_protein["C"] + df_features_single_protein["S"]
    df_features_single_protein["DE"] = df_features_single_protein["D"] + df_features_single_protein["E"]
    df_features_single_protein["KR"] = df_features_single_protein["K"] + df_features_single_protein["R"]
    df_features_single_protein["QN"] = df_features_single_protein["Q"] + df_features_single_protein["N"]
    df_features_single_protein["LIV"] = df_features_single_protein["L"] + df_features_single_protein["I"] + df_features_single_protein["V"]

    # round the relative positions so there are 10 options (to 10%)
    df_features_single_protein["RelPos_TMD"] = df_features_single_protein["RelPos_TMD"].round(1)
    df_features_single_protein["RelPos_fullseq"] = df_features_single_protein["RelPos_fullseq"].round(1)

    return df_features_single_protein



def normalise_number_of_homologues(x):
    """Convert non-linear number of homologues to an integer value.

    Parameters
    ----------
    x : int
        Number of homologues in fasta alignment of TMD region including gaps.
    """
    if x <= 75:
        return 1
    elif x <= 100:
        return 2
    elif x <= 200:
        return 3
    elif x <= 400:
        return 4
    elif x <= 800:
        return 5
    elif x <= 1600:
        return 6
    elif x <= 3200:
        return 7
    else:
        return 8


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
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        if database == "ETRA":
            #experimental_data_file = os.path.join(s["base_dir"], "data_xy", "Figure", "Show_interface", "Interface_xlsx", "{}.xlsx".format(acc))
            experimental_data_file = os.path.join(s["dropbox_dir"], "ETRA_data", "Average_with_interface", "{}_mul_scan_average_data.xlsx".format(acc))
        else:
            #experimental_data_file = os.path.join(s["thoipapy_data_folder"], "Features", 'Structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))
            experimental_data_file = os.path.join(s["thoipapy_data_folder"], "Features", 'Structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))

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

        # overwrite existing combined features file
        df_combined_plus_exp_data.to_csv(feature_combined_file)
        logging.info("{} add_experimental_data_to_combined_features_mult_prot finished ({})".format(acc, experimental_data_file))

    else:
        logging.warning("{} add_experimental_data_to_combined_features failed, {} not found".format(acc, experimental_data_file))


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
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        feature_combined_file_rand_int = os.path.join(s["thoipapy_data_folder"], "Features", "combined", "rand_int", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        thoipapy.utils.make_sure_path_exists(feature_combined_file_rand_int, isfile=True)

        if database == "ETRA":
            #experimental_data_file = os.path.join(s["base_dir"], "data_xy", "Figure", "Show_interface", "Interface_xlsx", "{}.xlsx".format(acc))
            experimental_data_file = os.path.join(s["dropbox_dir"], "ETRA_data", "Average_with_interface", "{}_mul_scan_average_data.xlsx".format(acc))
        else:
            #experimental_data_file = os.path.join(s["thoipapy_data_folder"], "Features", 'Structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))
            experimental_data_file = os.path.join(s["thoipapy_data_folder"], "Features", 'Structure', database, '{}.{}pairmax.bind.closedist.csv'.format(acc,s['inter_pair_max']))

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
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))

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


def add_physical_parameters_to_features_mult_prot(s, df_set, logging):
    """Run add_physical_parameters_to_features for multiple proteins.

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
    logging.info('adding physical parameters into traindata')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        #feature_combined_file_incl_phys_param = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
        #                                                     "{}.surr{}.gaps{}.combined_features_incl_phys_param.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        add_physical_parameters_to_features(acc, feature_combined_file, logging)


def add_physical_parameters_to_features(acc, feature_combined_file, logging):
    """Add physical parameters (e.g. PI of amino acid) to the features.

    Parameters
    ----------
	acc : str
        Protein accession (e.g. UniProt, PDB)
    feature_combined_file : str
        Path to csv with all features combined
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    feature_combined_file_incl_phys_param = feature_combined_file[:-4] + "_incl_phys_param.csv"
    thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    physical_parameter_file = os.path.join(thoipapy_module_path, "setting", "Physical_Property_csv.txt")

    dictionary = {}

    with open(physical_parameter_file, "r") as physical_parameter_file_handle:

        if os.path.isfile(feature_combined_file):
            with open(feature_combined_file, "r") as train_data_file_handle:
                # train_data_add_physical_parameter_file = os.path.join("/scratch/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.physipara.traindata1.csv") %acc
                train_data_add_physical_parameter_file_handle = open(feature_combined_file_incl_phys_param, "w")
                # train_data_physical_parameter_file_handle = open(r"/scratch/zeng/thoipapy/Features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
                writer = csv.writer(train_data_add_physical_parameter_file_handle, delimiter=',', quotechar='"',
                                    lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                for row in physical_parameter_file_handle:
                    if re.search("^Hydro", row):
                        continue
                    else:
                        array = row.split()
                        dictionary[array[0]] = array[1:16]
                for row1 in train_data_file_handle:
                    if re.search("residue_num", row1):
                        array1 = row1.rstrip().split(",")
                        array1[42:14] = ["Hydrophobicity_sAA", "Charge_sAA", "PI_sAA", "LIPSI_sAA", "LIPSM_sAA", "Hydrophobic_sAA", "Aliphatic_sAA", "Aromatic_sAA", "Polar_sAA", "Negative_sAA", "Positive_sAA", "Small_sAA", "branched",
                                         "mass", "Volume_sAA"]#Mass_sAA
                        # array2 = array1[0:31]
                        # array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                        writer.writerow(array1)
                    else:
                        array3 = row1.rstrip().split(",")
                        array3[42:14] = dictionary[array3[2]]
                        writer.writerow(array3)
                train_data_add_physical_parameter_file_handle.close()

        else:
            logging.warning("{} does not exist".format(feature_combined_file))

    df = pd.read_csv(feature_combined_file_incl_phys_param, index_col=0)
    cols_to_delete = ["Charge_sAA", "LIPSI_sAA", "LIPSM_sAA", "Hydrophobic_sAA", "Aliphatic_sAA", "Negative_sAA", "Positive_sAA", "Volume_sAA"]
    df.drop(cols_to_delete, axis=1, inplace=True)
    df.to_csv(feature_combined_file_incl_phys_param)

    # overwrite the original feature_combined_file, and delete feature_combined_file_incl_phys_param
    copyfile(feature_combined_file_incl_phys_param, feature_combined_file)
    try:
        os.remove(feature_combined_file_incl_phys_param)
    except:
        logging.warning("{} could not be removed".format(feature_combined_file_incl_phys_param))

    logging.info("{} add_physical_parameters_to_features_mult_prot finished. (updated {})".format(acc, feature_combined_file))

def remove_crystal_hetero_contact_residues_mult_prot(s, df_set, logging):
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
    sum_hetero_contact = 0
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        if database == "crystal":
            feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                             "{}.surr{}.gaps{}.combined_features.csv".format(acc,
                                                                                             s["num_of_sur_residues"],
                                                                                             s[
                                                                                                 "max_n_gaps_in_TMD_subject_seq"]))
            no_hetero_feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                 "{}.nohetero.surr{}.gaps{}.combined_features.csv".format(acc,
                                                                                                 s[
                                                                                                     "num_of_sur_residues"],
                                                                                                 s[
                                                                                                     "max_n_gaps_in_TMD_subject_seq"]))
            homo_hetero_contact_file = os.path.join(s["thoipapy_data_folder"], "Features", "Structure", database, "{}.homohetero.bind.closedist.csv".format(acc))

            hetero_contact_num = remove_crystal_hetero_contact_residues(acc, feature_combined_file, homo_hetero_contact_file,no_hetero_feature_combined_file, logging)
            sum_hetero_contact = sum_hetero_contact + hetero_contact_num
    logging.info(
        "there are in crystal data set in total hetero contact residues: {}  ".format( sum_hetero_contact))

def remove_crystal_hetero_contact_residues(acc, feature_combined_file, homo_hetero_contact_file,no_hetero_feature_combined_file, logging):
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
        raise FileNotFoundError("homo_hetero_contact_file NOT FOUND. hetero contact residues not calculated and added to combined file.\n({})".format(homo_hetero_contact_file))

    hetero_contact_num = 0
    if os.path.isfile(homo_hetero_contact_file):
        df_contacts = pd.read_csv(homo_hetero_contact_file, index_col=0)
        combined_interface = df_combined.interface
        homo_hetero_interface = df_contacts.bind
        hetero_inter = [1 if homo_hetero_interface.iloc[i] == 1 and combined_interface.iloc[i] == 0 else 0 for i in
                        homo_hetero_interface.index]
        hetero_contact_num = hetero_inter.count(1)
        # df_addhetreo_combined["hetero_interface"] = hetero_inter
        # df_addhetreo_combined.to_csv(add_hetero_feature_combined_file)
        hetero_inter_index = []
        for i in range(len(hetero_inter)):
            if hetero_inter[i] == 1:
                hetero_inter_index.append(i)
        df_combined = df_combined.drop(df_combined.index[hetero_inter_index])
        #df_combined["interface"] = hetero_inter
        df_combined.to_csv(no_hetero_feature_combined_file)

        logging.info("{} add_hetero_contact_to_crystal_combined_files finished ({})".format(acc, no_hetero_feature_combined_file))

    else:
        logging.warning(
            "{} add_hetero_contact_to_crystal_combined_file failed, {} not found".format(acc, feature_combined_file))
    return hetero_contact_num


             ###################################################################################################
             #                                                                                                 #
             #            combining test data and add physical parameters                                      #
             #                                                                                                 #
             ###################################################################################################


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

    train_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_train_data.csv".format(s["setname"]))

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
        logging.warning("No experimental data has been added to this dataset. Hope you're not trying to train with it!!!")

    # reset the index to be a range (0,...).
    df_all.index = range(df_all.shape[0])

    # reorder the columns
    column_list = ['acc_db', 'interface', 'interface_score', 'residue_num', 'residue_name', 'n_homologues']
    df_all = thoipapy.utils.reorder_dataframe_columns(df_all, column_list)

    # # remove crystal hetero_interface residues and drop "hetero_interface" column
    # hetero_inter_index = []
    # for i in range(df_all.shape[0]):
    #     if df_all.loc[i,"hetero_interface"] == 1:
    #         hetero_inter_index.append(i)
    # df_all = df_all.drop(df_all.index[hetero_inter_index])
    # df_all.drop(["hetero_interface"],axis=1, inplace=True)
    df_all.to_csv(train_data_csv)
    logging.info('Finished creating train or test data for machine learning.')


