import csv
import os
import pandas as pd
import scipy as sc
from pandas import Series
import scipy.stats
import glob
import re
import math
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from difflib import SequenceMatcher
import thoipapy.utils as utils
import sys
import matplotlib.pyplot as plt
from weighslide import calculate_weighted_windows
import thoipapy
from math import isclose

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


def mem_a3m_homologues_filter(set_,logging):

    tmp_list_loc = set_["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0][0:5]
        tm_start = int(row.strip().split(",")[4]) + 1
        tm_end = int(row.strip().split(",")[4]) + (int(row.strip().split(",")[3]) - int(row.strip().split(",")[2]) + 1)
        #return tm_end
        homo_a3m_file = os.path.join(set_["output_oa3m_homologues"], "SinglePassTmd/%s.surr20.parse.a3m") % acc
        if os.path.isfile(homo_a3m_file):

            homo_a3m_file_handle=open(homo_a3m_file,"r")
            homo_filter_file=os.path.join(set_["output_oa3m_homologues"], "SinglePassTmd/%s.surr20.a3m.mem.uniq.2gaps") %acc
            homo_filter_file_handle = open(homo_filter_file,"w")
            homo_mem_lips_input_file = os.path.join(set_["output_oa3m_homologues"], "SinglePassTmd/%s.surr20.mem.lips.input") %acc
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            logging.info("starting parsing a3m file: %s\n" %homo_filter_file)
            i = 0
            tm_query=""
            for line in homo_a3m_file_handle:
                tm_str = line[(tm_start - 1):tm_end]
                if i == 0:
                    tm_query = tm_str
                    i = i + 1
                    print("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                    print("{}".format(tm_str), file=homo_filter_file_handle)
                    continue
                mean_hydrophobicity = calc_lipophilicity(tm_str)
                ratio = SequenceMatcher(None, tm_query, tm_str).ratio()
                if not re.search("-", tm_str) and not re.search("X",
                                                                tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and ratio < set_["max_identity_of_TMD_seq"]and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    print("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                gap_num = tm_str.count("-")
                if (gap_num <= 3 and not re.search("X",
                                                   tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and ratio < set_["max_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    print("{}".format(tm_str), file=homo_filter_file_handle)
                    # homo_filter_file_handle.write(line)
                    # homo_mem_lips_input_file_handle.write(tm_str
            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()
        homo_a3m_file_handle.close()
    tmp_file_handle.close()



def create_PSSM_from_MSA_mult_prot(set_, df_set, logging):
    logging.info('start pssm calculation')
    #
    # tmp_list_loc = set_["list_of_tmd_start_end"]
    # tmp_file_handle = open(tmp_list_loc, 'r')
    # # skip header
    # next(tmp_file_handle)
    # for row in tmp_file_handle:
    #     protein_name = row.strip().split(",")[0][0:6]
    #     database = row.strip().split(",")[6]
    #     # if protein_name=="4cbj_E":
    #     tm_start = int(row.strip().split(",")[2]) - 5  ###for fullseq
    #     if (tm_start <= 0):
    #         tm_start = 1
    #     tm_end = int(row.strip().split(",")[3]) + 5  ###for fullseq
    #     if tm_end > int(row.strip().split(",")[1]):
    #         tm_end = int(row.strip().split(",")[1])

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir,"{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, set_["num_of_sur_residues"],set_["max_n_gaps_in_TMD_subject_seq"]))
        pssm_csv = os.path.join(set_["feature_pssm"], database, "{}.surr{}.gaps{}.pssm.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, logging)

        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.txt".format(acc, set_["max_n_gaps_in_TMD_subject_seq"]))

        pssm_csv_surr5 = os.path.join(set_["feature_pssm"], database, "{}.surr5.gaps{}.pssm.csv".format(acc, set_["max_n_gaps_in_TMD_subject_seq"]))
        create_PSSM_from_MSA(path_uniq_TMD_seqs_surr5_for_LIPO, pssm_csv_surr5, acc, logging)

    #tmp_file_handle.close()

def create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, logging):

    #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologues"], database, "%s.a3m.mem.uniq.2gaps%s") % (acc, set_["surres"])

    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        #try:
        thoipapy.utils.make_sure_path_exists(pssm_csv, isfile=True)
        with open(pssm_csv, 'w') as pssm_file_handle:
            mat = []
            # with csv.writer(pssm_file_handle, delimiter=',', quotechar='"',
            #                     lineterminator='\n',
            #                     quoting=csv.QUOTE_NONNUMERIC, doublequote=True) as writer:
            writer = csv.writer(pssm_file_handle)
            writer.writerow(['residue_num', 'residue_name', 'A', 'I', 'L', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'M', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P'])
            # if os.path.isfile(homo_filter_fasta_file):
            # for line in open(homo_filter_fasta_file).readlines():
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, "r") as f:
                for line in f.readlines():
                    if not re.search("^>", line):
                        mat.append(list(line))
            rowlen = len(mat[0])
            collen = len(mat)
            column = []
            # write 20 amino acids as the header of pssm output file
            # pssm_file_handle.write(
            # 'residue'+' '+'A' + ' ' + 'I' + ' ' + 'L' + ' ' + 'V' + ' ' + 'F' + ' ' + 'W' + ' ' + 'Y' + ' ' + 'N' + ' ' + 'C' + ' ' + 'Q' + ' ' + 'M' + ' ' + 'S' + ' ' + 'T' + ' ' + 'D' + ' ' + 'E' + ' ' + 'R' + ' ' + 'H' + ' ' + 'K' + ' ' + 'G' + ' ' + 'P' + '\n')
            for j in range(0, rowlen - 1):
                for i in range(0, collen):
                    column.append(mat[i][j])
                aa_num = [column.count('A') / collen, column.count('I') / collen, column.count('L') / collen, column.count('V') / collen, column.count('F') / collen,
                          column.count('W') / collen, column.count('Y') / collen, column.count('N') / collen, column.count('C') / collen, column.count('Q') / collen,
                          column.count('M') / collen, column.count('S') / collen, column.count('T') / collen, column.count('D') / collen, column.count('E') / collen,
                          column.count('R') / collen, column.count('H') / collen, column.count('K') / collen, column.count('G') / collen, column.count('P') / collen]
                aa_num.insert(0, mat[0][j])  ###adding the residue name to the second column
                aa_num.insert(0, j + 1)  ##adding the residue number to the first column
                writer.writerow(aa_num)
                column = []
            logging.info('{} pssm calculation finished ({})'.format(acc, pssm_csv))

    else:
        logging.warning("{} homo_filter_fasta_file does not exist({})".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))


def lipo_from_pssm_mult_prot(set_, df_set, logging):
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

    thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    hydrophob_scale_path = os.path.join(thoipapy_module_path, "setting", "hydrophobicity_scales.xlsx")

    # set name of hydrophobicity scale
    # current options are KyteDoolittle, Wimley, Hessa, Elazar, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, Engelman(GES)
    scalename = "Hessa"
    failed_acc_list = []
    plot_linechart = True

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]

        tm_surr_left = df_set.loc[i, "tm_surr_left"]
        tm_surr_right = df_set.loc[i, "tm_surr_right"]
        # reset to 5. AM NOT SURE WHY.
        if tm_surr_left >= 5:
            tm_surr_left = 5
        if tm_surr_right >= 5:
            tm_surr_right = 5

        # pssm_csv = os.path.join(set_["feature_pssm"], database, "{}.mem.2gap.pssm_surr5.csv".format(acc))
        pssm_csv_surr5 = os.path.join(set_["feature_pssm"], database, "{}.surr5.gaps{}.pssm.csv".format(acc, set_["max_n_gaps_in_TMD_subject_seq"]))
        lipo_csv = os.path.join(set_["feature_lipophilicity"], database, "{}_{}_lipo.csv".format(acc, scalename))

        result_tuple = lipo_from_pssm(acc, pssm_csv_surr5, lipo_csv, tm_surr_left, tm_surr_right, scalename, hydrophob_scale_path, logging, plot_linechart)

        if result_tuple[1] is False:
            failed_acc_list.append(acc)

    if len(failed_acc_list) > 0:
        sys.stdout.write("\n\nlipo_from_pssm_mult_prot failed for the following : {}".format(failed_acc_list))
    logging.info('~~~~~~~~~~~~                 finished lipo_from_pssm_mult_prot              ~~~~~~~~~~~~')


def lipo_from_pssm(acc, pssm_csv_surr5, lipo_csv, tm_surr_left, tm_surr_right, scalename, hydrophob_scale_path, logging, plot_linechart=False):
    """Calculates lipophilicity from a PSSM for a single protein.

    Takes a PSSM as an input. The PSSM should have the fractions of each residue at each position.
    1) calculates the lipophilicity for each position (e.g. lipo_hessa)
    2) uses the weighslide module to get the average lipophilicity of the 3 N-terminal aa to each position (e.g. lipo_Hessa_mean_3_N_pos)
       and the 3 C-terminal aa to each position (e.g. lipo_Hessa_mean_3_C_pos).


    Parameters
    ----------
    set_ : dict
        settings dictionary
    acc : str
        UniProt or PDB accession.
    scalename : str
        Name of hydrophobicity scale. Should match the excel file with the scales.
        Current options are KyteDoolittle, Wimley, Hessa, Elazar, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, Engelman(GES)
    hydrophob_scale_path : str
        Path to excel file with hydrophobicity scales
    plot_linechart : bool
        Whether to plot the linecharts showing the hydrophobicity.
        Plotting the figure slows down the processing time.

    Saved files
    -----------
    lipo_csv - csv
        csv output file with the lipophilicity values for each position
        index = ["A1", "I2", ....]
        columns = lipo_Hessa, lipo_Hessa_mean_3_N_pos, lipo_Hessa_mean_3_C_pos (e.g. for Hessa scale)

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

    df_hs = pd.read_excel(hydrophob_scale_path, skiprows=2)
    df_hs.set_index("1aa", inplace=True)
    df_hs.sort_index(inplace=True)
    hs_arr = df_hs[scalename].as_matrix()

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
    df_lipo["lipo_{}".format(scalename)] = dfh.sum(axis=1)
    #df_lipo["IND"] = df_lipo.index
    df_lipo.index = range(df_lipo.shape[0])
    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [1, 1, 1, "x", 0, 0, 0]
    lipo_mean_3_N_pos = calculate_weighted_windows(df_lipo["lipo_{}".format(scalename)], window, statistic="mean", full_output=False)
    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [0, 0, 0, "x", 1, 1, 1]
    lipo_mean_3_C_pos = calculate_weighted_windows(df_lipo["lipo_{}".format(scalename)], window, statistic="mean", full_output=False)

    # replace positions with nan that could not be properly calculated
    # this inlcludes the first 3 positions of lipo_mean_3_N_pos, and last 3 positions of lipo_mean_3_C_pos
    # these nan positions should be OUTSIDE THE TMD
    lipo_mean_3_N_pos.iloc[0:3] = np.nan
    lipo_mean_3_C_pos.iloc[-3:] = np.nan

    """df_lipo now looks like this.

              lipo_hessa  IND  lipo_mean_3_N_pos  lipo_mean_3_C_pos
    position                                                                   
    0          -0.339400   V1                      NaN                 1.231000
    1           0.890600   S2                      NaN                 0.632300
    2           2.073000   P3                      NaN                -0.019160
    3           0.729400   G4                 0.437367                -0.120600
    ....

    28          2.495600  R29                 0.131967                 0.197203
    29         -0.437800  L30                 0.530000                      NaN
    30         -0.298367  V31                 0.688400                      NaN
    31          1.919388  P32                 0.586478                      NaN

    """

    df_lipo["lipo_{}_mean_3_N_pos".format(scalename)] = lipo_mean_3_N_pos
    df_lipo["lipo_{}_mean_3_C_pos".format(scalename)] = lipo_mean_3_C_pos


    df_lipo["residue_num"] = df["aa_pos"].values
    df_lipo["residue_name"] = df["residue_name"].values
    #df_lipo.index = range(df_lipo.shape[0])
    #df_lipo.index=[int(i) -tm_surr_left for i in df_lipo.index]
    df_lipo["residue_num"] = [int(i) -tm_surr_left for i in df_lipo["residue_num"]]
    df_lipo.index =df_lipo["residue_num"]

    df_lipo = df_lipo[["residue_name","lipo_{}".format(scalename), "lipo_{}_mean_3_N_pos".format(scalename),"lipo_{}_mean_3_C_pos".format(scalename)]]
    #df_lipo.set_index("IND", inplace=True)
    df_lipo=df_lipo[tm_surr_left:-tm_surr_right]

    # print(acc, "\n", df_lipo,"\n\n")
    # with pd.ExcelWriter(lipo_excel) as writer:
    #     df.to_excel(writer, sheet_name="df")
    #     dfa.to_excel(writer, sheet_name="dfa")
    #     df_lipo.to_excel(writer, sheet_name="df_lipo")

    df_lipo.to_csv(lipo_csv)

    if plot_linechart:
        # plot a linechart with lipo, lipo_mean_3_N_pos, lipo_mean_3_C_pos
        fig, ax = plt.subplots()
        df_lipo.plot(ax=ax)
        ax.set_ylabel("lipophilicity, {} scale".format(scalename))
        ax.set_xticks(range(len(df_lipo)))
        ax.set_xticklabels(df_lipo.index, rotation=90)
        ax.set_xlabel("")
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(lipo_linechart, dpi=200)
        plt.close("all")

    logging.info("{} lipo_from_pssm_mult_prot finished using {} scale ({})".format(acc, scalename, lipo_csv))
    return acc, True, ""


def entropy_calculation_mult_prot(set_, df_set, logging):
    logging.info('start entropy calculation')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologues"],database,"%s.a3m.mem.uniq.2gaps%s") % (acc,set_["surres"])
        alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        entropy_file = os.path.join(set_["feature_entropy"], database, "{}.surr{}.gaps{}.uniq.entropy.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, entropy_file, logging)


def entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, entropy_file, logging):
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):

        if not os.path.isdir(os.path.dirname(entropy_file)):
            os.makedirs(os.path.dirname(entropy_file))

        with open(entropy_file, 'w') as entropy_file_handle:
            mat = []
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
                for line in f.readlines():
                    # for line in open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT).readlines():
                    if not re.search("^>", line):
                        mat.append(list(line))
                rowlen = len(mat[0])
                collen = len(mat)
                column = []
            writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "Entropy"])
            for j in range(0, rowlen - 1):
                for i in range(0, collen):
                    column.append(mat[i][j])
                column_serie = Series(column)
                p_data = column_serie.value_counts() / len(column_serie)  # calculates the probabilities
                entropy = sc.stats.entropy(p_data)  # input probabilities to get the entropy
                csv_header_for_ncbi_homologues_file = [j + 1, mat[0][j], entropy]
                writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_ncbi_homologues_file)
                # entropy_file_handle.write(mat[0][j]+' '+ str(entropy)+'\n')
                column = []
            # entropy_file_handle.close()
            logging.info('{} entropy_calculation finished ({})'.format(acc, entropy_file))
    else:
        logging.warning("{} entropy_calculation failed. {} input file not found".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))


def coevolution_calculation_with_freecontact_mult_prot(set_, df_set, logging):
    logging.info('start coevolution calculation using freecontact')
    freecontact_loc=set_["freecontact_dir"]

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        print(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT))

        coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, freecontact_loc, logging)

def coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, freecontact_loc, logging):
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        try:
            # freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], "zpro/NoRedundPro/%s.mem.2gap.freecontact") % acc
            thoipapy.utils.make_sure_path_exists(freecontact_file, isfile=True)
            # freecontact_file_handle = open(freecontact_file, 'w')
            exect_str = "grep -v '^>' {aln_file} |sed 's/[a-z]//g'|{freecontact} >{freecontact_output_file}".format(
                aln_file=path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact=freecontact_loc, freecontact_output_file=freecontact_file)

            command = utils.Command(exect_str)
            # command=Command(exect_str)
            command.run(timeout=400)

            logging.info("Output file: %s\n" % freecontact_file)
        except:
            logging.warning("freecontact gives an error")
    else:
        logging.warning("{} does not exist".format(path_uniq_TMD_seqs_for_PSSM_FREECONTACT))




def parse_freecontact_coevolution_mult_prot(set_, df_set, logging):
    """

    Parameters
    ----------

    pathdict : dict

    set_ : dict
        ......
    logging : logging.Logger
        ...

    Returns
    -------

    """
    logging.info('cumulative co-evolutionary strength parsing')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_parsed_csv = os.path.join(set_["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, logging)

def parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, logging):
    if os.path.isfile(freecontact_file):
        try:
            dict_di = {}
            dict_mi = {}
            dict_di_list = {}
            dict_mi_list = {}
            s = "-"
            # freecontact_parsed_csv = "/scratch/zeng/Q6ZRP7.csv"
            freecontact_parsed_csv_handle = open(freecontact_parsed_csv, 'w')
            freecontact_file_handle = open(freecontact_file, 'r')
            dict_residuenum_residuename = {}
            for row in freecontact_file_handle:  # get the last line of the *.freecontact file which contains the tmd length infor
                residue_pairs = row.strip().split()
                dict_di[s.join([residue_pairs[0], residue_pairs[2]])] = residue_pairs[4]
                dict_mi[s.join([residue_pairs[0], residue_pairs[2]])] = residue_pairs[5]
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

            tmd_length = int(row.strip().split()[2])
            CumDI4 = [0] * tmd_length
            CumMI4 = [0] * tmd_length
            CumDI8 = [0] * tmd_length
            CumMI8 = [0] * tmd_length
            CoevDI4 = [0] * tmd_length
            CoevMI4 = [0] * tmd_length
            CoevDI8 = [0] * tmd_length
            CoevMI8 = [0] * tmd_length
            CoevDImax = [0] * tmd_length  # the sum of the top 8
            CoevMImax = [0] * tmd_length
            for key in dict_di_list:
                CoevDI4[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:4])) / 4
                CoevMI4[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:4])) / 4
                CoevDI8[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:8])) / 8
                CoevMI8[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:8])) / 8
                CoevDImax[int(key) - 1] = sorted(dict_di_list[key], reverse=True)[0]
                CoevMImax[int(key) - 1] = sorted(dict_mi_list[key], reverse=True)[0]
                # print(str(key)+"corresponding to"+str(dict_di_list[key]))
            # print(residue_di_new)
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
            writer.writerow(["residue_num", "residue_name", "CoevDImax", "CoevDI4", "CoevDI8", "CumDI4", "CumDI8", "CoevMImax", "CoevMI4", "CoevMI8", "CumMI4", "CumMI8"])
            for index in range(len(CumDI8)):
                csv_header_for_cumulative_strength_file = [(index + 1), dict_residuenum_residuename[(index + 1)],
                                                           CoevDImax[index], CoevDI4[index], CoevDI8[index], CumDI4[index], CumDI8[index], CoevMImax[index], CoevMI4[index],
                                                           CoevMI8[index], CumMI4[index], CumMI8[index]]
                # writer = csv.writer(freecontact_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
                # quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_cumulative_strength_file)
            freecontact_file_handle.close()
            freecontact_parsed_csv_handle.close()
            # freecontact_parsed_csv_handle.write(str(residue_di[index])+"\t"+str(residue_mi[index])+"\n")
        except:
            print("coevolution parsing gave an error")
    else:
        logging.warning("{} parse_freecontact_coevolution failed, {} not found.".format(acc, freecontact_file))


def relative_position_calculation(set_, df_set, logging):
    """calculate the residue relative position on the TMD

    Parameters
    ----------
    set_
    logging

    Returns
    -------
    Rp1 : float
    position relative to the tmd
    Rp2 : float
    position relative to the full sequence

    """
    logging.info('start to calculate the relative positions')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #tm_start = int(row.strip().split(",")[2])
        #seq_len = int(row.strip().split(",")[1])
        tm_start = df_set.loc[i, "TMD_start"]
        tm_seq = df_set.loc[i, "full_seq"]
        seq_len = df_set.loc[i, "seqlen"]

        relative_position_file = os.path.join(set_["feature_relative_position"], database, "%s.relative_position%s.csv") % (acc, set_["surres"])
        thoipapy.utils.make_sure_path_exists(relative_position_file, isfile=True)
        alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
            relative_position_file_handle = open(relative_position_file, 'w')
            mat = []
            writer = csv.writer(relative_position_file_handle, delimiter=',', quotechar='"',
                                lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(
                ['residue_num', 'residue_name', 'RelPos_TMD', 'RelPos_fullseq'])
            # if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
            # for line in open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT).readlines():

            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
                mat = []
                for line in f.readlines():
                    if not re.search("^>", line):
                        mat.append(list(line))
                tm_seq = mat[0]
                tm_len=len(tm_seq)
                for i in range(1,tm_len):
                    rp1=i/tm_len
                    rp2=(i+tm_start-1)/seq_len
                    writer.writerow([i,tm_seq[i-1],rp1,rp2])
            relative_position_file_handle.close()
            logging.info('{} relative position calculation finished ({})'.format(acc, relative_position_file))


def LIPS_score_calculation_mult_prot(set_, df_set, logging):
    logging.info('start lips score calculation')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #LIPS_input_file = os.path.join(set_["output_oa3m_homologues"],database, "%s.mem.lips.input%s") % (acc,set_["surres"])
        alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, set_["num_of_sur_residues"]))

        if os.path.isfile(path_uniq_TMD_seqs_no_gaps_for_LIPS):
            # LIPS_output_file = os.path.join(set_["feature_lips_score"], "zpro/NoRedundPro/%s.mem.lips.output") % acc
            #path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, set_["num_of_sur_residues"]))

            #LIPS_output_file = os.path.join(set_["feature_lips_score"], database, "%s.mem.lips.output%s") % (acc, set_["surres"])

            LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, set_["num_of_sur_residues"]))

            LIPS_score_calculation(path_uniq_TMD_seqs_no_gaps_for_LIPS, LIPS_output_file)
        else:
            logging.warning("{} path_uniq_TMD_seqs_no_gaps_for_LIPS not found".format(acc))


def LIPS_score_calculation(input_seq_file, LIPS_output_file):

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
            print("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
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
                print("%3s" % rn, res, "%6.3f" % prop,
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
                        print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
                              file=LIPS_output_file_handle)
                        # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
                j = j + 7
        for i in range(4, 7):  # for surfaces from 4 to 6
            print("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
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
                print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j], file=LIPS_output_file_handle)
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
                        print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
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
                    print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k], file=LIPS_output_file_handle)
                    # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                k = k + 1
        print("SURFACE LIPOPHILICITY ENTROPY   LIPS", file=LIPS_output_file_handle)
        # LIPS_output_file_handle.write("SURFACE LIPOPHILICITY ENTROPY   LIPS")
        for i in sumim.keys():
            avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
            avpim = avpim * 2
            ave = sume[i] / aanum[i]  # average entropy for surface i
            peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
            print("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
                  "%8.3f" % peim,
                  file=LIPS_output_file_handle)  # print seven surfaces and see which surface with lowewst LIPS score
            # LIPS_output_file_handle.write("%s" % i, "%10.3f" % avpim, "%8.3f" % ave, "%8.3f" % peim)

def parse_LIPS_score_mult_prot(set_, df_set, logging):
    logging.info('start parsing lips output to cons and lips scores')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #LIPS_output_file = os.path.join(set_["feature_lips_score"],database,"%s.mem.lips.output%s") % (acc,set_["surres"])
        alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, set_["num_of_sur_residues"]))
        LIPS_parsed_csv = os.path.join(set_["feature_lips_score"], database, "{}.surr{}.LIPS_score_parsed.csv".format(acc, set_["num_of_sur_residues"]))
        parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging)

def parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging):

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

                # LIPS_output_handle.close()

                writer = csv.writer(LIPS_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num", "residue_name", "LIPS_lipo", "LIPS_entropy", "LIPS_surface"])
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


def convert_bind_data_to_csv(set_, df_set, logging):
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        bind_file = os.path.join(set_["structure_bind"], "%s.4.0closedist") % acc
        csv_output_file = os.path.join(set_["structure_bind"], "%s.4.0closedist.csv") % acc
        if os.path.isfile(bind_file):
            try:
                with open(bind_file,"r") as bind_file_handle:
                    #csv_output_file=os.path.join(set_["structure_bind"],"NoRedundPro/%s.csv") %acc
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
                print("bind file parsing occures errors")
        else:
            sys.stdout.write("{} convert_bind_data_to_csv failed. {} not found".format(acc, bind_file))


                         ###################################################################################################
                         #                                                                                                 #
                         #            combining train data and add physical parameters                                     #
                         #                                                                                                 #
                         ###################################################################################################

def combine_csv_files_with_features(set_, df_set, logging):
    logging.info('Combining features into traindata')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        scalename = set_["lipophilicity_scale"]
        lipo_csv = os.path.join(set_["feature_lipophilicity"], database, "{}_{}_lipo.csv".format(acc, scalename))
        relative_position_file = os.path.join(set_["feature_relative_position"], database, "%s.relative_position%s.csv") % (acc, set_["surres"])
        LIPS_parsed_csv = os.path.join(set_["feature_lips_score"], database, "{}.surr{}.LIPS_score_parsed.csv".format(acc, set_["num_of_sur_residues"]))
        pssm_csv = os.path.join(set_["feature_pssm"], database, "{}.surr{}.gaps{}.pssm.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        entropy_file = os.path.join(set_["feature_entropy"], database, "{}.surr{}.gaps{}.uniq.entropy.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        freecontact_parsed_csv = os.path.join(set_["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        feature_combined_file = os.path.join(set_["RF_features"], "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        thoipapy.utils.make_sure_path_exists(feature_combined_file, isfile=True)


        for n, filepath in enumerate([entropy_file, pssm_csv, lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv]):
            if not os.path.isfile(filepath):
                raise FileNotFoundError("combine_csv_files_with_features failed. {} not found".format(filepath))

        entropy_file_df = pd.read_csv(entropy_file)


        pssm_csv_df = pd.read_csv(pssm_csv)
        lipophilicity_file_df=pd.read_csv(lipo_csv)
        freecontact_parsed_csv_df = pd.read_csv(freecontact_parsed_csv)
        relative_position_file_df = pd.read_csv(relative_position_file)
        LIPS_parsed_csv_df=pd.read_csv(LIPS_parsed_csv)
        merge1=entropy_file_df.merge(pssm_csv_df,on=['residue_num','residue_name'])
        merge2=merge1.merge(lipophilicity_file_df,on=['residue_num','residue_name'])
        merge3=merge2.merge(freecontact_parsed_csv_df,on=["residue_num","residue_name"])
        merge4=merge3.merge(relative_position_file_df, on=["residue_num","residue_name"])
        merge5=merge4.merge(LIPS_parsed_csv_df,on=["residue_num","residue_name"])
        merge5.to_csv(feature_combined_file)


        file_list = ["entropy_file", "pssm_csv", "lipo_csv", "freecontact_parsed_csv", "relative_position_file", "LIPS_parsed_csv"]
        df_list = ["entropy_file_df", "pssm_csv_df", "lipophilicity_file_df", "freecontact_parsed_csv_df",  "relative_position_file_df", "LIPS_parsed_csv_df"]

        test_indexing = False
        if test_indexing:
            sys.stdout.write("{}".format(entropy_file_df))
            #entropy_file_df.loc[0, "residue_name"] = "X"
            #entropy_file_df.loc[0, "Entropy"] = 9999
            #entropy_file_df.loc[10, "Entropy"] = 7777
            #entropy_file_df["residue_num"] = range(3, entropy_file_df.shape[0] + 3)
            for df in [merge1, merge2, merge3, merge4, merge5]:
                print(df.shape)
            for n, df in enumerate([entropy_file_df, pssm_csv_df, lipophilicity_file_df,freecontact_parsed_csv_df,  relative_position_file_df, LIPS_parsed_csv_df]):
                dfname = df_list[n]
                TMD_seq_this_df = df.residue_name.str.cat()
                sys.stdout.write("\n{} ({})".format(TMD_seq_this_df, dfname))

        # Raise an error if the TMD sequence does not match original seq in settings file
        TMD_seq_in_merged_file = merge5.residue_name.str.cat()
        if TMD_seq != TMD_seq_in_merged_file:
            sys.stdout.write("acc = {}\noriginal = {}".format(acc, TMD_seq))
            sys.stdout.write("\nmerged   = {}".format(TMD_seq_in_merged_file))
            raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

        logging.info("combine_csv_files_with_features finished ({})".format(feature_combined_file))


def add_bind_data_to_combined_features(set_, df_set, logging):
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        feature_combined_file_incl_phys_param = os.path.join(set_["RF_features"], "combined", database,
        "{}.surr{}.gaps{}.combined_features_incl_phys_param.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))

        df_combined = pd.read_csv(feature_combined_file_incl_phys_param,index_col=0)
        bind_status_file = os.path.join(set_["RF_features"],'Structure/%s/%s.bind.closedist.csv') %(database,acc)
        if os.path.isfile(bind_status_file):
            df_bind = pd.read_csv(bind_status_file)
            df_combined_plus_bind = df_bind.merge(df_combined, on=["residue_num", "residue_name"])

            TMD_seq_in_merged_file = df_combined_plus_bind.residue_name.str.cat()
            if TMD_seq != TMD_seq_in_merged_file:
                sys.stdout.write("acc = {}\noriginal = {}".format(acc, TMD_seq))
                sys.stdout.write("\nmerged   = {}".format(TMD_seq_in_merged_file))
                raise IndexError("TMD_seq in original settings file and final merged features dataframe does not match.")

            # overwrite existing combined features file
            df_combined_plus_bind.to_csv(feature_combined_file_incl_phys_param)
            logging.info("{} add_bind_data_to_combined_features finished ({})".format(acc, feature_combined_file_incl_phys_param))

        else:
            logging.warning("{} add_bind_data_to_combined_features failed, {} not found".format(acc, bind_status_file))

def adding_physical_parameters_to_train_data(set_, df_set, logging):
    logging.info('adding physical parameters into traindata')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        dict = {}
        thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
        physical_parameter_file = os.path.join(thoipapy_module_path, "setting", "Physical_Property_csv.txt")

        with open(physical_parameter_file,"r") as physical_parameter_file_handle:
            feature_combined_file = os.path.join(set_["RF_features"], "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
            feature_combined_file_incl_phys_param = os.path.join(set_["RF_features"], "combined", database, "{}.surr{}.gaps{}.combined_features_incl_phys_param.csv".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
            if os.path.isfile(feature_combined_file):
                #try:
                with open(feature_combined_file,"r") as train_data_file_handle:
                    #train_data_add_physical_parameter_file = os.path.join("/scratch/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.physipara.traindata1.csv") %acc
                    train_data_add_physical_parameter_file_handle=open(feature_combined_file_incl_phys_param,"w")
                    #train_data_physical_parameter_file_handle = open(r"/scratch/zeng/thoipapy/Features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
                    writer = csv.writer(train_data_add_physical_parameter_file_handle, delimiter=',', quotechar='"',
                                        lineterminator='\n',
                                        quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                    for row in physical_parameter_file_handle:
                        if re.search("^Hydro", row):
                            continue
                        else:
                            array = row.split()
                            # print(array[1])
                            dict[array[0]] = array[1:16]
                            # for k, v in dict.items():
                            # print(k,v)
                    for row1 in train_data_file_handle:
                        # print(row1)
                        if re.search("residue_num", row1):
                            array1 = row1.rstrip().split(",")
                            array1[42:14]=["Hydrophobicity_sAA", "Charge_sAA", "PI_sAA", "LIPSI_sAA", "LIPSM_sAA", "Hydrophobic_sAA", "Aliphatic_sAA", "Aromatic_sAA", "Polar_sAA","Negative_sAA", "Positive_sAA", "Small_sAA", "Cbbranched_sAA", "Mass_sAA", "Volume_sAA"]
                            #array2 = array1[0:31]
                            #array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                            writer.writerow(array1)
                        else:
                            array3 = row1.rstrip().split(",")
                            array3[42:14] = dict[array3[2]]
                            writer.writerow(array3)
                    train_data_add_physical_parameter_file_handle.close()
                    logging.info("adding_physical_parameters_to_train_data finished. ({})".format(feature_combined_file_incl_phys_param))

                df = pd.read_csv(feature_combined_file_incl_phys_param, index_col=0)
                cols_to_delete = ["Charge_sAA","LIPSI_sAA","LIPSM_sAA","Hydrophobic_sAA", "Aliphatic_sAA", "Negative_sAA", "Positive_sAA", "Volume_sAA"]
                df.drop(cols_to_delete, axis=1, inplace=True)
                df.to_csv(feature_combined_file_incl_phys_param)

            else:
                logging.warning("{} does not exist".format(feature_combined_file))

                                             ###################################################################################################
                                             #                                                                                                 #
                                             #            combining test data and add physical parameters                                      #
                                             #                                                                                                 #
                                             ###################################################################################################

# DEPRECATED. USE SAME FUNCTION FOR BOTH TRAIN AND TEST DATASETS
# def adding_physical_parameters_to_test_data(set_, df_set, logging):
#     logging.info('adding physical parameters into testdata')
#
#     for i in df_set.index:
#         acc = df_set.loc[i, "acc"]
#         database = df_set.loc[i, "database"]
#         dict = {}
#         physical_parameter_file=os.path.join(set_["feature_physical_parameters"],"PhysicalProperty.txt")
#         physical_parameter_file_handle=open(physical_parameter_file,"r")
#         train_data_file=os.path.join(set_["RF_features"],database,"%s.mem.2gap.testdata%s.csv") %(acc,set_["surres"])
#         if os.path.isfile(train_data_file):
#             train_data_file_handle=open(train_data_file,"r")
#             train_data_add_physical_parameter_file=os.path.join(set_["RF_features"],database,"%s.mem.2gap.physipara.testdata%s.csv") %(acc,set_["surres"])
#             train_data_add_physical_parameter_file_handle=open(train_data_add_physical_parameter_file,"w")
#             writer = csv.writer(train_data_add_physical_parameter_file_handle, delimiter=',', quotechar='"',
#                                 lineterminator='\n',
#                                 quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
#             for row in physical_parameter_file_handle:
#                 if re.search("^Hydro", row):
#                     continue
#                 else:
#                     array = row.split()
#                     # print(array[1])
#                     dict[array[0]] = array[1:16]
#                     # for k, v in dict.items():
#                     # print(k,v)
#             for row1 in train_data_file_handle:
#                 # print(row1)
#                 if re.search("residue_num", row1):
#                     array1 = row1.rstrip().split(",")
#                     array1[42:14]=["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volume"]
#                     #array2 = array1[0:31]
#                     #array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
#                     writer.writerow(array1)
#                 else:
#                     array3 = row1.rstrip().split(",")
#                     array3[42:14] = dict[array3[2]]
#                     writer.writerow(array3)
#             train_data_add_physical_parameter_file_handle.close()
#             train_data_file_handle.close()
#         physical_parameter_file_handle.close()



def combine_all_train_data_for_random_forest(set_,logging):
    logging.info('creating train and test data for random forest')

    crystal_train_data_files = glob.glob(os.path.join(set_["RF_features"], "combined/crystal/*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))
    nmr_train_data_files = glob.glob(os.path.join(set_["RF_features"], "combined/NMR/*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))
    train_data_add_physical_parameter_files = glob.glob(os.path.join(set_["RF_features"], "combined/NMR/*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))) + glob.glob(os.path.join(set_["RF_features"], "Crystal/*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))
    etra_test_data_files = glob.glob(os.path.join(set_["RF_features"], "combined/ETRA/*.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))
    delete_crystal_traindata_files = []

    delete_acc_lists = ["4gyc_B","4hod_A","2hyn_A","4qe9_A","4roc_A"]
    for acc in delete_acc_lists:
        delete_crystal_traindata_files.append(os.path.join(set_["RF_features"],"combined/crystal","%s.surr{}.gaps{}.combined_features.csv".format( set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"])))

    header_saved = False

    crystal_train_data_after_delete_files = [item for item in crystal_train_data_files if item not in delete_crystal_traindata_files]
    train_data_after_delete_files = [item for item in train_data_add_physical_parameter_files if item not in delete_crystal_traindata_files]

    RF_dir = set_["RF_loc"]
    Crystal_Traindata_csv = os.path.join(RF_dir, 'Crystal_Traindata.csv')
    Nmr_Traindata_csv = os.path.join(RF_dir, 'Nmr_Traindata.csv')
    Crystal_Nmr_Traindata_csv = os.path.join(RF_dir, 'Crystal_Nmr_Traindata.csv')
    Etra_Testdata_csv = os.path.join(RF_dir, 'Etra_Testdata.csv')

    with open(Crystal_Traindata_csv , 'w')   as fout :

        for filename in crystal_train_data_after_delete_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)

    header_saved = False
    #with open('/scratch/zeng/thoipapy/RandomForest/Nmr_Traindata.csv' , 'w')   as fout :
    with open(Nmr_Traindata_csv, 'w')   as fout :

        for filename in nmr_train_data_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)

    header_saved = False
    with open(Crystal_Nmr_Traindata_csv, 'w')   as fout :

        for filename in train_data_after_delete_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)

    header_saved = False
    with open(Etra_Testdata_csv, 'w')   as fout :

        for filename in etra_test_data_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)