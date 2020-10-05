import math
import os
import sys
from math import isclose

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from weighslide import calculate_weighted_windows


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

        # pssm_csv = os.path.join(s["thoipapy_data_folder"], "features", "pssm", database, "{}.mem.2gap.pssm_surr5.csv".format(acc))
        pssm_csv_surr5 = os.path.join(s["thoipapy_data_folder"], "features", "pssm", database, "{}.surr5.gaps{}.pssm.csv".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        lipo_csv = os.path.join(s["thoipapy_data_folder"], "features", "lipophilicity", database, "{}_{}_lipo.csv".format(acc, scalename))

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

    # FLEXIBLE HYDROBICITY SCALE CURRENTLY NOT WORKING DUE TO PROBLEMS WITH LOCATING SETTINGS FILE IN DOCKER
    # thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    # hydrophob_scale_path = os.path.join(thoipapy_module_path, "setting", "hydrophobicity_scales.xlsx")
    # df_hs = pd.read_excel(hydrophob_scale_path, skiprows=2)
    # df_hs.set_index("1aa", inplace=True)
    # df_hs.sort_index(inplace=True)
    # hs_arr = df_hs[scalename].to_numpy()

    # hard-coded Engelman (GES) hydrophobicity scale
    # if re-implementing flexible scale, use the csv instead, or hard-code the various array values into python
    # "FAULTS: error: can't copy 'setting\hydrophobicity_scales.xlsx': doesn't exist or not a regular file"
    assert scalename == "Engelman(GES)"
    hs_arr = np.array([1.6,2.,-9.2,-8.2,3.7,1.,-3.,3.1,-8.8,2.8,3.4,-4.8,-0.2,-4.1,-12.3,0.6,1.2,2.6,1.9,-0.7])

    """The series should look like this for the Hessa scale. 
    For speed, this is typically converted to a numpy array, sorted alphabetically according to the residue.

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

    df_lipo: pd.DataFrame = pd.DataFrame()
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
        df_lipo: pd.DataFrame = - df_lipo
    else:
        raise ValueError("Panic. The scalename is wrong somehow.")

    # add the lowest value, so polarity starts at zero and gets larger
    # e.g. lowest value in Hessa scale is -0.6 for Ile
    min_value_dict = {"Hessa" : 0.6, "Elazar" : 1.92, "Hopp - Woods" : 3.4, "KyteDoolittle" : 4.5, "Wimley" : 1.85,
                      "Cornette" : 5.7, "Eisenberg" : 1.38, "Rose" : 0.91, "Janin" : 0.9, "Engelman(GES)" : 3.7}
    lowest_polarity_value = min_value_dict[scalename]

    # normalise the data by adding the lowes polarity value, so that all values are at least zero
    df_lipo: pd.DataFrame = df_lipo + lowest_polarity_value

    # normalise the data by applying a square root (usually to power of 1/4). Values must be positive.
    assert df_lipo.min().min() >= 0
    df_lipo: pd.DataFrame = df_lipo.applymap(lambda x: math.pow(x, 1/4))

# lowest_polarity_value = 0.6
    # columns = ["polarity", "polarity3Nmean", "polarity3Cmean", "polarity1mean"]
    # for col in columns:
    #     df_lipo[col] = df_lipo[col] + lowest_polarity_value

    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [1, 1, 1, "x", 0, 0, 0]
    polarity_i1_i3_N: pd.Series = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    # take mean over a window that includes the 3 N-terminal residues to the original position
    window = [0, 0, 0, "x", 1, 1, 1]
    polarity_i1_i3_C: pd.Series = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    window = [1, 1, 1]
    polarity1mean: pd.Series = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    # calculate polarity of central position relative to 6 surrounding residues
    window = [1, 1, 1, "x", 1, 1, 1]
    mean_polarity_surr_6_res: pd.Series = calculate_weighted_windows(df_lipo["polarity"], window, statistic="mean", full_output=False)
    relative_polarity: pd.Series = df_lipo["polarity"] / mean_polarity_surr_6_res

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
    #    lipo_excel = lipo_csv[:-4] + ".xlsx"
    #    lipo_linechart = lipo_csv[:-4] + "_linechart.png"
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