import csv
import os
import sys

import numpy as np
import pandas as pd

import thoipapy
from thoipapy import utils as utils
from thoipapy.utils import normalise_0_1


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

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_file = os.path.join(s["thoipapy_data_folder"], "features", "coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, logging)


def coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, logging):
    """Runs freecontact from command line on a multiple sequence alignment in FASTA format.

    Parameters
    ----------
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    freecontact_file : str
        Path to freecontact output file.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        try:
            thoipapy.utils.make_sure_path_exists(freecontact_file, isfile=True)
            exect_str = "grep -v '^>' {aln_file} |sed 's/[a-z]//g'|freecontact >{freecontact_output_file}".format(
                aln_file=path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_output_file=freecontact_file)

            command = utils.Command(exect_str)
            command.run(timeout=400, log_stderr=False)

            logging.info(f"coevolution_calculation_with_freecontact finished ({freecontact_file})")
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
    logging.info('starting parse_freecontact_coevolution_mult_prot')

    for i in df_set.index:
        sys.stdout.write(".")
        sys.stdout.flush()
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_start = int(df_set.loc[i, "TMD_start"])
        TMD_end = int(df_set.loc[i, "TMD_end"])
        freecontact_file = os.path.join(s["thoipapy_data_folder"], "features", "coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        freecontact_parsed_csv = os.path.join(s["thoipapy_data_folder"], "features", "coevolution", database, "{}.surr{}.gaps{}.freecontact_parsed.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, TMD_start, TMD_end, logging)
    sys.stdout.write("\n")
    logging.info('finished parse_freecontact_coevolution_mult_prot')


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
        NOTE! For the DI4cum etc, originally the dictionary for mi and di was switched.
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
    DI4cum = [0] * tmd_length
    DI8cum = [0] * tmd_length
    DItop4mean = [0] * tmd_length
    MItop4mean = [0] * tmd_length
    DItop8mean = [0] * tmd_length
    MItop8mean = [0] * tmd_length
    DImax_1 = [0] * tmd_length  # the sum of the top 8
    MImax_1 = [0] * tmd_length
    for key in dict_di_list:
        DItop4mean[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:4])) / 4
        MItop4mean[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:4])) / 4
        DItop8mean[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:8])) / 8
        MItop8mean[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:8])) / 8
        # MT note: as far as I can see, the DImax is not limited to [0:8], and is therefore the max DI value between
        # the residue of interest, and any of the residues in the TMD
        DImax_1[int(key) - 1] = sorted(dict_di_list[key], reverse=True)[0]
        MImax_1[int(key) - 1] = sorted(dict_mi_list[key], reverse=True)[0]
        # sys.stdout.write(str(key)+"corresponding to"+str(dict_di_list[key]))
    dict_di_value_sort = sorted(dict_di.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]
    dict_mi_value_sort = sorted(dict_mi.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]

    for i in range(0, 4):
        res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
        res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
        DI4cum[res0] = DI4cum[res0] + float(dict_di_value_sort[i][1])
        DI4cum[res1] = DI4cum[res1] + float(dict_di_value_sort[i][1])

    for i in range(0, 8):
        res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
        res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
        DI8cum[res0] = DI8cum[res0] + float(dict_di_value_sort[i][1])
        DI8cum[res1] = DI8cum[res1] + float(dict_di_value_sort[i][1])

    writer = csv.writer(freecontact_parsed_csv_handle, delimiter=',', quotechar='"',
                        lineterminator='\n',
                        quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
    writer.writerow(["residue_num", "residue_name", "DImax", "DItop4mean", "DItop8mean", "DI4cum", "DI8cum", "MImax", "MItop4mean", "MItop8mean"])
    for index in range(len(DI8cum)):
        csv_header_for_cumulative_strength_file = [(index + 1), dict_residuenum_residuename[(index + 1)],
                                                   DImax_1[index], DItop4mean[index], DItop8mean[index], DI4cum[index], DI8cum[index], MImax_1[index], MItop4mean[index],
                                                   MItop8mean[index]]
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

    # open up the existing parsed freecontact file, with DI4cum, DI8cum, etc
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
        # position_list = range(min_ - padding, max_ + padding + 1)
        position_list = range(TMD_start - padding, TMD_end + padding + 1)
        dfp = dfp.reindex(index=position_list, columns=position_list)
        # put data on both sides of the table for easy indexing
        for col in dfp.columns:
            start = col + 1
            dfp.loc[start:, col] = dfp.loc[col, start:]
        # drop rows with only nan
        # dfp.dropna(how="all", inplace=True)
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
        df_out["XIall_mean"] = dfp.mean(axis=1)

        for pos in range(TMD_start, TMD_end + 1):
            # iterate from i-1 and i+1 to i-5 and i+5
            for n in range(1, 6):
                i_minus = pos - n if pos - n in dfp.columns else dfp.columns.min()
                i_plus = pos + n if pos + n in dfp.columns else dfp.columns.max()
                # select the two datapoints (e.g. i-4 and i+4 relative to i)
                sel_XI_ser = dfp.loc[pos, [i_minus, i_plus]]
                df_out.loc[pos, "XI{}mean".format(n)] = sel_XI_ser.mean()

                if n == 4:
                    # add the direct connection between i-4 and i+4 (excluding i)
                    sel_XI_ser["m4_to_p4_value"] = dfp.loc[i_plus, i_minus]
                    df_out.loc[pos, "XI4value"] = sel_XI_ser.mean()

                    # calculate mean and max of all pairwise values between i and from i-4 to i+4 (total of 8 values)
                    m4_to_p4_ser = dfp.loc[pos, i_minus:i_plus]
                    df_out.loc[pos, "XI4_inclusive_mean"] = m4_to_p4_ser.mean()
                    df_out.loc[pos, "XI4_inclusive_max"] = m4_to_p4_ser.max()

        # HEPTAD MOTIF
        a, b, c, d, e, f, g = 1, np.nan, np.nan, 1, 1, np.nan, 1
        # extend the list longer than any TMD
        # e.g. [1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, 1, 1, nan, nan, 1, 1, nan, ......
        hep_list = [a, b, c, d, e, f, g] * 10

        highest_XI_face_value = 0

        # iterate through heptad faces, assuming perfect helix
        for face in range(7):
            # truncate heptad [1,nan...] list to match the length of the TMD
            end = dfp.shape[0] + face
            hep_list_trunc = hep_list[face: end]
            hep_list_trunc = np.array(hep_list_trunc)
            # get indices of the residues corresponding to that heptad motif
            # e.g. [ 88  91  92  94  95  98  99 101 102 105 106 108 109 112 113 115 116]
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
        df_out.loc[index_positions_highest_face, "XI_highest_face"] = 1
        df_out["XI_highest_face"] = df_out["XI_highest_face"].fillna(0).astype(int)

        # replace the XI in the column names with either MI or DI
        new_column_names = pd.Series(df_out.columns).str.replace("XI", XI)
        df_out.columns = new_column_names

    # normalise all columns except for residue_num and residue_name
    column_list = ["residue_num", "residue_name"]
    coev_colname_list = df_out.columns.tolist()[len(column_list):]

    # Specifically overwrite normalised values for Cum metrics. Convert to 0 or 1.
    coev_cum_colname_list = ["DI4cum", "DI8cum"]
    for col in coev_cum_colname_list:
        df_out[col] = df_out[col].apply(lambda x: 1 if x > 0 else 0)

    # rename so that original data is "NAME_raw" and normalised data is "NAME"
    for col in coev_colname_list:

        already_normalised = "cum" in col or "highest_face" in col

        if not already_normalised:
            df_out["{}_raw".format(col)] = df_out[col]
            df_out[col] = normalise_0_1(df_out[col])[0]

    df_out.to_csv(freecontact_parsed_csv)

    logging.info(f"parse_freecontact_coevolution finished ({freecontact_parsed_csv})")
