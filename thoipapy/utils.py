#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities file containing useful functions.
More recent functions are at the top.
Authors: Mark Teese, Rimma Jenske
Created on Fri Nov  8 15:45:06 2013
"""
import csv
import errno
import glob
import logging
import os
import re as re
import subprocess
import sys
import tarfile
import threading
import numpy as np
import pandas as pd
from shutil import copyfile
from time import strftime
import platform
import matplotlib.colors as colors
import ctypes

from scipy.special import comb


class Command(object):
    '''
    subprocess for running shell commands in win and linux
    This will run commands from python as if it was a normal windows console or linux terminal.
    taken from http://stackoverflow.com/questions/17257694/running-jar-files-from-python)'
    '''
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout, log_stderr=True):
        def target():
            #logging.info('Thread started')
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #self.process.communicate()
            stdout, stderr = self.process.communicate() # from http://stackoverflow.com/questions/14366352/how-to-capture-information-from-executable-jar-in-python
            # Thus far, SIMAP has only ever given java faults, never java output. Don't bother showing.
            # if the console prints anything longer than 5 characters, log it
            if len(stderr.decode("utf-8")) > 5:
                if log_stderr:
                    logging.warning('FAULTS: %s' % stderr.decode("utf-8"))
            #logging.info('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            logging.info('Terminating process')
            self.process.terminate()
            thread.join()
        # simply returns 0 every time it works. Waste of logging space! :)
        #logging.info(self.process.returncode)

def run_command(command):
    #this stopped working for some reason. Did I mess up a path variable?
    p = subprocess.Popen(command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')



def aaa(df_or_series):
    """ Function for use in debugging.
    Saves pandas Series or Dataframes to a user-defined csv file.
    """
     # convert any series to dataframe
    if isinstance(df_or_series, pd.Series):
        df_or_series = df_or_series.to_frame()
    csv_out = r"D:\data\000_aaa_temp_df_out.csv"
    df_or_series.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)

def make_sure_path_exists(path, isfile=False):
    """ If path to directory or folder doesn't exist, creates the necessary folders.

    Parameters
    ----------
    path : str
        Path to desired directory or file.
    isfile :
        If True, the path is to a file, and the subfolder will be created if necessary
    """
    if isfile:
        path = os.path.dirname(path)
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def reorder_dataframe_columns(dataframe, cols, front=True):
    '''Takes a dataframe and a subsequence of its columns,
       returns dataframe with seq as first columns if "front" is True,
       and seq as last columns if "front" is False.
       taken from https://stackoverflow.com/questions/12329853/how-to-rearrange-pandas-column-sequence

    Parameters
    ----------
    df : pd.DataFrame
        original dataframe
    cols : list
        List of cols to put first or last.
    front : bool
        Whether to put the columns at the front or the back.

    Usage
    -----
    df = reorder_dataframe_columns(df, ["TMD_start", "database", "whatever other col I want first"])
    '''
    for col in cols:
        if col not in dataframe.columns:
            cols.remove(col)
    cols_to_place_first = cols[:] # copy so we don't mutate seq
    for existing_col in dataframe.columns:
        if existing_col not in cols_to_place_first:
            if front: #we want "seq" to be in the front
                #so append current column to the end of the list
                cols_to_place_first.append(existing_col)
            else:
                #we want "seq" to be last, so insert this
                #column in the front of the new column list
                #"cols" we are building:
                cols_to_place_first.insert(0, existing_col)
    return dataframe[cols_to_place_first]

def create_tarballs_from_xml_files_in_folder(xml_dir, download_date="2017.11.02"):
    """script to create .tar.gz compressed files from a folder of xml files

    COPY SECTIONS INTO JUPYTER NOTEBOOK AND MODIFY AS NEEDED

    to add today's date
    date = strftime("%Y.%m.%d")

    THIS SCRIPT DOES NOT DELETE THE ORIGINAL XML FILES

    """

    xml_list = glob.glob(os.path.join(xml_dir, "*.xml"))

    for xml in xml_list:
        xml_renamed = xml[:-4] + ".BLAST.xml"
        xml_tar_gz = xml[:-4] + ".BLAST.xml.tar.gz"
        xml_txt = xml[:-4] + "_details.txt"
        #xml_txt = xml[:-4] + ".BLAST_details.txt"

        if not os.path.isfile(xml_tar_gz):
            copyfile(xml, xml_renamed)
            acc = os.path.basename(xml).split(".")[0]
            sys.stdout.write("{}, ".format(acc))
            sys.stdout.flush()
            # create an empty text file with the download date
            date = strftime("%Y%m%d")
            with open(xml_txt, "w") as f:
                f.write("acc\t{}\ndownload_date\t{}\ndatabase\tncbi_nr\ne_value\t1\n".format(acc, download_date))

            with tarfile.open(xml_tar_gz, mode='w:gz') as tar:
                # add the files to the compressed tarfile
                tar.add(xml_renamed, arcname=os.path.basename(xml_renamed))
                tar.add(xml_txt, arcname=os.path.basename(xml_txt))

            # delete the original files
            try:
                os.remove(xml_renamed)
                os.remove(xml_txt)
            except:
                sys.stdout.write("{} could not be deleted".format(xml_renamed))

def delete_BLAST_xml(blast_xml_file):
    """Small function to remove files that are already compressed in a tarball.

    Also deletes the associated text file with the date.

    Parameters
    ----------
    blast_xml_file : str
        Path to BLAST xml file
    """
    xml_txt = blast_xml_file[:-4] + "_details.txt"

    # delete the original files
    try:
        os.remove(blast_xml_file)
    except:
        sys.stdout.write("{} could not be deleted".format(blast_xml_file))
    try:
        os.remove(xml_txt)
    except:
        sys.stdout.write("{} could not be deleted".format(xml_txt))

def setup_biopol_plotly(username, api_key):
    import plotly
    plotly.tools.set_config_file(world_readable=False, sharing='private')
    plotly.tools.set_credentials_file(username=username, api_key=api_key)

def get_n_of_gaps_at_start_and_end_of_seq(seq):
    start = 0
    end = 0
    for aa in seq:
        if aa == "-":
            start += 1
        else:
            break
    for aa in seq[::-1]:
        if aa == "-":
            end += 1
        else:
            break
    return start, end


def get_list_residues_in_motif(seq, motif_ss, motif_len):
    # list of positions that matched the start of a motif
    match_start_list = []
    match_end_list = [0] * motif_len
    # counter for the matches
    match_number = 0
    result_dict = {}

    for start in range(len(seq) - motif_len):
        # for a SmallxxxSmall motif, the end is 4 residues later
        end = start + motif_len
        # get the matched segment
        segment = seq[start:end + 1]
        # check if the segment contains a motif
        match = re.match(motif_ss, segment)
        if match:
            # classify position as start of a motif
            match_start_list.append(1)
            match_end_list.append(1)
        else:
            match_start_list.append(0)
            match_end_list.append(0)
    # add the final zeros on the end of the list, so it's length matches the original sequence
    match_start_list = match_start_list + [0] * motif_len

    match_start_arr = np.array(match_start_list)
    match_end_arr = np.array(match_end_list)

    list_residues_in_motif = match_start_arr + match_end_arr
    list_residues_in_motif[list_residues_in_motif > 1] = 1

    # sys.stdout.write(seq, "seq")
    # sys.stdout.write("".join([str(x) for x in match_start_list]), "match_start_list")
    # sys.stdout.write("".join([str(x) for x in match_end_list]), "match_end_list")
    # sys.stdout.write("".join([str(x) for x in list_residues_in_motif]), "list_residues_in_motif")

    return list_residues_in_motif

def slice_TMD_seq_pl_surr(df_set):
    # note that due to uniprot-like indexing, the start index = start-1
    return df_set['full_seq'][int(df_set['TMD_start_pl_surr'] - 1):int(df_set['TMD_end_pl_surr'])]

def create_column_with_TMD_plus_surround_seq(df_set, num_of_sur_residues):
    df_set["TMD_start_pl_surr"] = df_set.TMD_start - num_of_sur_residues
    df_set.loc[df_set["TMD_start_pl_surr"] < 1, "TMD_start_pl_surr"] = 1
    df_set["TMD_end_pl_surr"] = df_set.TMD_end + num_of_sur_residues
    for i in df_set.index:
        #acc = df_set.loc[i, "acc"]
        if df_set.loc[i, "TMD_end_pl_surr"] > df_set.loc[i, "seqlen"]:
            df_set.loc[i, "TMD_end_pl_surr"] = df_set.loc[i, "seqlen"]
    TMD_seq_pl_surr_series = df_set.apply(slice_TMD_seq_pl_surr, axis=1)
    return df_set, TMD_seq_pl_surr_series


def create_namedict(names_excel_path, style="shortname [acc-db]"):
    """ Create protein name dictionary from an excel file with detailed protein info.

    e.g. namedict[P02724-NMR

    Parameters
    ----------
    names_excel_path : str
        Path to excel file with manually edited protein names
    style : str
        Style of protein name in output dictionary. Current options are "shortname [acc-db]" or "shortname [acc]"
        "shortname [acc-db]" = 'Q9Y286-ETRA': 'Siglec7 [Q9Y286-ETRA]'
        "shortname [acc]" = 'Q9Y286-ETRA': 'Siglec7 [Q9Y286]'

    Returns
    -------
    namedict : dict
        Dictionary in format namedict[acc_db] = "formatted protein name"
    """
    #################################################################
    #             EXTRACT NAMES FROM NAMES EXCEL FILE               #
    #################################################################
    df_names = pd.read_excel(names_excel_path, index_col=0)
    # restrict names dict to only that database
    df_names["acc"] = df_names.index
    df_names["acc_db"] = df_names.acc + "-" + df_names.database
    df_names.set_index("acc_db", inplace=True, drop=False)
    #df_names = df_names.loc[df_names.database == database]
    if style == "shortname [acc-db]":
        df_names["label"] = df_names.shortname + " [" + df_names.acc_db + "]"
    elif style == "shortname [acc]":
        df_names["label"] = df_names.shortname + " [" + df_names.acc + "]"
    else:
        raise ValueError("other styles not implemented")
    namedict = df_names["label"].to_dict()
    return namedict

def pdf_subpath(png_path):
    """Create a subfolder "pdf" where the png is saved.

    Also checks that the path exists.

    Parameters
    ----------
    png_path : str
        Path to png file

    Returns
    -------

    """
    pdf_path = os.path.join(os.path.dirname(png_path), "pdf", os.path.basename(png_path)[:-4] + ".pdf")
    make_sure_path_exists(pdf_path, isfile=True)
    return pdf_path

def drop_redundant_proteins_from_list(df_set, logging):
    """Simply drops the proteins labeled "False" as CD-HIT cluster representatives.

    Relies on the dataframe containing the cdhit_cluster_rep column.

    Parameters
    ----------
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Returns
    -------
    df_set_nonred : pd.DataFrame
        df_set with only non-redundant proteins, based on redundancy settings
    """
    if "no_cdhit_results" in df_set.cdhit_cluster_rep:
        logging.warning("No CD-HIT results were used to remove redundant seq,  but model is being trained anyway.")

    n_prot_initial = df_set.shape[0]
    # create model only using CD-HIT cluster representatives
    df_set_nonred = df_set.loc[df_set.cdhit_cluster_rep != False]
    n_prot_final = df_set_nonred.shape[0]
    logging.info("CDHIT redundancy reduction : n_prot_initial = {}, n_prot_final = {}, n_prot_dropped = {}".format(n_prot_initial, n_prot_final, n_prot_initial - n_prot_final))
    return df_set_nonred

def add_res_num_full_seq_to_df(acc, df, TMD_seq, full_seq, prediction_name, file):
    """

    Parameters
    ----------
    df : pd.DataFrame
    TMD_seq : str
        TMD sequence according to df_set
    full_seq : str
        Full protein sequence according to df_set
    Returns
    -------

    """
    # use regex to get indices for start and end of TMD in seq
    m = re.search(TMD_seq, full_seq)
    if m:
        # convert from python indexing to unprot indexing
        TMD_start = m.start() + 1
        TMD_end = m.end()
        df["res_num_full_seq"] = np.array(range(df.shape[0])) + TMD_start
    else:
        raise IndexError("TMD seq not found in full_seq.\nacc = {}\nTMD_seq = {}\nfull_seq = {}\n"
                         "prediction_name={},file={}".format(acc, TMD_seq, full_seq, prediction_name, file))
    return df

def calculate_identity(sequenceA, sequenceB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned.
    """

    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
    matches = [sa[i] == sb[i] for i in range(sl)]
    seq_id = (100 * sum(matches)) / sl

    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    gap_num = sum([1 for i in range(sl) if (sa[i] == '-' or sb[i] == '-')])
    if gapless_sl > 0:
        gap_id = (100 * sum(matches)) / gapless_sl
    else:
        gap_id = 0
    return (seq_id, gap_id, gap_num)

def join_two_algined_seqences(aligned_A, aligned_B):
    '''
    combine two aligned sequences, take the non-gap residue for the combined seqence
    Parameters
    ----------
    aligned_A
    aligned_B

    Returns
    -------

    '''
    aligned_AB = ""
    for j in range(len(aligned_A)):
        if aligned_A[j] == '-' and aligned_B[j] == '-':
            sys.stdout.write(
                "this homo pair should be considered removed:\n check {},{}".format(aligned_A,aligned_B))
        if aligned_A[j] != '-':
            aligned_AB = aligned_AB + aligned_A[j]
            continue
        if aligned_B[j] != '-':
            aligned_AB = aligned_AB + aligned_B[j]
            continue
    return aligned_AB

def shorten(x):
    '''
    convert 3-letter amino acid name to 1-letter form
    Parameters
    ----------
    x : str , thrre letter aa sequence

    Returns
    y : str, one-letter aa sequence
    -------

    '''
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
         'UNK' : 'U'}
    if len(x) % 3 != 0:
        raise ValueError('Input length should be a multiple of three')

    y = ''
    for i in range(int(len(x)/3)):
            y += d[x[3*i:3*i+3]]
    return y

def Get_Closedist_between_ChianA_ChainB(hashclosedist):
    i = 0
    j = 0
    hashA = {}
    hashB = {}
    closest_dist_arr = []

    jk = ""
    for k, v in sorted(hashclosedist.items()):
        if re.search('NEN', k) or re.search('CEN', k):
            continue
        k = k.split(':')
        k1 = '_'.join(k)
        k2 = '_'.join([k1, str(v)])
        jk = '+'.join([jk,  k2])
    return jk

def add_mutation_missed_residues_with_na(s,acc,database,df):
    acc_combind_feature_file = os.path.join(s['features_folder'],"combined",database,"{}.surr{}.gaps{}.combined_features.csv".format(acc,s["num_of_sur_residues"],s["max_n_gaps_in_TMD_subject_seq"]))
    df_feature = pd.read_csv(acc_combind_feature_file,engine="python",index_col=0)
    not_df_index = [element for element in df_feature["residue_num"].values if element not in df.index.values ]
    for element in not_df_index:
        df.loc[element] = [df_feature.loc[element -1,"residue_name"], np.nan]
    df = df.sort_index()
    return df

def normalise_0_1(arraylike):
    """ Normalise an array to values between 0 and 1.

    The following linear formula is used.
    norm_array = (orig_array - array_min)/(array_max - array_min)

    The use of this simple linear formula allows the normalised data to be "denormalised" later, so long as
    the min and max values of the original array are known.

    Parameters
    ----------
    arraylike : array
        Numpy array (or other arraylike) dataset of floats or ints to be normalised.

    Returns
    -------
    normalised : array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # or, if denormalisation is not necessary
    normalised_array = normalise_0_1(original_array)[0]
    # for further usage examples, see the docstring for denormalise_0_1
    """
    array_min = np.min(arraylike)
    array_max = np.max(arraylike)
    normalised = (arraylike - array_min)/(array_max - array_min)
    # convert to float
    normalised = np.array(normalised).astype(float)
    return normalised, array_min, array_max

def denormalise_0_1(value_or_array, array_min, array_max):
    """ Denormalise a value or array to orig values.

    For use after normalisation between 0 and 1 with the normalise_0_1 function.

    The normalisation formula (normalise_0_1):
        norm_array = (orig_array - array_min)/(array_max - array_min)

    The denormalisation formula (denormalise_0_1):
        denormalised_array = norm_array*(array_max - array_min) + array_min

    Parameters
    ----------
    value_or_array : int, float or arraylike
        Int or float to be denormalised.
        Numpy array (or other arraylike) of data (float, int, etc) to be denormalised.

    Returns
    -------
    normalised : float, or numpy array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    from eccpy.tools import normalise_0_1, denormalise_0_1
    import numpy as np
    original_array = np.linspace(10,130,10)
    original_array[2], original_array[4] = 3, 140
    # normalise original array
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # do stuff to normalised array (e.g., multiply by 0.5)
    normalised_array_halved = normalised_array * 0.5
    # denormalise values to match the equivalents in the original array.
    # Note that the min value (3) was normalised to zero, and was therefore not affected by multiplication.
    normalised_array_halved_denorm = denormalise_0_1(normalised_array_halved, min_, max_)
    # now calculate average values, and check that they match
    norm_array_mean = np.mean(normalised_array)
    norm_array_mean_denormalised = denormalise_0_1(norm_array_mean, min_, max_)
    orig_array_mean = np.mean(original_array)
    # print the two mean values. They should be equal.
    norm_array_mean_denormalised == orig_array_mean
    """
    if isinstance(value_or_array, list):
        raise ValueError('this function accepts arraylike data, not a list. '
                         'Please check data or convert list to numpy array')
    elif isinstance(value_or_array, float):
        denormalised = value_or_array*(array_max - array_min) + array_min
    elif isinstance(value_or_array, np.ndarray):
        denormalised = value_or_array*(array_max - array_min) + array_min
    elif isinstance(value_or_array, pd.Series):
        denormalised = value_or_array*(array_max - array_min) + array_min
    else:
        sys.stdout.write("Unknown datatype. denormalise_0_1 has been given an input that does not appear to be "
              "an int, float, np.ndarray or pandas Series\n"
              "Attempting to process as if it is arraylike.....")
    return denormalised


def normalise_between_2_values(arraylike, min_value, max_value, invert=False):
    """Normalises an array of data between two desired values.

    Any values below min_value will be converted to 0.
    Any values above max_value will be converted to 1.
    Optionally, the normalised array can be inverted, so that the original highest
    values are 0, and the original lowest values are now 1.

    Parameters
    ----------
    arraylike : np.ndarray
        Arraylike original data (numpy array or pandas Series)
    min_value : float
        Desired minimum value for normalisation
    max_value : float
        Desired max value for normalisation
    invert : bool
        If True, normalised data will be inverted (former highest value = 0)

    Returns
    -------
    normalised : np.ndarray
        Normalised array of data

    Usage
    -----
    from eccpy.tools import normalise_between_2_values
    # for array
    orig_array = np.array(range(0, 15))
    norm_array = normalise_between_2_values(orig_array, 3, 10)
    # for pandas Dataframe
    df["norm_data"] = normalise_between_2_values(df["orig_data"], 3, 10)
    """
    # normalise array between min and max values
    normalised = (arraylike - min_value)/(max_value - min_value)
    # replace anything above 1 with 1
    normalised[normalised > 1] = 1
    # replace anything below 0 with 0
    normalised[normalised < 0] = 0
    # if desired, invert the normalised values
    if invert:
        normalised = abs(normalised - 1)
    return normalised

def create_colour_lists():
    '''
    Converts several lists of rgb colours to the python format (normalized to between 0 and 1)
    Returns a dictionary that contains dictionaries of palettes with named colours (eg. TUM blues)
    and also lists of unnamed colours (e.g. tableau20)
    (copied from tlabtools 2016.08.08)
    '''
    output_dict = {}

    matplotlib_150 = list(colors.cnames.values())
    output_dict['matplotlib_150'] = matplotlib_150

    #define colour dictionaries. TUM colours are based on the style guide.
    colour_dicts = {
                    'TUM_colours' : {
                                    'TUMBlue':(34,99,169),
                                    'TUM1':(100,160,200),
                                    'TUM2':(1,51,89),
                                    'TUM3':(42,110,177),
                                    'TUM4':(153,198,231),
                                    'TUM5':(0,82,147)
                                    },
                    'TUM_oranges': {
                        'TUM0': (202, 101, 10),
                        'TUM1': (213, 148, 96),
                        'TUM2': (102, 49, 5),
                        'TUM3': (220, 108, 11),
                        'TUM4': (247, 194, 148),
                        'TUM5': (160, 78, 8)
                    },
                    'TUM_accents' : {
                                    'green':(162,183,0),
                                    'orange':(227,114,34),
                                    'ivory':(218,215,203),
                                    }
                    }

    #convert the nested dicts to python 0 to 1 format
    for c_dict in colour_dicts:
        for c in colour_dicts[c_dict]:
            #define r, g, b as ints
            r, g, b = colour_dicts[c_dict][c]
            #normalise r, g, b and add to dict
            colour_dicts[c_dict][c] = (r / 255., g / 255., b / 255.)
        #add normalised colours to output dictionary
        output_dict[c_dict] = colour_dicts[c_dict]

    #define colour lists
    colour_lists = {
                    'tableau20' : [
                                 (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
                                    ],
                    'tableau20blind' : [
                                         (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
                                         (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
                                         (255, 188, 121), (207, 207, 207)
                                          ]
                    }
    #normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list])/255.
        colour_array_tup = tuple(map(tuple,colour_array))
        colour_lists[rgb_list] = colour_array_tup
        #add normalised colours to output dictionary
        output_dict[rgb_list] = colour_lists[rgb_list]
    #create a mixed blue/grey colour list, with greys in decreasing darkness
    TUM_colours_list_with_greys = []
    grey = 0.7
    for c in colour_dicts['TUM_colours'].values():
        TUM_colours_list_with_greys.append('%0.2f' % grey)
        TUM_colours_list_with_greys.append(c)
        grey -= 0.1
    output_dict['TUM_colours_list_with_greys'] = TUM_colours_list_with_greys

    output_dict['HTML_list01'] = ['#808080', '#D59460', '#005293', '#A1B11A', '#9ECEEC', '#0076B8', '#454545', "#7b3294", "#c2a5cf", "#008837", "#a6dba0"]
    return output_dict

def get_free_space(folder, format="MB"):
    """
        Return folder/drive free space
    """
    fConstants = {"GB": 1073741824,
                  "MB": 1048576,
                  "KB": 1024,
                  "B": 1
                  }
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return (int(free_bytes.value/fConstants[format.upper()]), format)
    else:
        return (int(os.statvfs(folder).f_bfree*os.statvfs(folder).f_bsize/fConstants[format.upper()]), format)

def create_regex_string(inputseq):
    ''' adds '-*' between each aa or nt/aa in a DNA or protein sequence, so that a particular
    aligned sequence can be identified via a regex search, even if it contains gaps
    inputseq : 'LQQLWNA'
    output   : 'L-*Q-*Q-*L-*W-*N-*A'
    '''
    search_string = ''
    for letter in inputseq:
        letter_with_underscore = letter + '-*'
        search_string += letter_with_underscore
    return search_string[:-2]

def convert_truelike_to_bool(input_item, convert_int=False, convert_float=False, convert_nontrue=False):
    """Converts true-like values ("true", 1, True", "WAHR", etc) to python boolean True.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "true", 1, "WAHR" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "1.0" will be converted to True
    convert_nontrue : bool
        If True, the output for input_item not recognised as "True" will be False.
        If True, the output for input_item not recognised as "True" will be the original input_item.

    Returns
    -------
    return_value : True, or input_item
        If input_item is True-like, returns python bool True. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_truelike_to_bool("true")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_truelike_to_bool)
    """
    list_True_items = [True, 'True', "true","TRUE","T","t",'wahr', 'WAHR', 'prawdziwy', 'verdadeiro', 'sann', 'istinit',
                       'veritable', 'Pravda', 'sandt', 'vrai', 'igaz', 'veru', 'verdadero', 'sant', 'gwir', 'PRAWDZIWY',
                       'VERDADEIRO', 'SANN', 'ISTINIT', 'VERITABLE', 'PRAVDA', 'SANDT', 'VRAI', 'IGAZ', 'VERU',
                       'VERDADERO', 'SANT', 'GWIR', 'bloody oath', 'BLOODY OATH', 'nu', 'NU','damn right','DAMN RIGHT']

    # if you want to accept 1 or 1.0 as a true value, add it to the list
    if convert_int:
        list_True_items += ["1"]
    if convert_float:
        list_True_items += [1.0, "1.0"]
    # check if the user input string is in the list_True_items
    input_item_is_true = input_item in list_True_items
    # if you want to convert non-True values to "False", then nontrue_return_value = False
    if convert_nontrue:
        nontrue_return_value = False
    else:
        # otherwise, for strings not in the True list, the original string will be returned
        nontrue_return_value = input_item
    # return True if the input item is in the list. If not, return either False, or the original input_item
    return_value = input_item_is_true if input_item_is_true == True else nontrue_return_value
    # special case: decide if 1 as an integer is True or 1
    if input_item == 1:
        if convert_int == True:
            return_value = True
        else:
            return_value = 1
    return return_value

def convert_falselike_to_bool(input_item, convert_int=False, convert_float=False):
    """Converts false-like values ("false", 0, FALSE", "FALSCH", etc) to python boolean False.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "FALSE", 0, "FALSCH" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "0.0" will be converted to True

    Returns
    -------
    return_value : False, or input_item
        If input_item is False-like, returns python bool False. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_falselike_to_bool("false")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_falselike_to_bool)
    """
    list_False_items = [False, "False", "false", "FALSE", "F", "f", "falsch", "FALSCH", "valse", "lažna", "fals",
                        "NEPRAVDA", "falsk", "vals", "faux", "pa vre", "tsis tseeb", "hamis", "palsu", "uongo", "ngeb",
                        "viltus", "klaidinga", "falz", "falso", "USANN", "wartosc false", "falošné", "falskt", "yanlis",
                        "sai", "ffug", "VALSE", "LAŽNA", "FALS", "FALSK", "VALS", "FAUX", "PA VRE", "TSIS TSEEB",
                        "HAMIS", "PALSU", "UONGO", "NGEB", "VILTUS", "KLAIDINGA", "FALZ", "FALSO", "WARTOSC FALSE",
                        "FALOŠNÉ", "FALSKT", "YANLIS", "SAI", "FFUG"]

    # if you want to accept 0 or 0.0 as a false value, add it to the list
    if convert_int:
        list_False_items += [0, "0"]
    if convert_float:
        list_False_items += [0.0,"0.0"]
    # return boolean False if the input item is in the list. If not, return the original input_item
    return_value = False if input_item in list_False_items else input_item

    return return_value

class HardDriveSpaceException(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        # name changed to allow p-r( to be unique to the print function
        canonical_string_representation = repr
        return canonical_string_representation(self.parameter)

class Log_Only_To_Console(object):
    def __init__(self):
        pass
    def info(self, message):
        sys.stdout.write("\n{}".format(message))
    def warning(self, message):
        sys.stdout.write("\n{}".format(message))
    def critical(self, message):
        sys.stdout.write("\n{}".format(message))



def calc_rand_overlap_DEPRECATED_METHOD(sequence_length, sample_size):
    """Calculate expected random overlap between two random selections from a sample.

    You have a bowl of 20 balls, either blue or red.
    sample_size = 1
     - There is only one red ball. You only pick out one ball.
     - you do this 1000 times
     - the average number of red (correct) balls 0.05
    sample_size = 2
     - there are 2 red balls. You can pick out 2 balls.
     - you do this 1000 times
     - what is the average number of red (correct) balls?
     [This function calculates the answer, which is 0.2, corresponding to a red ball 10% of the time (0.2/2))
    sample_size = 10
     - there are 10 red balls. You can pick out 10 balls.
     - you do this 1000 times
     - what is the average number of red (correct) balls?
     [This function calculates the answer, which is 5, corresponding to a red ball 50% of the time (5/10))
     The answer depends on the number of balls in the bowl. If there are 24 balls, the average number of red balls is 4.16,
     corresponding to a red ball 41.6% of the time (4.16/10))

    This is a mathematical question related to sampling without replacement.

    Parameters
    ----------
    sequence_length : int
        Number of samples in total, from which a subset is randomly selected.
        E.g. a bowl with 20 numbered balls.
        In our case, TMD_length.
    sample_size : int
        Size of the selected sample (e.g. 3, for 3 balls taken from the bowl)
        In our case, number of interface residues considered.

    Returns
    -------
    inter_num_random : float
        Expected overlap seen by random chance.
        E.g., when 5 balls from 20 are randomly selected without replacement.
        The expected overlap between two random selections is 1.25

    Usage
    -----
    from thoipapy.utils import calc_rand_overlap
    # size of bowl (in our case, transmembrane domain length)
    TMD_length = 22
    # number of balls withdrawn from bowl (in our case, top 4 residues predicted)
    sample_size = 4
    # expected overlap is the number predicted to be correct, simply by random chance
    expected_overlap = calc_rand_overlap(TMD_length, sample_size)
    percentage_randomly_correct = expected_overlap / sample_size * 100
    """
    # start the [[write description]] at 0
    inter_num_random = 0
    # iterate through [[write description]]
    for random_inter_num_ober in range(1, sample_size + 1):
        # random_inter_num_ober = j

        # write description
        random_non_inter_num = sequence_length - sample_size
        # write description
        variable_1 = comb(sample_size, random_inter_num_ober)
        # write description
        variable_2 = comb(random_non_inter_num, sample_size - random_inter_num_ober)
        # write description
        variable_3 = comb(sequence_length, sample_size)
        inter_num_random = inter_num_random + variable_1 * variable_2 / variable_3 * random_inter_num_ober
        # inter_num_random = inter_num_random + comb(sample_size, random_inter_num_ober) * comb(random_non_inter_num, sample_size - random_inter_num_ober) / comb(tm_len, sample_size) * random_inter_num_ober

    return inter_num_random

def rename_features(df_features_single_protein):
    """rename selected features.

    coev_all_top4_MI -> MItop4mean
    coev_all_top4_DI -> DItop4mean
    coev_all_top8_MI -> MItop8mean
    coev_all_top8_DI -> DItop8mean
    coev_i1_MI -> MI1mean
    coev_i1_DI -> DI1mean
    coev_i3_MI -> MI3mean
    coev_i3_DI -> DI3mean
    coev_i4_MI -> MI4mean
    coev_i4_DI -> DI4mean
    coev_all_max_MI -> MImax
    coev_all_max_DI -> DImax
    coev_i1-i4_max_MI -> MI4max
    coev_i1-i4_max_DI -> DI4max
    CumMI4 -> MI4cum
    CumDI4 -> DI4cum
    highest_face_MI -> MI_highest_face
    highest_face_DI -> DI_highest_face


    Parameters
    ----------
    df_features_single_protein : pd.DataFrame
        Dataframe with all features for a protein
        Index : range index
        Columns : "residue_num", "residue_name", "Entropy", etc
    """

    # rename features
    df_features_single_protein = df_features_single_protein.rename(
        columns={"coev_all_top4_MI": "MItop4mean", "coev_all_top4_DI": "DItop4mean", "coev_all_top8_MI": "MItop8mean", "coev_all_top8_DI": "DItop8mean",
                 "coev_i1_MI": "MI1mean", "coev_i1_DI": "DI1mean", "coev_i3_MI": "MI3mean","coev_i3_DI": "DI3mean",
                 "coev_i4_MI": "MI4mean", "coev_i4_DI": "DI4mean", "coev_all_max_MI": "MImax","coev_all_max_DI": "DImax",
                 "coev_i1-i4_max_MI": "MI4max", "coev_i1-i4_max_DI": "DI4max", "CumMI4": "MI4cum","CumDI4": "DI4cum",
                 "highest_face_MI": "MI_highest_face", "highest_face_DI": "DI_highest_face"})

    return df_features_single_protein