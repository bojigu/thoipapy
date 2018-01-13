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
import plotly

class Command(object):
    '''
    subprocess for running shell commands in win and linux
    This will run commands from python as if it was a normal windows console or linux terminal.
    taken from http://stackoverflow.com/questions/17257694/running-jar-files-from-python)'
    '''
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            #logging.info('Thread started')
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #self.process.communicate()
            stdout, stderr = self.process.communicate() # from http://stackoverflow.com/questions/14366352/how-to-capture-information-from-executable-jar-in-python
            # Thus far, SIMAP has only ever given java faults, never java output. Don't bother showing.
            # if the console prints anything longer than 5 characters, log it
            if len(stderr.decode("utf-8")) > 5:
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
                print("{} could not be deleted".format(xml_renamed))

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
        print("{} could not be deleted".format(blast_xml_file))
    try:
        os.remove(xml_txt)
    except:
        print("{} could not be deleted".format(xml_txt))

def setup_biopol_plotly(username, api_key):
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

def add_res_num_full_seq_to_df(acc, df, TMD_seq, full_seq):
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
        raise IndexError("TMD seq not found in full_seq.\nacc = {}\nTMD_seq = {}\nfull_seq = {}".format(acc, TMD_seq, full_seq))
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
