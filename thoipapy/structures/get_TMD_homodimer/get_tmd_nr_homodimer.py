import os
import sys
import thoipapy
import thoipapy.utils as utils
import urllib.request
import xml.etree.ElementTree as ET
import csv
import pandas as pd
from Bio import pairwise2
import shutil
import glob
import gzip
import re
import math
import numpy as np
import tarfile
from matplotlib import pyplot as plt
from scipy.stats import ttest_ind
from scipy.optimize import leastsq
from pytoxr.mathfunctions import sine_perfect_helix, residuals
from thoipapy.utils import create_colour_lists
from korbinian.utils import pn, aaa

def download_xml_get_alphahelix_get_homo_pair(s, logging):
    """
    download xml for all tmps
    Parameters
    ----------
    s : dict. setting file dictionary
    logging

    Returns
    -------

    """
    logging.info('****start to download xml file for all alpha helix proteins from database pdbtm*********')
    thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    pdbtm_alpha_list = os.path.join(thoipapy_module_path, "get_TMD_homodimer", "pdbtm_alpha_list","pdbtm_alpha_20170616.list")
    pdbtm_xml_path = os.path.join(s["pdbtm_homodimer_folder"],"xml")
    with open(pdbtm_alpha_list , 'r' ) as f:
        for line in f:
            pdb_id = line[0:4]
            #get_xml_file(pdb_id, pdbtm_xml_path, logging)
            #get_tmd_alpha_helix_infor(pdb_id, pdbtm_xml_path, logging)
            get_inter_pair_homodimer(pdb_id, pdbtm_xml_path, logging)


def download_trpdb_calc_inter_rr_pairs(s, logging):
    tm_homo_pair_file_lists = glob.glob(os.path.join(s["pdbtm_homodimer_folder"],"xml","*.interhelix.csv" ))
    for tm_homo_pair_file in tm_homo_pair_file_lists:
        pdb_id = tm_homo_pair_file.split('\\')[-1][0:4]
        print(pdb_id)
        tr_pdb_path = os.path.join(s["pdbtm_homodimer_folder"],"xml")
        utils.make_sure_path_exists(tr_pdb_path)
        pdbtm_trpdb_file = os.path.join(tr_pdb_path, "{}.trpdb.gz".format(pdb_id))
        if not os.path.isabs(pdbtm_trpdb_file):
            get_trpdb_file(pdb_id, pdbtm_trpdb_file, logging)
        else:
            inter_def_cutoff = 4.5
            #if pdb_id == "1jb0":
            multiple_add_closedist_rr_inter_pairs(s,pdbtm_trpdb_file,tm_homo_pair_file,inter_def_cutoff)

def create_redundant_interact_homodimer_rm_shorttm(s, logging):
    tm_homo_pair_file_lists = glob.glob(os.path.join(s["pdbtm_homodimer_folder"], "xml", "*.interhelix.csv"))
    redundant_interact_homodimer_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","redundant_interact_homodimer.csv")
    redundant_interact_homodimer_targz_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","redundant_interact_homodimer.tar.gz")
    utils.make_sure_path_exists(redundant_interact_homodimer_file, isfile=True)
    redundant_interact_homodimer_file_handle =open(redundant_interact_homodimer_file,'w')
    dfc = pd.DataFrame()
    for tm_homo_pair_file in tm_homo_pair_file_lists:
        pdb_id = tm_homo_pair_file.split('\\')[-1][0:4]
        #if pdb_id == "1ap9":
        df_homo = pd.read_csv(tm_homo_pair_file, engine = "python", index_col=0)
        cols = [c for c in df_homo.columns if c[:7] != 'Unnamed']
        df_homo = df_homo[cols]
        df_homo.rename(columns={"interpari_4.5": "interpair_4.5"}, inplace=True)
        df_homo_matrix = df_homo.as_matrix()
        index_drop = []
        for i in range(df_homo.shape[0]):
            ## filter1: remove all pairs without interact pair
            if pd.isnull(df_homo.iloc[i]["interpair_4.5"]):
                index_drop.append(i)
                continue
            ##filter2: remove all pairs that with tm length <15 or >30
            if len(df_homo.iloc[i]["aligned_AB"]) <15 or len(df_homo.iloc[i]["aligned_AB"]) > 30:
                index_drop.append(i)
        df_homo.drop(df_homo.index[index_drop], inplace=True)
        if dfc.empty:
            dfc = df_homo
        else:
            dfc = pd.concat([dfc, df_homo])
    dfc.reset_index(drop=True,inplace=True)
    dfc.to_csv(redundant_interact_homodimer_file_handle)
    redundant_interact_homodimer_file_handle.close()
    if os.path.isfile(redundant_interact_homodimer_file):
        with tarfile.open(redundant_interact_homodimer_targz_file, mode='w:gz') as tar:
            # add the files to the compressed tarfile
            tar.add(redundant_interact_homodimer_file, arcname=os.path.basename(redundant_interact_homodimer_file))
    ##delete the reduandant inter homodimer csv file, since it was saved in targz file
    try:
        os.remove(redundant_interact_homodimer_file)
    except:
        print("{} could not be deleted".format(redundant_interact_homodimer_file))
    sys.stdout.write("the header of the redundant inter homodimers I: {}\n".format(dfc.head()))
    sys.stdout.flush()
    sys.stdout.write("the dimension of the redundant inter homodimers is: {} \n".format(dfc.shape))
    sys.stdout.flush()

def extract_crystal_resolv035_interact_pairs_and_create_fasta_file(s,logging):
    redundant_interact_homodimer_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","redundant_interact_homodimer.csv")
    crystal035_redundant_interact_homodimer_file = os.path.join(s["pdbtm_homodimer_folder"], "xml", "homodimer",
                                                     "crystal3.5_redundant_interact_homodimer.csv")
    redundant_interact_homodimer_targz_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","redundant_interact_homodimer.tar.gz")
    crystal0_35_redundant_interact_homodimer_fasta_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","crystal3.5_redundant_interact_homodimer.fasta")
    cdhit60_crystal035_clr_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","cdhit_0.6_crystal3.5_redundant_interact_homodimer.clstr")
    cdhit60_crystal035_fasta_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","cdhit_0.6_crystal3.5_redundant_interact_homodimer.fasta")
    cdhit60_nr_represent228_pair = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","cdhit_0.6_nr_represent_228_interpair.csv")
    cdhit60_nr_represent228_pair_handle = open(cdhit60_nr_represent228_pair, 'w')
    crystal0_35_redundant_interact_homodimer_fasta_file_handle = open(crystal0_35_redundant_interact_homodimer_fasta_file,'w')
    with tarfile.open(redundant_interact_homodimer_targz_file, 'r:gz') as tar:
        tar.extractall(os.path.dirname(redundant_interact_homodimer_targz_file))
    df_rd = pd.read_csv(redundant_interact_homodimer_file)
    resolve_list = []
    for i in range(df_rd.shape[0]):
        pdb_id = str(df_rd.iloc[i]["pdb_id"])
        resolv = get_pdb_experimental_resolution(s, pdb_id)
        resolve_list.append(resolv)
    df_rd["resolv"] = resolve_list
    print(df_rd.shape)
    df_rd =df_rd[df_rd["resolv"].astype(float) <= 3.5]
    print(df_rd.shape)

    #save tmseq into one fasta file for cd-hit use
    for i in range(df_rd.shape[0]):
        fasta_header = ">" + df_rd.iloc[i]["pdb_id"] + "_" + df_rd.iloc[i]["tm_numA"] + "_" + df_rd.iloc[i]["tm_numB"]
        fasta_tmseq = df_rd.iloc[i]["aligned_AB"]
        crystal0_35_redundant_interact_homodimer_fasta_file_handle.write(fasta_header + "\n" +fasta_tmseq + "\n")
    crystal0_35_redundant_interact_homodimer_fasta_file_handle.close()

    ##add cdhit cluster number to dataframe
    cluster_number_list = []
    dict_cluster_pair = get_interpair_cdhist_cluster_number(df_rd,cdhit60_crystal035_clr_file)
    for i in range(df_rd.shape[0]):
        inter_pair= df_rd.iloc[i]["pdb_id"] + "_" + df_rd.iloc[i]["tm_numA"] + "_" + df_rd.iloc[i]["tm_numB"]
        cluster_number = dict_cluster_pair[inter_pair]
        cluster_number_list.append(cluster_number)
    df_rd["cdhit0.6_clusternum"] = cluster_number_list

    ###get toatl rr inter pair number and the non-self inter pair rr number and add to dataframe
    non_self_inter_pair_num_list, all_inter_pair_num_list = get_inter_pair_rr_num(df_rd)
    df_rd["total_intpair_num"] = all_inter_pair_num_list
    df_rd["nonself_intpair_num"] = non_self_inter_pair_num_list
    print(df_rd[["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]])
    a = df_rd["cdhit0.6_clusternum"].unique()
    # for j in sorted(a):
    #     #print(df_rd[df_rd["cdhit0.6_clusternum"].astype(int) == j,["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]])
    #     print(df_rd[["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]].loc[df_rd['cdhit0.6_clusternum'] == j])

    ##update redundant_interact_homodimer_file only keep crystal with resolution less than 0.35, and add 'cdhit60_clusternum
    crystal035_redundant_interact_homodimer_file_handle = open(crystal035_redundant_interact_homodimer_file,'w')
    df_rd.to_csv(crystal035_redundant_interact_homodimer_file_handle)
    crystal035_redundant_interact_homodimer_file_handle.close()

    ##save cdhit60 represent 228 interpairs
    df_new = get_cdhit60_represent_interpair_df(df_rd, cdhit60_crystal035_fasta_file)

    cols = [c for c in df_new.columns if c[:7] != 'Unnamed']
    df_new = df_new[cols]
    df_new.to_csv(cdhit60_nr_represent228_pair_handle)
    cdhit60_nr_represent228_pair_handle.close()

    ##create set228 file used for thoipapy , the inter pair with non-self intor rr pair number bigger than 1
    set_file = os.path.join(s["sets_folder"],'set228_crystal.xlsx')
    create_inter_rr_number_bt1_set228_file(df_new, set_file)

def create_multiple_bind_closedist_file(s,df_set, logging):

    non_redundant_homodimer_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","homodimer","cdhit_0.6_nr_represent_228_interpair.csv")
    df_homo = pd.read_csv(non_redundant_homodimer_file, engine="python")
    cols = [c for c in df_homo.columns if c[:7] != "Unnamed"]
    df_homo = df_homo[cols]
    inter_pair_max = s['inter_pair_max']
    for i in range(df_homo.shape[0]):
        create_single_bind_closedist_file(s, df_homo,i, inter_pair_max)

def create_single_bind_closedist_file(s, df_homo,i, inter_pair_max):
        pdb_id_chain = df_homo.iloc[i]['pdb_id'] + df_homo.iloc[i]['tm_numA']
        row = []
        #bind_closedist_file = os.path.join(s['structure_bind'], "crystal", "{}.{}pairmax.bind.closedist.csv".format(pdb_id_chain,inter_pair_max))
        bind_closedist_file = os.path.join(s['structure_bind'], "crystal", "{}.{}pairmax.bind.closedist.csv".format(pdb_id_chain,inter_pair_max))
        tm_seq = df_homo.iloc[i]['aligned_AB']
        full_seq = df_homo.iloc[i]['full_seqA']
        aligned_AB_startA = df_homo.iloc[i]['full_seqA'].index(tm_seq) + 1
        aligned_AB_startB = df_homo.iloc[i]['full_seqB'].index(tm_seq) + 1
        closedistA = df_homo.iloc[i]['closedist_A']
        closedistB = df_homo.iloc[i]['closedist_B']
        closedistA_dict = {}
        closedistB_dict = {}
        closedistAB_dict = {}
        for j in range(1, len(closedistA.split('+'))):
            closedistA_dict[int(closedistA.split('+')[j].split('_')[1]) - aligned_AB_startA] = \
            float(closedistA.split('+')[j].split('_')[4])
        for j in range(1, len(closedistB.split('+'))):
            closedistB_dict[int(closedistB.split('+')[j].split('_')[1]) - aligned_AB_startB] = \
            float(closedistB.split('+')[j].split('_')[4])
        closedistA_list = [int(closedistA.split('+')[i].split('_')[1]) - aligned_AB_startA for i in
                           range(1, len(closedistA.split('+')))]
        closedistB_list = [int(closedistB.split('+')[i].split('_')[1]) - aligned_AB_startB for i in
                           range(1, len(closedistB.split('+')))]
        closedistAB_list = set(closedistA_list + closedistB_list)
        for j in closedistAB_list:
            if j in closedistA_dict and j not in closedistB_dict:
                closedistAB_dict[j] = closedistA_dict[j]
            if j in closedistB_dict and i not in closedistA_dict:
                closedistAB_dict[j] = closedistB_dict[j]
            if j in closedistA_dict and j in closedistB_dict:
                if closedistA_dict[j] <= closedistB_dict[j]:
                    closedistAB_dict[j] = closedistA_dict[j]
                else:
                    closedistAB_dict[j] = closedistB_dict[j]

        interpair = df_homo.iloc[i]['interpair_4.5']
        interpair_dict = {}
        interpair_dict_new = {}
        for j in range(1, len(interpair.split('+'))):
            if (str(int(interpair.split('+')[j].split('_')[5]) - aligned_AB_startB +1) +'_'+
                str(int(interpair.split('+')[j].split('_')[1]) - aligned_AB_startA +1)
                           ) not in interpair_dict:
                interpair_dict[str(int(interpair.split('+')[j].split('_')[1]) - aligned_AB_startA + 1) + '_' +
                           str(int(interpair.split('+')[j].split('_')[5]) - aligned_AB_startB + 1)] = \
            interpair.split('+')[j].split('_')[8]
        interpair_num = len(interpair_dict)
        if interpair_num <= inter_pair_max:
            interpair_dict_new = interpair_dict
        else:
            n = 0
            sorted_keys = sorted(interpair_dict, key=interpair_dict.get, reverse=True)
            for r in sorted_keys:
                if n < inter_pair_max:
                    interpair_dict_new[r] = interpair_dict[r]
                    n = n + 1
        tm_inter_list = []
        for key, value in interpair_dict_new.items():
            tm_inter_list.append(key.split('_')[0])
            tm_inter_list.append(key.split('_')[1])
        tm_inter_list = set(tm_inter_list)
        for j in range(len(tm_seq)):
            if j in closedistAB_dict:
                row.append([j + 1, tm_seq[j], '1' if str(j + 1) in tm_inter_list else '0', closedistAB_dict[j]])
            else:
                print("residue{}not in {} closedist".format(j, pdb_id_chain))
                if j - 1 in closedistAB_dict:
                    closedistAB_dict[j] = closedistAB_dict[j - 1]
                if j + 1 in closedistAB_dict:
                    closedistAB_dict[j] = closedistAB_dict[j + 1]
                row.append([j + 1, tm_seq[j], '1' if str(j + 1) in tm_inter_list else '0', closedistAB_dict[j]])
        closest_dist_df = pd.DataFrame.from_records(row, columns=["residue_num", "residue_name", 'bind', "closedist"])
        closest_dist_df.set_index('residue_num', drop=True, inplace=True)
        closest_dist_df = closest_dist_df.sort_index()
        closest_dist_df.to_csv(bind_closedist_file)

def return_set_id_list(set_file):
    df_set = pd.read_excel(set_file)
    id_list = df_set.loc[:, "acc"].values
    return id_list


def get_pivot_table_coev_data(s, i, XI, df_set):
    acc = df_set.loc[i, "acc"]
    database = df_set.loc[i, "database"]
    TMD_start = int(df_set.loc[i, "TMD_start"])
    TMD_end = int(df_set.loc[i, "TMD_end"])
    freecontact_file = os.path.join(s["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

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

    dfp = df.pivot_table(index="n1", columns="n2", values=XI)

    """ asymmetrical pivoted data

         235       236       237       238       239       240 ...        252       253       254       255       256       
    n1                                                     ...                                                              
    235   0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523 
    236        NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577 
    237        NaN       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414 
    238        NaN       NaN       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822 
    """
    # get full list of residues
    position_list = range(TMD_start, TMD_end + 1)
    dfp = dfp.reindex(index=position_list, columns=position_list)
    # put data on both sides of the table for easy indexing
    for col in dfp.columns:
        start = col + 1
        dfp.loc[start:, col] = dfp.loc[col, start:]

    return dfp

def calc_coev_vs_res_dist(s, df_set, logging):
    """Calculate mean coevolution scores for each residue distance

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

    logging.info('calc_coev_vs_res_dist starting')
    coev_vs_res_dist_xlsx = os.path.join(s["set_results_folder"], "{}_coev_vs_res_dist.xlsx".format(s["setname"]))
    writer = pd.ExcelWriter(coev_vs_res_dist_xlsx)

    for XI in ["MI", "DI"]:
        nested_coev_dist_dict = {}
        nested_Cterm_dist_dict = {}

        for i in df_set.index:
            sys.stdout.write(".")
            sys.stdout.flush()
            acc_db = df_set.loc[i, "acc_db"]
            TMD_start = int(df_set.loc[i, "TMD_start"])
            TMD_end = int(df_set.loc[i, "TMD_end"])

            dfp = get_pivot_table_coev_data(s, i, XI, df_set)

            """dfp pivot table has a symmetric coevolution values between all residues
            
                    307       308       309       310       311       312       313       314       315       316    ...          326       327       328       329       330       331       332       333       334       335
            n1                                                                                                         ...                                                                                                       
            307       NaN  0.233388  0.251910  0.257193  0.365270  0.468933  0.313458  0.253943  0.278989  0.297606    ...     0.206634  0.153770  0.271118  0.364004  0.185186  0.286575  0.166321  0.313355  0.269962  0.307910
            308  0.233388       NaN  0.194896  0.290690  0.230827  0.300403  0.371423  0.149162  0.250657  0.283342    ...     0.130392  0.135035  0.317070  0.266557  0.134991  0.244770  0.164932  0.211624  0.185211  0.155684
            309  0.251910  0.194896       NaN  0.376993  0.401137  0.480202  0.313931  0.298846  0.336291  0.317149    ...     0.253053  0.229490  0.359081  0.366537  0.203667  0.264654  0.221240  0.373255  0.300027  0.240920
            310  0.257193  0.290690  0.376993       NaN  0.563525  0.651108  0.411248  0.353645  0.366177  0.455358    ...     0.264179  0.280039  0.398423  0.492760  0.226311  0.401377  0.242655  0.401358  0.316718  0.305772
            311  0.365270  0.230827  0.401137  0.563525       NaN  0.686111  0.367423  0.446554  0.457673  0.488534    ...     0.327568  0.286183  0.489552  0.468001  0.259249  0.447306  0.250619  0.455895  0.444784  0.345875

            """

            # set the "buffer", which controls the distance between the residues
            # will also be used to buffer the dataframe
            buffer = 20
            # extend axes by adding buffer
            rng = range(TMD_start - buffer, TMD_end + buffer)
            dfp = dfp.reindex(index=rng, columns=rng)

            coev_dist_dict = {}
            N_coev_dist_dict = {}

            # iterate through each distance, b, in the buffer
            for b in range(1, buffer + 1):

                ###############################################
                #         Method 1 : mean of i+1 and i-1      #
                ###############################################
                single_dist_list = []
                for i in range(TMD_start, TMD_end + 1):
                    seln = dfp.loc[i, [i - b, i + b]]
                    """ The selection contains the two values, whose mean is taken.
                    Often, as b increases, one of the values is a nan.
                    [0.582511, nan]
                    [nan, 0.612737]
                    [nan, 0.5906779999999999]
                    """
                    mean_ = seln.mean()
                    if not np.isnan(mean_):
                       single_dist_list.append(mean_)
                    #single_dist_list.append(mean_)

                    # if np.isnan(mean_):
                    #     sys.stdout.write("({}.{}.{}.{})".format(b,i,seln, mean_))
                    #     sys.stdout.flush()

                ###############################################
                #     Method 2 : collect only i+1, i+2, etc   #
                ###############################################
                i_plus_b_coev_list = []
                for i in range(TMD_start, TMD_end - b + 1):
                    #i_plus_b is always a single coevolution value, between residue i, and residue i+b
                    i_plus_b = dfp.loc[i, i + b]
                    if not np.isnan(i_plus_b):
                        i_plus_b_coev_list.append(i_plus_b)
                    else:
                        # there should be no nans here, but check anyway.
                        sys.stdout.write("*")

                # get mean for both methods, for this TMD
                mean_for_this_dist = np.mean(single_dist_list)
                mean_Cterm = np.mean(i_plus_b_coev_list)

                # add mean for this distance(b) to summary dictionary
                if not np.isnan(mean_for_this_dist):
                    coev_dist_dict[b] = float("{:.03f}".format(mean_for_this_dist))
                if not np.isnan(mean_Cterm):
                    N_coev_dist_dict[b] = float("{:.03f}".format(mean_Cterm))

            # add summary values for this TMD to the output dicts
            nested_coev_dist_dict[acc_db] = coev_dist_dict
            nested_Cterm_dist_dict[acc_db] = N_coev_dist_dict

        # save output dicts with all TMD values to excel
        pd.DataFrame(nested_coev_dist_dict).to_excel(writer, sheet_name="coev_{}".format(XI))
        pd.DataFrame(nested_Cterm_dist_dict).to_excel(writer, sheet_name="C_{}".format(XI))

    writer.close()
    logging.info('calc_coev_vs_res_dist finished')


def plot_coev_vs_res_dist(s, logging):
    """Plot figure comparing coevolution values against residue distance.

    Uses excel output file from calc_coev_vs_res_dist.

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
    logging.info('plot_coev_vs_res_dist starting')
    plt.rcParams["font.family"] = "Verdana"
    plt.rcParams["font.family"] = "Verdana"
    colour_dict = create_colour_lists()
    blue1 = colour_dict["TUM_colours"]['TUM1']
    blue5 = colour_dict["TUM_colours"]['TUM5']
    TUMblue = colour_dict["TUM_colours"]['TUMBlue']

    fontsize = 9
    coev_vs_res_dist_xlsx = os.path.join(s["set_results_folder"], "{}_coev_vs_res_dist.xlsx".format(s["setname"]))

    fig, ax = plt.subplots(figsize=(4.5, 3.42))


    #######################################################################################################
    #                                                                                                     #
    #                                    add alpha helical sine wave                                      #
    #                                                                                                     #
    #######################################################################################################
    # [  0. ,   3.6,   7.2,  10.8,  14.4,  18. ,  21.6] for desired length
    x_helical = np.ogrid[0:23.4:14j]
    # [1,0, etc to train fitted sine as desired
    y_helical = [1,0]*7
    # fit to a perfect helix starting at 0 using leastsq
    sine_constants_guess_perfhelix = [1.575, 0.5]
    sine_constants_perfhelix1, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = leastsq(residuals,
                                                                        sine_constants_guess_perfhelix,
                                                                        args=(sine_perfect_helix,x_helical,y_helical),
                                                                        full_output=1)

    # create smooth sine curve x-points
    x_rng = np.linspace(0, 20, 200)
    logging.info("sine_constants_perfhelix1 : {}".format(sine_constants_perfhelix1))
    # adjust the height as desired
    sine_constants_fixed_centre = (sine_constants_perfhelix1[0], 0)
    # fit and plot
    yvalues_fitted_perfhelix1 = sine_perfect_helix(sine_constants_fixed_centre, x_rng)
    ax.plot(x_rng, yvalues_fitted_perfhelix1, color="0.8", linestyle="--", label=r"$\alpha$-helical periodicity")

    # only to show that heptad periodicity is not desirable
    plot_heptad_periodicity = False
    if plot_heptad_periodicity:
        # test of heptad motif
        x = range(0,21)
        y = np.array([1,0,0,1,1,0,1] * 3)/4
        ax.plot(x, y, color="r", alpha=0.5, linestyle=":", label="heptad periodicity")

    #######################################################################################################
    #                                                                                                     #
    #                                     add mean MI and DI values                                       #
    #                                                                                                     #
    #######################################################################################################

    # tab is "coev_" or "C_"
    # this corresponds to Method 1 or Method 2 above
    # they both give similar values, but i+b method is less likely to count values twice and is preferred
    excel_tab = "C_"
    #excel_tab = "coev_"

    # plot MI on secondary axis
    ax2 = ax.twinx()
    df = pd.read_excel(coev_vs_res_dist_xlsx, sheetname="{}MI".format(excel_tab))
    mean_ser = df.mean(axis=1)
    mean_ser.plot(ax=ax2, label="mutual information (MI)", fontsize=fontsize, color=colour_dict["TUM_accents"]['orange'])

    df = pd.read_excel(coev_vs_res_dist_xlsx, sheetname="{}DI".format(excel_tab))
    mean_ser = df.mean(axis=1)
    mean_ser.plot(ax=ax, label="direct information (DI)", fontsize=fontsize, color=TUMblue)

    # axis colours
    ax2.tick_params("y", colors=colour_dict["TUM_accents"]['orange'])
    ax.tick_params("y", colors=TUMblue)
    # axis labels
    ax.set_ylabel("mean DI coevolution score", labelpad=-3, fontsize=fontsize, color = TUMblue)
    ax2.set_ylabel("mean MI coevolution score", labelpad=1, fontsize=fontsize, color=colour_dict["TUM_accents"]['orange'])

    ax.set_xlabel("residue distance", fontsize=fontsize)

    ax.set_xticks(range(0, df.index.max()))
    # data gets messy due to low numbers after 15 residues.
    ax.set_xlim(0, 15)
    # add a bit of height, so the legend does not overlap data
    ax.set_ylim(mean_ser.min(), mean_ser.max() + 0.1)
    figpath = coev_vs_res_dist_xlsx[:-5] + "_coev" + ".png"
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    ax.legend(handles, labels, ncol=1, loc=1, fontsize=fontsize, frameon=True, facecolor='white')
    #fig.legend(fontsize=fontsize, loc="upper right", bbox_to_anchor=[0.85, 0.95])
    fig.tight_layout()
    fig.savefig(figpath, dpi=240)
    fig.savefig(figpath[:-4] + ".pdf")

def calc_retrospective_coev_from_list_interf_res(s, df_set, logging):

    logging.info('calc_retrospective_coev_from_list_interf_res starting')
    retrospective_coev_xlsx = os.path.join(s["set_results_folder"], "{}_retrospective_coev.xlsx".format(s["setname"]))
    writer = pd.ExcelWriter(retrospective_coev_xlsx)

    randomise_int_res = False
    remove_residues_outside_interface_region = False
    logging.info("randomise_int_res = {}, remove_residues_outside_interface_region = {}".format(randomise_int_res, remove_residues_outside_interface_region))
    InterResList_of_last_TMD = None
    NoninterResList_of_last_TMD = None
    TMD_start_of_last_TMD = None

    for XI in ["MI", "DI"]:
        sub_dict = {}

        for i in df_set.index:
            sys.stdout.write(".")
            sys.stdout.flush()
            if i == 0:
                is_first_TMD = True
            else:
                is_first_TMD = False
            sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD = calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set,
                                                                                                                                                              randomise_int_res, InterResList_of_last_TMD,
                                                                                                                                                              NoninterResList_of_last_TMD, TMD_start_of_last_TMD,
                                                                                                                                                              remove_residues_outside_interface_region)

        if randomise_int_res == True:
            # need to add the data for the first TMD, which was skipped above
            i = 0
            is_first_TMD = False
            sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD = calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set,
                                                                                                                                                              randomise_int_res, InterResList_of_last_TMD,
                                                                                                                                                              NoninterResList_of_last_TMD, TMD_start_of_last_TMD,
                                                                                                                                                              remove_residues_outside_interface_region)
        # save each MI and DI separately
        df_retro = pd.DataFrame(sub_dict).T
        df_retro.to_excel(writer, sheet_name = XI)

        create_quick_plot = True
        if create_quick_plot:
            retrospective_coev_plot = retrospective_coev_xlsx[:-5] + XI + ".png"
            df_retro["inter_larger"] = df_retro.AverageInter > df_retro.AverageNoninter
            fig, ax = plt.subplots()
            df_retro[["AverageInter", "AverageNoninter"]].plot(kind="bar", ax=ax)
            fig.tight_layout()
            fig.savefig(retrospective_coev_plot)

        # drop any (rare) rows without data, where the interface region was outside the length of the TMD?
        df_retro.dropna(how="any", inplace=True)

        vc = df_retro["inter_larger"].value_counts()
        if True in vc.index.tolist():
            n_TMDs_with_higher_int = vc[True]
        else:
            n_TMDs_with_higher_int = 0

        perc_higher_int = n_TMDs_with_higher_int / df_retro.shape[0]
        logging.info("\n{:.2f} % ({}/{}) of TMDs have higher {} of interface than non-interface".format(perc_higher_int*100, n_TMDs_with_higher_int, df_retro.shape[0], XI))

        logging.info("\n\nmean values\n{}\n".format(df_retro.mean()))

        t_value, p_value = ttest_ind(df_retro.AverageInter, df_retro.AverageNoninter)

        #logging.info("remove_residues_outside_interface_region = {}".format(remove_residues_outside_interface_region))
        logging.info("\np-value for average {} coevolution of interface vs non-interface = {:.03f}".format(XI, p_value))

    writer.close()
    sys.stdout.write("\n")
    logging.info('calc_retrospective_coev_from_list_interf_res finished')


def calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set, randomise_int_res, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region):
    """Calculate average fraction of DI for a single protein.

    PANDAS METHOD USED FOR ETRA DATASET.

    SPLIT INTO A SUB-FUNCTION FOR USE DURING RANDOMISATION.

    Parameters
    ----------
    sub_dict : dict
        dictionary for each TMD. Separate dicts are made for MI and DI.
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    i : int
        TMD number
    XI : str
        "MI" or "DI"
    is_first_TMD : bool
        whether the TMD is the first one, without initial randomised data
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    randomise_int_res : bool
        whether the interface residues shold be randomised
    InterResList_of_last_TMD : list
        List of interface residues from the last TMD. To be used during randomisation.
    NoninterResList_of_last_TMD : list
        List of non-interface residues from the last TMD. To be used during randomisation.
    TMD_start_of_last_TMD : int
        TMD start used to convert residue positions to range index

    Returns
    -------
    sub_dict : dict
        dictionary for each TMD. Separate dicts are made for MI and DI.
    InterResList_of_last_TMD : list
        List of interface residues from the last TMD. To be used during randomisation.
    NoninterResList_of_last_TMD : list
        List of non-interface residues from the last TMD. To be used during randomisation.
    TMD_start_of_last_TMD : int
        TMD start used to convert residue positions to range index
    """
    acc = df_set.loc[i, "acc"]
    database = df_set.loc[i, "database"]
    TMD_start = int(df_set.loc[i, "TMD_start"])
    TMD_end = int(df_set.loc[i, "TMD_end"])
    TMD_len = TMD_end - TMD_start
    freecontact_file = os.path.join(s["feature_cumulative_coevolution"], database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
    feature_combined_file = os.path.join(s["features_folder"], "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

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

    dfp = df.pivot_table(index="n1", columns="n2", values=XI)

    """ asymmetrical pivoted data

         235       236       237       238       239       240 ...        252       253       254       255       256       
    n1                                                     ...                                                              
    235   0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523 
    236        NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577 
    237        NaN       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414 
    238        NaN       NaN       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822 
    """
    # get full list of residues
    position_list = range(TMD_start, TMD_end + 1)
    dfp = dfp.reindex(index=position_list, columns=position_list)
    # put data on both sides of the table for easy indexing
    for col in dfp.columns:
        start = col + 1
        dfp.loc[start:, col] = dfp.loc[col, start:]


    """dfp now contains the coevolution data as a symmetrical dataframe, from each residue to the other

              192       193       194       195       196       197       198       199       200       201    ...          210       211       212       213       214       215       216       217       218       219
    n1                                                                                                         ...                                                                                                       
    192       NaN  0.092237  0.026186  0.126701  0.108622  0.107383  0.075048  0.070287  0.084822  0.037957    ...     0.152848  0.074908  0.073767  0.159693  0.044335  0.092576  0.057039  0.176549  0.076715  0.066157
    193  0.092237       NaN  0.089528  0.137392  0.112203  0.153103  0.114659  0.173971  0.134006  0.091982    ...     0.237441  0.107704  0.097004  0.216488  0.146309  0.100271  0.101273  0.301949  0.105543  0.193257
    194  0.026186  0.089528       NaN  0.102470  0.078647  0.138274  0.141817  0.142261  0.133799  0.079009    ...     0.172375  0.111071  0.121039  0.171232  0.106160  0.095982  0.188747  0.230212  0.093526  0.217379
    195  0.126701  0.137392  0.102470       NaN  0.162021  0.124095  0.131162  0.248673  0.167416  0.094939    ...     0.179470  0.168239  0.139384  0.193543  0.102942  0.172607  0.153524  0.289339  0.113594  0.181711
    196  0.108622  0.112203  0.078647  0.162021       NaN  0.147395  0.106920  0.186598  0.170876  0.074893    ...     0.152920  0.130958  0.104620  0.165248  0.071461  0.117822  0.113831  0.243438  0.097208  0.153550
    197  0.107383  0.153103  0.138274  0.124095  0.147395       NaN  0.185372  0.300418  0.254464  0.135116    ...     0.294558  0.214323  0.237466  0.396039  0.111643  0.203568  0.221890  0.442481  0.167183  0.255704
    198  0.075048  0.114659  0.141817  0.131162  0.106920  0.185372       NaN  0.188028  0.174667  0.158833    ...     0.145839  0.134066  0.147938  0.256873  0.098789  0.146614  0.202526  0.266566  0.114003  0.211277
    199  
    """

    # open combined file with interface definitions
    dfc = pd.read_csv(feature_combined_file, index_col=0)

    # set the residue numbering as the index
    dfc.index = dfc.res_num_full_seq.astype(int)
    # get list of interface and noninterface residue positions
    InterResList = dfc.loc[dfc.interface == 1].index
    NoninterResList = list(dfc.loc[dfc.interface == 0].index)

    # get dataframe of distances between the residues


    if remove_residues_outside_interface_region:
        #logging.info("orig NoninterResList = {}".format(NoninterResList))
        lowest_interface_res = InterResList.min()
        highest_interface_res = InterResList.max()
        #logging.info("interface residue range = {} to {}".format(lowest_interface_res, highest_interface_res))
        NoninterResList = [x for x in NoninterResList if lowest_interface_res < x < highest_interface_res]
        #logging.info("final NoninterResList = {}".format(NoninterResList))

    if randomise_int_res == False:
        # calculate mean coevolution values for the desired selection of the dataframe.
        # Note that the values are symmetric and doubled ([235,236] and also [236,235])
        # but the mean will be unaffected
        mean_XI_interface = dfp.loc[InterResList, InterResList].mean().mean()
        mean_XI_noninterface = dfp.loc[NoninterResList, NoninterResList].mean().mean()

        sub_dict[acc] = {"AverageInter": mean_XI_interface, "AverageNoninter": mean_XI_noninterface}
    elif randomise_int_res == True and is_first_TMD == True:
        # can't used the interface residues from previous protein
        pass
    elif randomise_int_res == True and is_first_TMD != True:
        # FOR RANDOMISATION, THE ORIGINAL INDEXING BASED ON AA NUMBERS MUST BE REPLACED BY RANGE INDEXING

        # convert from amino acid numbering to a range index by subtracting the TMD_start
        InterResList_of_last_TMD = pd.Series(InterResList_of_last_TMD) - TMD_start_of_last_TMD
        NoninterResList_of_last_TMD = pd.Series(NoninterResList_of_last_TMD) - TMD_start_of_last_TMD
        # drop any residue positions that are longer than the TMD (i.e., where index TMD is longer than TMD supplying coevolution data)
        InterResList_of_last_TMD = InterResList_of_last_TMD[InterResList_of_last_TMD <= TMD_len].tolist()
        NoninterResList_of_last_TMD = NoninterResList_of_last_TMD[NoninterResList_of_last_TMD <= TMD_len].tolist()

        # convert pivoted data to range index
        dfp.index = range(dfp.shape[0])
        dfp.columns = range(dfp.shape[1])

        # slice out the interface and non-interface residues, as above
        mean_XI_interface = dfp.loc[InterResList_of_last_TMD, InterResList_of_last_TMD].mean().mean()
        mean_XI_noninterface = dfp.loc[NoninterResList_of_last_TMD, NoninterResList_of_last_TMD].mean().mean()

        sub_dict[acc] = {"AverageInter": mean_XI_interface, "AverageNoninter": mean_XI_noninterface}

    InterResList_of_last_TMD = InterResList
    NoninterResList_of_last_TMD = NoninterResList
    TMD_start_of_last_TMD = TMD_start

    return sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD

def get_array_dist_separating_res_in_list(pos_list):
    """Get array of the distance separating the residues in a list.

    Parameters
    ----------
    pos_list : list
        list of positions that are interacting, e.g [28,29,32,33,34,36]

    Returns
    -------
    dist_2D_arr np.ndarray
        numpy array showing the distances between all combinations in the list of residues.
        E.g.
        array([[ nan,   1.,   4.,   5.,   6.,   8.],
               [  1.,  nan,   3.,   4.,   5.,   7.],
               [  4.,   3.,  nan,   1.,   2.,   4.],
               [  5.,   4.,   1.,  nan,   1.,   3.],
               [  6.,   5.,   2.,   1.,  nan,   2.],
               [  8.,   7.,   4.,   3.,   2.,  nan]])

        the dataframe equivalent looks like this:
             28   29   32   33   34   36
        28  NaN  1.0  4.0  5.0  6.0  8.0
        29  1.0  NaN  3.0  4.0  5.0  7.0
        32  4.0  3.0  NaN  1.0  2.0  4.0
        33  5.0  4.0  1.0  NaN  1.0  3.0
        34  6.0  5.0  2.0  1.0  NaN  2.0
        36  8.0  7.0  4.0  3.0  2.0  NaN

    """
    # make empty array
    le = len(pos_list)
    dist_2D_arr = np.zeros(le*le).reshape(le, le)

    # calculate distances for each position separately. Add to array
    for i, pos in enumerate(pos_list):
        arr = np.array(pos_list)
        arr = abs(arr - pos)
        dist_2D_arr[i, :] = arr
    # replace 0 with nan, so that it doesn't count in the averages
    dist_2D_arr[dist_2D_arr == 0] = np.nan
    return dist_2D_arr

def calc_retrospective_coev_from_struct_contacts(s, dfset, logging):
    """Calculate retrospective coevolution from structural data with interpair contacts

    Calculates coevolution of interface (contacting) residues from homodimer structures
    using the method of Wang and Barth (2015).

    Takes csv with list of interacting residues as an input:
    acc	    inter1	inter2
    1orqC4	14	    20
    1orqC4	14	    23
    1orqC4	14	    24
    1orqC4	18	    24
    1xioA4	1	    2
    1xioA4	5	    6

    Converts this to a list of interacting pairs:
     - saved as InterPairList by the calc_retrospective_coev_from_struct_contacts_single_prot function
     - this is named as "last TMD" because the randomisation takes took? the positions of first TMD and applied them to the last

    Takes pairwise coevolution values from FreeContact output file.
    E.g. D:\Dropbox\tm_homodimer_dropbox\THOIPA_data\Features\cumulative_coevolution\ NMR\O15455.surr20.gaps5.freecontact.csv
    Data looks like this:
    1 F 2 F 0.195863 0.552187
    1 F 3 M 0.172853 -0.530669
    1 F 4 I 0.406909 0.445122
    1 F 5 N 0.245805 2.64414

    Finds coevolution scores for all interacting residue pairs.
     - calculates average

    Finds coevolution scores for all non-interacting pairs, separated by no longer than 8 residues.
     - creates list of non-interacting pairs (NonInterPairList)

    REPEATS THE ABOVE WITH RANDOM POSITIONS, TAKEN FROM A DIFFERENT TMD IN THE LIST OF PROTEINS/TMDS TO TEST



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
    logging.info('calc_retrospective_coev_from_list_interf_res starting')
    #retrospective_coev_xlsx = os.path.join(s["set_results_folder"], "set04_retrospective_coev.xlsx")
    randomise_int_res = True
    if randomise_int_res == True:
        retrospective_coev_xlsx = os.path.join(s["set_results_folder"], "set04_randomise_retrospective_coev.xlsx")
        writer = pd.ExcelWriter(retrospective_coev_xlsx)
    else:
        retrospective_coev_xlsx = os.path.join(s["set_results_folder"], "set04_retrospective_coev.xlsx")
        writer = pd.ExcelWriter(retrospective_coev_xlsx)
    remove_residues_outside_interface_region = False
    logging.info("randomise_int_res = {}, remove_residues_outside_interface_region = {}".format(randomise_int_res, remove_residues_outside_interface_region))
    InterPairList_of_last_TMD = None
    NonInterPairList_of_last_TMD = None
    TMD_start_of_last_TMD = None
    crystal_NMR_interpair_file = os.path.join(s["set_results_folder"], "Average_Fraction_DI","Crystal_NMR_interpair.csv")
    pd_int = pd.read_csv(crystal_NMR_interpair_file, engine="python")
    """pd_int looks like this    
    acc	    inter1	inter2
    1orqC4	14	    20
    1orqC4	14	    23
    1orqC4	14	    24
    1orqC4	18	    24
    1xioA4	1	    2
    1xioA4	5	    6
    """

    # get list of proteins from set04, sort by TMD length
    # MT NOTE: Set04 (NMR and crystal) is currently hard-coded here
    set04_file = os.path.join(s["sets_folder"],"set04_crystal_NMR.xlsx")
    df_set = pd.read_excel(set04_file,sheetname="proteins")

    print(df_set['TMD_seq'].str.len())
    print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz")

    df_set.index = df_set['TMD_seq'].str.len()
    df_set = df_set.sort_index(ascending=True).reset_index(drop=True)

    print(df_set['TMD_seq'].str.len())
    print("yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy")
    sub_dict = {}
    NoninterPairList_dict = {}
    InterPairList_dict = {}
    for i in df_set.index:
        acc = df_set.at[i, "acc"]
        TMD_len = len(df_set.at[i, "TMD_seq"])
        sys.stdout.write(".")
        sys.stdout.flush()
        if i == 0:
            is_first_TMD = True
        else:
            is_first_TMD = False
            InterPairList_of_last_TMD_max = np.array(InterPairList_of_last_TMD).max()
            NonInterPairList_of_last_TMD_max = np.array(NonInterPairList_of_last_TMD).max()
            print("InterPairList_of_last_TMD_max", InterPairList_of_last_TMD_max, "NonInterPairList_of_last_TMD_max", NonInterPairList_of_last_TMD_max, "TMD_len", TMD_len)
            TMD_len = len(df_set.at[i, "TMD_seq"])
            print("TMD_len", TMD_len)
        #if not all([randomise_int_res == True and i == 0]):

        #is_first_TMD = False
        sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD = calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, df_set, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region)
        # save all lists of interacting and non-interacting pairs
        InterPairList_dict[acc] = str(InterPairList)
        NoninterPairList_dict[acc] = str(NoninterPairList)

        InterPairList_of_last_TMD = InterPairList
        NonInterPairList_of_last_TMD = NoninterPairList


    if randomise_int_res == True:
        # need to add the data for the first TMD, which was skipped above
        i = 0
        acc = df_set.at[i, "acc"]
        is_first_TMD = False
        TMD_len = len(df_set.at[i, "TMD_seq"])
        print("TMD_len", TMD_len)
        print("InterPairList_of_last_TMD", InterPairList_of_last_TMD)
        print("starting first TMD")
        sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD = calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, df_set, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region)
        print("real InterPairList of the first TMD", InterPairList)



    df_retro = pd.DataFrame(sub_dict).T
    df_retro.to_excel(writer, sheet_name="DI")

    """Save lists of interacting or non-interacting residue pairs:
    E.g. 
            InterPair	                                 NonInterPair
    1orqC4	[[14, 20], [14, 23], [14, 24], [18, 24]]	[[1, 2], [2, 1], [1, 3], [3, 1], [1, 4], [4, 1], [1, 5], [5, 1], ..........
    1xioA4	[[1, 2], [5, 6], [12, 12], [16, 16]]	    [[1, 3], [3, 1], [1, 4], [4, 1], [1, 5], [5, 1], [1, 6], [6, 1], .............
    """
    InterPair_ser = pd.Series(InterPairList_dict)
    NonInterPair_ser = pd.Series(NoninterPairList_dict)
    df_res_list = pd.DataFrame([InterPair_ser, NonInterPair_ser], index=["InterPair", "NonInterPair"]).T
    df_res_list.to_excel(writer, sheet_name="PairLists")

    create_quick_plot = True
    if create_quick_plot:
        retrospective_coev_plot = retrospective_coev_xlsx[:-5]  + "DI.png"
        df_retro["inter_larger"] = df_retro.AverageInter > df_retro.AverageNoninter
        fig, ax = plt.subplots()
        df_retro[["AverageInter", "AverageNoninter"]].plot(kind="bar", ax=ax)
        fig.tight_layout()
        fig.savefig(retrospective_coev_plot)

    # drop any (rare) rows without data, where the interface region was outside the length of the TMD?
    df_retro.dropna(how="any", inplace=True)

    vc = df_retro["inter_larger"].value_counts()
    if True in vc.index.tolist():
        n_TMDs_with_higher_int = vc[True]
    else:
        n_TMDs_with_higher_int = 0

    perc_higher_int = n_TMDs_with_higher_int / df_retro.shape[0]
    logging.info(
        "\n{:.2f} % ({}/{}) of TMDs have higher DI of interface than non-interface".format(perc_higher_int * 100,
                                                                                           n_TMDs_with_higher_int,
                                                                                           df_retro.shape[0]))

    logging.info("\n\nmean values\n{}\n".format(df_retro.mean()))

    t_value, p_value = ttest_ind(df_retro.AverageInter, df_retro.AverageNoninter)

    # logging.info("remove_residues_outside_interface_region = {}".format(remove_residues_outside_interface_region))
    logging.info("\np-value for average DI coevolution of interface vs non-interface = {:.03f}".format( p_value))

    writer.close()
    sys.stdout.write("\n")
    logging.info('calc_retrospective_coev_from_list_interf_res finished')

def calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, df_set, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region):
    """Calculate retrospective coevolution from structural data with interpair contacts, for a single protein.

    see docstring above for calc_retrospective_coev_from_struct_contacts

    Parameters
    ----------
    sub_dict
    s
    pd_int
    i
    logging
    is_first_TMD
    df_set
    randomise_int_res
    InterPairList
    NonInterPairList
    TMD_start_of_last_TMD
    remove_residues_outside_interface_region

    Returns
    -------

    """
    acc = df_set.loc[i,'acc']
    database = df_set.loc[i,'database']
    InterPairList= pd_int[["inter1", "inter2"]][pd_int["acc"] == acc].values.tolist()
    interlist = []
    for x in InterPairList:
        interlist.extend(x)
    lowest_interface_res = min(interlist)
    highest_interface_res = max(interlist)
    NoninterPairList = []
    freecontact_file = os.path.join(s["features_folder"], "cumulative_coevolution", database, "{}.surr20.gaps5.freecontact.csv".format(acc))
    inter_within8_dict = {}

    if os.path.isfile(freecontact_file):
        with open(freecontact_file, 'r') as f:
            for line in f:
                arr = line.strip().split()
                tmd_len = int(arr[2])
                if int(arr[2]) - int(arr[0]) <= 8:
                    inter_within8_dict[arr[0] + '_' + arr[2]] = arr[5]
                    inter_within8_dict[arr[2] + '_' + arr[0]] = arr[5]
                    if [int(arr[0]), int(arr[2])] not in InterPairList and [int(arr[0]), int(
                        arr[2])] not in InterPairList:
                        if remove_residues_outside_interface_region:
                            if lowest_interface_res < int(arr[0]) < highest_interface_res and lowest_interface_res < int(arr[2]) < highest_interface_res:
                                NoninterPairList.append([int(arr[0]), int(arr[2])])
                                NoninterPairList.append([int(arr[2]), int(arr[0])])
                        else:
                            NoninterPairList.append([int(arr[0]), int(arr[2])])
                            NoninterPairList.append([int(arr[2]), int(arr[0])])
        f.close()

    if randomise_int_res == False:
        inter_DI = []
        ninter_DI = []
        for key, value in inter_within8_dict.items():
            inter_pair = [int(x) for x in key.split('_')]
            if inter_pair in InterPairList:
                inter_DI.append(float(inter_within8_dict[key]))
            elif inter_pair in NoninterPairList:
                ninter_DI.append(float(inter_within8_dict[key]))
        average_inter = np.mean(inter_DI)
        average_ninter = np.mean(ninter_DI)
        #print(acc,average_inter,average_ninter, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)
        sub_dict[acc] = {"AverageInter": average_inter, "AverageNoninter": average_ninter}

    # Randomisation of first TMD not possible. Skip and return the true interface only.
    elif randomise_int_res == True and is_first_TMD == True:
        #sub_dict = "First TMD when randomise_int_res == True"
        #print("first TMD here")
        #return sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD
        pass

    elif randomise_int_res == True and is_first_TMD != True:
        # random interpair lists taken from a previous TMD
        InterPairList_rand =[x for x in InterPairList_of_last_TMD if x[0] <= tmd_len and x[1] <= tmd_len ]
        NonInterPairList_rand = [x for x in NonInterPairList_of_last_TMD if x[0] <= tmd_len and x[1] <= tmd_len ]

        orig_max = np.array(InterPairList_of_last_TMD).max()
        final_max_after_filtering_out_res_too_long = np.array(InterPairList_rand).max()
        print(orig_max == final_max_after_filtering_out_res_too_long)
        if orig_max != final_max_after_filtering_out_res_too_long:
            print("must be the last TMD, if sorting is done correctly")

        #print(acc,InterPairList_rand,NonInterPairList_rand, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)
        inter_DI = []
        ninter_DI = []
        for key, value in inter_within8_dict.items():
            inter_pair = [int(x) for x in key.split('_')]
            if inter_pair in InterPairList_rand:
                inter_DI.append(float(inter_within8_dict[key]))
            elif inter_pair in NonInterPairList_rand:
                ninter_DI.append(float(inter_within8_dict[key]))
        average_inter = np.mean(inter_DI)
        average_ninter = np.mean(ninter_DI)
        #print(acc,average_inter,average_ninter, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)
        sub_dict[acc] = {"AverageInter": average_inter, "AverageNoninter": average_ninter}

    #InterPairList = InterPairList
    #NonInterPairList = NoninterPairList


    return sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD




def create_inter_rr_number_bt1_set228_file(df_new, set_file):

    writer = pd.ExcelWriter(set_file)
    acc_list = []
    notes_list =[]
    database_list = []
    TMD_seq_list =[]
    full_seq_list =[]
    for i in range(df_new.shape[0]):
        if df_new.iloc[i]["nonself_intpair_num"].astype(int) >=1:
            acc = df_new.iloc[i]["pdb_id"] + df_new.iloc[i]["tm_numA"]
            note = acc
            database = "crystal"
            TMD_seq = df_new.iloc[i]["aligned_AB"]
            full_seq = df_new.iloc[i]["full_seqA"]
            acc_list.append(acc)
            notes_list.append(note)
            database_list.append(database)
            TMD_seq_list.append(TMD_seq)
            full_seq_list.append(full_seq)
    col_list = [acc_list, notes_list,database_list, TMD_seq_list, full_seq_list]
    df = pd.DataFrame.from_records(list(map(list, zip(*col_list))), columns = ["acc", "notes","database", "TMD_seq", "full_seq"])
    df.to_excel(writer, sheet_name = "proteins")
    writer.save()





def get_inter_pair_rr_num(df_new):
    non_self_inter_pair_num_list = []
    all_inter_pair_num_list = []
    for i in range(df_new.shape[0]):
        inter_pair_arr = df_new.iloc[i]["interpair_4.5"].split('+')
        all_inter_pair_num = len(inter_pair_arr) - 1
        non_self_inter_pair_num = all_inter_pair_num
        all_inter_pair_num_list.append(all_inter_pair_num)
        for j in range(1,len(inter_pair_arr)):
            A_seq_res = inter_pair_arr[j].split('_')[1]
            B_seq_res = inter_pair_arr[j].split('_')[5]
            if A_seq_res == B_seq_res:
                non_self_inter_pair_num = non_self_inter_pair_num -1
        non_self_inter_pair_num_list.append(non_self_inter_pair_num)
    return non_self_inter_pair_num_list, all_inter_pair_num_list




def get_cdhit60_represent_interpair_df(df_rd, cdhit60_crystal035_fasta_file):
    inter_pair_dict = {}
    with open(cdhit60_crystal035_fasta_file, 'r') as f:
        for line in f:
            if re.search(r'^>',line.strip()):
                inter_pair= line.strip()[1:]
                inter_pair_dict[inter_pair] = 1
    drop_list = []
    for i in range(df_rd.shape[0]):
        pair = df_rd.iloc[i]["pdb_id"] + "_" + df_rd.iloc[i]["tm_numA"] + "_" + df_rd.iloc[i]["tm_numB"]
        if pair not in inter_pair_dict:
            drop_list.append(i)
    df_rd.drop(df_rd.index[drop_list], inplace=True)
    df_rd.reset_index(drop=True,inplace=True)
    return df_rd



def get_interpair_cdhist_cluster_number(df_rd,cdhit60_crystal035_clr_file):
    dict_cluster_pair = {}
    with open(cdhit60_crystal035_clr_file, 'r') as f:
        for line in f:
            m = re.search(r'^>Cluster\s*(\d+).*', line)
            if m:
                cluster_number = m.group(1)
            else:
                pair = re.search(r'\d+\s*\d+aa,\s*>(.*)\.\.\.\s*.*',line)
                if pair:
                    dict_cluster_pair[pair.group(1)] = cluster_number
    return dict_cluster_pair

    #print(cdhit60_crystal035_clr_file)
    #return df_rd

def get_pdb_experimental_resolution(s, pdb_id):
    trpdb_file = os.path.join(s["pdbtm_homodimer_folder"], "xml","{}.trpdb.gz".format(pdb_id))
    if os.path.isfile(trpdb_file):
        with gzip.open(trpdb_file, "rt") as f:
            for line in f:
                m = re.search(r'^REMARK\s*2\s*RESOLUTION.\s*([0-9]\.[0-9][0-9])\s*.*', line)
                if m:
                    return m.group(1)
    else:
        return None




def multiple_add_closedist_rr_inter_pairs(s,pdbtm_trpdb_file,tm_homo_pair_file,inter_def_cutoff):
    df_homo = pd.read_csv(tm_homo_pair_file, engine = "python" )
    index_len =len(df_homo.index)
    aligned_AB_list = []
    seq_minus_pdb_A_list = []
    seq_minus_pdb_B_list = []
    closedist_A_list = []
    closedist_B_list = []
    inter_pair_cutoff_list = []
    for i in range(index_len):
        pdb_id = df_homo.iloc[i]["pdb_id"]
        chainA = df_homo.iloc[i]["chainA"]
        chainB = df_homo.iloc[i]["chainB"]
        tm_numA = df_homo.iloc[i]["tm_numA"]
        tm_numB = df_homo.iloc[i]["tm_numB"]
        seq_begA = df_homo.iloc[i]["seq_begA"]
        seq_endA = df_homo.iloc[i]["seq_endA"]
        seq_begB = df_homo.iloc[i]["seq_begB"]
        seq_endB = df_homo.iloc[i]["seq_endB"]
        pdb_begA = df_homo.iloc[i]["pdb_begA"]
        pdb_endA = df_homo.iloc[i]["pdb_endA"]
        pdb_begB = df_homo.iloc[i]["pdb_begB"]
        pdb_endB = df_homo.iloc[i]["pdb_endB"]
        full_seqA = df_homo.iloc[i]["full_seqA"]
        full_seqB = df_homo.iloc[i]["full_seqB"]
        aligned_A = df_homo.iloc[i]["aligned_A"]
        aligned_B = df_homo.iloc[i]["aligned_B"]
        aligned_AB = utils.join_two_algined_seqences(aligned_A, aligned_B)
        aligned_AB_len = len(aligned_AB)
        seq_minus_pdb_A = seq_begA - pdb_begA
        seq_minus_pdb_B = seq_begB - pdb_begB

        aligned_AB_list.append(aligned_AB)
        seq_minus_pdb_A_list.append(seq_minus_pdb_A)
        seq_minus_pdb_B_list.append(seq_minus_pdb_B)
        try:
            aligned_AB_seq_startA = full_seqA.index(aligned_AB) + 1
            aligned_AB_seq_startB = full_seqB.index(aligned_AB) + 1
            closedist_A, closedist_B, inter_pair_cutoff = single_add_closedist_rr_inter_pairs( pdbtm_trpdb_file, chainA,
                                                                                     chainB, aligned_AB_len,
                                                                                     aligned_AB_seq_startA,
                                                                                     aligned_AB_seq_startB,
                                                                                     seq_minus_pdb_A, seq_minus_pdb_B,
                                                                                     inter_def_cutoff)
            closedist_A_list.append(closedist_A)
            closedist_B_list.append(closedist_B)
            inter_pair_cutoff_list.append(inter_pair_cutoff)
        except:
            sys.stdout.write("tm_seq_A and tm_seq_B has variant residues, skip this tm\n")
            closedist_A, closedist_B, inter_pair_cutoff = np.nan, np.nan,np.nan
            closedist_A_list.append(closedist_A)
            closedist_B_list.append(closedist_B)
            inter_pair_cutoff_list.append(inter_pair_cutoff)

    df_homo["aligned_AB"] = aligned_AB_list
    df_homo["seq_minus_pdb_A"] = seq_minus_pdb_A_list
    df_homo["seq_minus_pdb_B"] = seq_minus_pdb_B_list
    df_homo["closedist_A"] = closedist_A_list
    df_homo["closedist_B"] = closedist_B_list
    inter_pair_cutoff = "interpari_{}".format(inter_def_cutoff)
    df_homo[inter_pair_cutoff] = inter_pair_cutoff_list
    ##updatae *.interhelix.csv by adding closedist and rr interact pair clolumns
    tm_homo_pair_file_handle = open(tm_homo_pair_file,'w')
    df_homo.to_csv(tm_homo_pair_file_handle)
    tm_homo_pair_file_handle.close()


def single_add_closedist_rr_inter_pairs(pdbtm_trpdb_file,chainA, chainB, aligned_AB_len, aligned_AB_seq_startA, aligned_AB_seq_startB, seq_minus_pdb_A, seq_minus_pdb_B,inter_def_cutoff):
    '''
    calculate closedist for each chain, and get the interact rrs from two chains under the interact defination cutoff
    Parameters
    ----------
    pdbtm_trpdb_file
    chainA
    chainB
    aligned_AB_len
    aligned_AB_seq_startA
    aligned_AB_seq_startB
    seq_minus_pdb_A
    seq_minus_pdb_B
    inter_def_cutoff

    Returns
    -------
    closedist, interact_rr_pairs as str
    '''
    chain1_pdb_tm_start = aligned_AB_seq_startA - seq_minus_pdb_A
    chain1_pdb_tm_end = chain1_pdb_tm_start + aligned_AB_len - 1
    chain2_pdb_tm_start = aligned_AB_seq_startB - seq_minus_pdb_B
    chain2_pdb_tm_end = chain2_pdb_tm_start + aligned_AB_len - 1
    hash1arrayx = {}
    hash1arrayy = {}
    hash1arrayz = {}
    hash2arrayx = {}
    hash2arrayy = {}
    hash2arrayz = {}
    hashclosedist1 = {}
    hashclosedist2 = {}
    hash_inter_pair = {}
    with gzip.open(pdbtm_trpdb_file, "rt") as f:
        for line in f:
            if re.search("^ATOM", line):
                atom = line[12:16]
                if not re.search("^H", atom):  # non-H atom distance
                    index = line[6:11]
                    x = line[30:38]
                    y = line[38:46]
                    z = line[46:54]
                    chain = line[21:22]
                    residue_num_pdb = int(line[22:26])
                    residue_name = line[17:20]
                    residue_name = utils.shorten(residue_name)
                    if chain == chainA  and residue_num_pdb >= chain1_pdb_tm_start and residue_num_pdb <= chain1_pdb_tm_end:
                        residue_num_seq = residue_num_pdb + seq_minus_pdb_A
                        kk = ':'.join([chain, str(residue_num_seq), str(residue_num_pdb), residue_name])
                        if kk not in hash1arrayx:
                            hash1arrayx[kk] = x
                            hash1arrayy[kk] = y
                            hash1arrayz[kk] = z
                        else:
                            hash1arrayx[kk] = hash1arrayx[kk] + '_' + x
                            hash1arrayy[kk] = hash1arrayy[kk] + '_' + y
                            hash1arrayz[kk] = hash1arrayz[kk] + '_' + z
                    if chain == chainB and residue_num_pdb >= chain2_pdb_tm_start and residue_num_pdb <= chain2_pdb_tm_end:
                        residue_num_seq = residue_num_pdb + seq_minus_pdb_B
                        kk = ':'.join([chain,str(residue_num_seq), str(residue_num_pdb),residue_name])
                        if kk not in hash2arrayx:
                            hash2arrayx[kk] = x
                            hash2arrayy[kk] = y
                            hash2arrayz[kk] = z
                        else:
                            hash2arrayx[kk] = hash2arrayx[kk] + '_' + x
                            hash2arrayy[kk] = hash2arrayy[kk] + '_' + y
                            hash2arrayz[kk] = hash2arrayz[kk] + '_' + z

    for key1, count in hash1arrayx.items():
        arr_xvalue1 = hash1arrayx[key1].split('_')
        arr_yvalue1 = hash1arrayy[key1].split('_')
        arr_zvalue1 = hash1arrayz[key1].split('_')
        size1 = len(arr_xvalue1)
        for key2, count2 in hash2arrayx.items():
            for i in range(0, size1):
                arr_xvalue2 = hash2arrayx[key2].split('_')
                arr_yvalue2 = hash2arrayy[key2].split('_')
                arr_zvalue2 = hash2arrayz[key2].split('_')
                size2 = len(arr_xvalue2)
                for j in range(0, size2):
                    dist = math.sqrt((float(arr_xvalue1[i]) - float(arr_xvalue2[j])) ** 2 + (
                    float(arr_yvalue1[i]) - float(arr_yvalue2[j])) ** 2 + (
                                     float(arr_zvalue1[i]) - float(arr_zvalue2[j])) ** 2)
                    if dist <= inter_def_cutoff:
                        inter_pair = '_'.join([key1, key2])
                        if inter_pair not in hash_inter_pair:
                            hash_inter_pair[inter_pair] = dist
                        else:
                            if dist < hash_inter_pair[inter_pair]:
                                hash_inter_pair[inter_pair] = dist

                    if key1 not in hashclosedist1:
                        hashclosedist1[key1] = dist
                    else:
                        if dist < hashclosedist1[key1]:
                            hashclosedist1[key1] = dist

                    if key2 not in hashclosedist2:
                        hashclosedist2[key2] = dist
                    else:
                        if dist < hashclosedist2[key2]:
                            hashclosedist2[key2] = dist

    closedistA = utils.Get_Closedist_between_ChianA_ChainB(hashclosedist1)
    closedistB= utils.Get_Closedist_between_ChianA_ChainB(hashclosedist2)
    inter_pair_cutoff = utils.Get_Closedist_between_ChianA_ChainB(hash_inter_pair)
    return closedistA, closedistB ,inter_pair_cutoff



def get_inter_pair_homodimer(pdb_id, pdbtm_xml_path, logging):
    xml_file = os.path.join(pdbtm_xml_path, "{}.xml".format(pdb_id))
    tm_alpha_helix_file = os.path.join(pdbtm_xml_path, "{}.alpha_helix.csv".format(pdb_id))
    tm_inter_helix_file = os.path.join(pdbtm_xml_path, "{}.interhelix.csv".format(pdb_id))
    no_homodimer_pdb_folder = os.path.join(pdbtm_xml_path, "NoHomodimerPdb")
    utils.make_sure_path_exists(no_homodimer_pdb_folder)
    if not os.path.isfile(tm_alpha_helix_file):
        sys.stdout.write("{}\n,{}\n,{}\n could have been moved to folder:\n{} because this pdb contains no homodimer\n".format(xml_file,tm_alpha_helix_file,tm_inter_helix_file,no_homodimer_pdb_folder))
        sys.stdout.flush()
        return None
    #if not os.path.isfile(tm_inter_helix_file ):
    #if pdb_id == "1ap9":
    inter_pair_num = 0
    tm_inter_helix_file_handle = open(tm_inter_helix_file, 'w', newline='')
    writer = csv.writer(tm_inter_helix_file_handle)
    writer.writerow(['pdb_id', 'chainA', 'tm_numA', "seq_begA", "seq_endA", "pdb_begA", "pdb_endA", "tm_seqA","full_seqA", 'chainB', 'tm_numB', "seq_begB", "seq_endB", "pdb_begB", "pdb_endB", "tm_seqB", "full_seqB"
                        , "aligned_A", "aligned_B","seq_identity","gap_seq_identity","gap_num"])
    df = pd.read_csv(tm_alpha_helix_file, engine="python",index_col = "tm_num")
    index_len = len(df.index)
    for i in range(index_len):
        SequenceA = df.iloc[i]["tm_seq"]
        SeqA_UX_count = sum([1 for i in range(len(SequenceA)) if (SequenceA[i] == 'U' or SequenceA[i] == 'X')])
        ChainA = df.iloc[i]["chain"]
        for j in range(i+1, index_len):
            SequenceB = df.iloc[j]["tm_seq"]
            SeqB_UX_count = sum([1 for i in range(len(SequenceB)) if (SequenceB[i] == 'U' or SequenceB[i] == 'X')])
            ChainB = df.iloc[j]["chain"]
            if not ChainA == ChainB:
                alns = pairwise2.align.globalxx(SequenceA, SequenceB, one_alignment_only=True)
                best_aln = alns[0]
                aligned_A, aligned_B, score, begin, end = best_aln
                seq_id, g_seq_id, gap_num = utils.calculate_identity(aligned_A, aligned_B)
                if gap_num <=4 and seq_id >= 80 and SeqA_UX_count<=4 and SeqB_UX_count <=4:
                    inter_pair_num = inter_pair_num + 1
                    writer.writerow([df.iloc[i]["pdb_id"], df.iloc[i]["chain"], df.index[i],df.iloc[i]["seq_beg"],df.iloc[i]["seq_end"],df.iloc[i]["pdb_beg"],df.iloc[i]["pdb_end"],df.iloc[i]["tm_seq"], df.iloc[i]["full_seq"]
                          , df.iloc[j]["chain"], df.index[j], df.iloc[j]["seq_beg"], df.iloc[j]["seq_end"],
                          df.iloc[j]["pdb_beg"], df.iloc[j]["pdb_end"], df.iloc[j]["tm_seq"], df.iloc[j]["full_seq"], aligned_A, aligned_B, seq_id, g_seq_id, gap_num])
    del writer
    tm_inter_helix_file_handle.close()
    if inter_pair_num == 0:
        try:
            shutil.move(xml_file, no_homodimer_pdb_folder)
            shutil.move(tm_alpha_helix_file , no_homodimer_pdb_folder)
            shutil.move(tm_inter_helix_file , no_homodimer_pdb_folder)
            #os.remove(tm_inter_helix_file)
        except:
            sys.stdout.write("moving {} to folder:{} occures error, check".format(tm_inter_helix_file,no_homodimer_pdb_folder))
            sys.stdout.flush()




def get_tmd_alpha_helix_infor(pdb_id, pdbtm_xml_path, logging):
    '''
    extract alpha helix information from xml file
    Parameters
    ----------
    pdb_id
    pdbtm_xml_path
    logging

    Returns
    -------

    '''

    pdb_xml_file = os.path.join(pdbtm_xml_path, "{}.xml".format(pdb_id))
    if not os.path.isfile(pdb_xml_file):
        sys.stdout.write("{} could have been moved to nohomodimer folder\n".format(pdb_xml_file))
        sys.stdout.flush()
        return None
    else:
        tm_alpha_helix_file = os.path.join(pdbtm_xml_path, "{}.alpha_helix.csv".format(pdb_id))
        if not os.path.isfile(tm_alpha_helix_file):
            tm_alpha_helix_file_handle = open(tm_alpha_helix_file, 'w',newline='')
            writer = csv.writer(tm_alpha_helix_file_handle)
            writer.writerow(['pdb_id', 'chain', 'tm_num', "seq_beg", "seq_end", "pdb_beg", "pdb_end", "tm_seq", "full_seq"])
            #if pdb_id == "5a3q":
            tree = ET.parse(pdb_xml_file)
            root = tree.getroot()
            for chain in root.findall('{http://pdbtm.enzim.hu}CHAIN'):
                ChainId = chain.attrib["CHAINID"]
                seq = chain.find('{http://pdbtm.enzim.hu}SEQ').text
                seq1 = ''.join(seq.split())
                i = 0
                for reg in chain.findall('{http://pdbtm.enzim.hu}REGION'):
                    seq_beg = int(reg.attrib["seq_beg"])
                    pdb_beg = int(reg.attrib["pdb_beg"])
                    seq_end = int(reg.attrib["seq_end"])
                    pdb_end = int(reg.attrib["pdb_end"])
                    type = reg.attrib["type"]
                    if type == "H":
                        i = i + 1
                        TMnum = ChainId + str(i)
                        TMseq = seq1[(seq_beg - 1): seq_end]
                        writer.writerow([pdb_id, ChainId, TMnum,str(seq_beg),str(seq_end),str(pdb_beg),str(pdb_end), TMseq, seq1])
            del writer
            tm_alpha_helix_file_handle.close()



            # for Chain in root.findall('CHAIN'):
            #     print(Chain)
            #     print(Chain.tag, Chain.attrib)

def get_xml_file(pdb_id, pdbtm_xml_path, logging):
    """
    download pdbtm xml file for single alpha helix tmps
    Parameters
    ----------
    pdb_id: str. pdb ID
    pdbtm_xml_path: tmp xml file output path
    logging

    Returns
    -------

    """
    utils.make_sure_path_exists(pdbtm_xml_path)
    pdb_xml_file = os.path.join(pdbtm_xml_path,  "{}.xml".format(pdb_id))
    if not os.path.isfile(pdb_xml_file):
        xml_file_url = 'http://pdbtm.enzim.hu/data/database/{}/{}.xml'.format(pdb_id[1:3],pdb_id)
        print("downloading xml file with urllib")
        urllib.request.urlretrieve(xml_file_url,pdb_xml_file )
        logging.info("Output file:     %s\n" % pdb_xml_file)
        if not os.path.exists(pdb_xml_file):
            logging.info('********************xml download failed for : %s***************' % pdb_xml_file)

def get_trpdb_file(pdb_id, pdbtm_trpdb_file, logging):
    '''
    download single trpdb file
    Parameters
    ----------
    pdb_id
    tr_pdb_path
    logging

    Returns
    -------

    '''
    trpdb_file_url = 'http://pdbtm.enzim.hu/data/database/{}/{}.trpdb.gz'.format(pdb_id[1:3],pdb_id)
    print(trpdb_file_url)
    try:
        urllib.request.urlretrieve(trpdb_file_url,pdbtm_trpdb_file )
        logging.info("Output file:     %s\n" % pdbtm_trpdb_file)
    except:
        logging.info('********trpdb download failed:{}\n'.format(pdb_id))
    if not os.path.exists(pdbtm_trpdb_file):
        logging.info('********************trpdb download failed for : %s***************' % pdbtm_trpdb_file)





