import csv
import glob
import gzip
import math
import os
import re
import shutil
import sys
import tarfile
import urllib.request
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from Bio import pairwise2

import thoipapy
import thoipapy.utils as utils


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
    pdbtm_alpha_list = os.path.join(thoipapy_module_path, "deprecated", "pdbtm_alpha_list","pdbtm_alpha_20170616.list")
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
        sys.stdout.write(pdb_id)
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
        df_homo_matrix = df_homo.to_numpy()
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
        sys.stdout.write("{} could not be deleted".format(redundant_interact_homodimer_file))
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
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar, os.path.dirname(redundant_interact_homodimer_targz_file))
    df_rd = pd.read_csv(redundant_interact_homodimer_file)
    resolve_list = []
    for i in range(df_rd.shape[0]):
        pdb_id = str(df_rd.iloc[i]["pdb_id"])
        resolv = get_pdb_experimental_resolution(s, pdb_id)
        resolve_list.append(resolv)
    df_rd["resolv"] = resolve_list
    df_rd =df_rd[df_rd["resolv"].astype(float) <= 3.5]

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
    sys.stdout.write(df_rd[["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]])
    a = df_rd["cdhit0.6_clusternum"].unique()
    # for j in sorted(a):
    #     #sys.stdout.write(df_rd[df_rd["cdhit0.6_clusternum"].astype(int) == j,["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]])
    #     sys.stdout.write(df_rd[["pdb_id", "cdhit0.6_clusternum", "total_intpair_num", "nonself_intpair_num"]].loc[df_rd['cdhit0.6_clusternum'] == j])

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
    set_file = os.path.join(s["base_dir"], "sets",'set228_crystal.xlsx')
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
        bind_closedist_file = os.path.join(s['data_dir'], "features", "structure", "crystal", "{}.{}pairmax.bind.closedist.csv".format(pdb_id_chain,inter_pair_max))
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
                sys.stdout.write("residue{}not in {} closedist".format(j, pdb_id_chain))
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

    #sys.stdout.write(cdhit60_crystal035_clr_file)
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
            if re.search(r"^ATOM", line):
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
            #     sys.stdout.write(Chain)
            #     sys.stdout.write(Chain.tag, Chain.attrib)

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
        logging.info("downloading xml file with urllib")
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
    sys.stdout.write(trpdb_file_url)
    try:
        urllib.request.urlretrieve(trpdb_file_url,pdbtm_trpdb_file )
        logging.info("Output file:     %s\n" % pdbtm_trpdb_file)
    except:
        logging.info('********trpdb download failed:{}\n'.format(pdb_id))
    if not os.path.exists(pdbtm_trpdb_file):
        logging.info('********************trpdb download failed for : %s***************' % pdbtm_trpdb_file)





