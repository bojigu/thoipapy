import re
import os
import sys
import pandas as pd
import numpy as np
import math
import thoipapy
import glob

cutoff = 4.0
def calc_closedist_from_NMR_best_model(s):
    NMR_data_folder = r"D:\THOIPA_data\NMR_data"
    nmr_data_infor_file = os.path.join(NMR_data_folder,"zprotein_pdb_chain_tm_pdbseq_uniprotseq_full")
    with open(nmr_data_infor_file,'r') as f:
        for line in f:
            arr = line.strip().split()
            chain1 = "A"
            chain2 = "B"
            if(arr[1] == chain1):
                chain1_pdb_tm_start = int(arr[3]) - (int(arr[7]) -int(arr[5]))
                chain1_pdb_tm_end = int(arr[4]) - (int(arr[7]) -int(arr[5]))
                if chain1_pdb_tm_start <= 0:
                    chain1_pdb_tm_end = chain1_pdb_tm_end +1
                # if pdb == "2JWA":
                #     chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) -141)
                #     chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) - 141)
                # if pdb == "2MK9":
                #     chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) -102)
                #     chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) - 102)
                # if pdb == "2MJO":
                #     chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) - 101)
                #     chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) - 101)
                # if pdb == "2N90":
                #     chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) - 102)
                #     chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) - 102)
                # if pdb == "2m20":
                #     chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) - 61)
                #     chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) - 61)
            if(arr[1] == chain2):
                hash1arrayx = {}
                hash1arrayy = {}
                hash1arrayz = {}
                hash2arrayx = {}
                hash2arrayy = {}
                hash2arrayz = {}
                protein = arr[2]
                pdb = arr[0]
                chain2_pdb_tm_start = int(arr[3]) - (int(arr[7]) -int(arr[5]))
                chain2_pdb_tm_end = int(arr[4]) - (int(arr[7]) -int(arr[5]))
                if chain2_pdb_tm_start <= 0:
                    chain2_pdb_tm_end = chain2_pdb_tm_end + 1
                num = 0
                pdb_file_NMR = os.path.join(NMR_data_folder, "{}.pdb".format(pdb.upper()))
                NMR_closedist_file = os.path.join(NMR_data_folder, "{}_{}_{}.closedist.csv".format(protein,pdb,cutoff))
                hashclosedist = {}
                print(pdb)
                #if pdb =="1afo":
                with open(pdb_file_NMR, 'r') as file:
                    for row in file:

                        # if re.search('^MODEL\s+2\s+', row):
                        #     break
                        if re.search("^MODEL\s+2",row):
                            break
                        if re.search("^ATOM", row):
                            atom = row[12:16]
                            if not re.search("^\s*H", atom):  # non-H atom distance
                                index = row[6:11]
                                x = row[30:38]
                                y = row[38:46]
                                z = row[46:54]
                                chain = row[21:22]
                                residue_num = int(row[22:26])
                                residue_name = row[17:20]
                                if chain == chain1 and chain1_pdb_tm_start <= residue_num and residue_num <= chain1_pdb_tm_end:
                                    kk = ':'.join([pdb, chain, str(residue_num - chain1_pdb_tm_start), residue_name])
                                    if kk not in hash1arrayx:
                                        hash1arrayx[kk] = x
                                        hash1arrayy[kk] = y
                                        hash1arrayz[kk] = z
                                    else:
                                        hash1arrayx[kk] = hash1arrayx[kk] + ':' + x
                                        hash1arrayy[kk] = hash1arrayy[kk] + ':' + y
                                        hash1arrayz[kk] = hash1arrayz[kk] + ':' + z

                                if chain == chain2 and chain2_pdb_tm_start <= residue_num and residue_num <= chain2_pdb_tm_end:
                                    kk = ':'.join([pdb, chain, str(residue_num - chain2_pdb_tm_start), residue_name])
                                    if kk not in hash2arrayx:
                                        hash2arrayx[kk] = x
                                        hash2arrayy[kk] = y
                                        hash2arrayz[kk] = z
                                    else:
                                        hash2arrayx[kk] = hash2arrayx[kk] + ':' + x
                                        hash2arrayy[kk] = hash2arrayy[kk] + ':' + y
                                        hash2arrayz[kk] = hash2arrayz[kk] + ':' + z

                for key1, count in hash1arrayx.items():
                    arr_xvalue1 = hash1arrayx[key1].split(':')
                    arr_yvalue1 = hash1arrayy[key1].split(':')
                    arr_zvalue1 = hash1arrayz[key1].split(':')
                    size1 = len(arr_xvalue1)
                    for key2, count2 in hash2arrayx.items():
                        for i in range(0, size1):
                            arr_xvalue2 = hash2arrayx[key2].split(':')
                            arr_yvalue2 = hash2arrayy[key2].split(':')
                            arr_zvalue2 = hash2arrayz[key2].split(':')
                            size2 = len(arr_xvalue2)
                            for j in range(0, size2):
                                dist = math.sqrt((float(arr_xvalue1[i]) - float(arr_xvalue2[j])) ** 2 + (
                                float(arr_yvalue1[i]) - float(arr_yvalue2[j])) ** 2 + (
                                                 float(arr_zvalue1[i]) - float(arr_zvalue2[j])) ** 2)
                                if key1 not in hashclosedist:
                                    hashclosedist[key1] = dist
                                else:
                                    if dist < hashclosedist[key1]:
                                        hashclosedist[key1] = dist

                                if key2 not in hashclosedist:
                                    hashclosedist[key2] = dist
                                else:
                                    if dist < hashclosedist[key2]:
                                        hashclosedist[key2] = dist
                    closest_dist_arr = get_closedist_between_chiana_chainb(hashclosedist)
                    closest_dist_df = pd.DataFrame.from_records(closest_dist_arr,
                                                                columns=["residue_num", "residue_name", "closedist"])
                    interface = [1 if float(i) < cutoff else 0 for i in closest_dist_df["closedist"].values]
                    closest_dist_df["bind"] = interface
                    closest_dist_df.residue_num = closest_dist_df.residue_num.astype(int)
                    closest_dist_df = closest_dist_df.sort_values(by=["residue_num"])
                    closest_dist_df.reset_index(inplace=True, drop=True)
                    closest_dist_df.to_csv(NMR_closedist_file)

def get_closedist_between_chiana_chainb(hashclosedist):
    i = 0
    j = 0
    hashA = {}
    hashB = {}
    closest_dist_arr = []

    for k, v in sorted(hashclosedist.items()):
        if re.search('NEN', k) or re.search('CEN', k):
            continue
        k = k.split(':')
        chain = k[1]
        if chain == 'A':
            i = i + 1
            k1 = ':'.join([str(k[2]), k[3]])
            hashA[k1] = v
        if chain == 'B':
            j = j + 1
            k2 = ':'.join([str(k[2]), k[3]])
            hashB[k2] = v

    for k, v in sorted(hashA.items()):
        if hashB[k] < v:
            k = k.split(':')
            k[1] = shorten(k[1])
            k.extend([v])
            closest_dist_arr.append(k)
        else:
            k = k.split(':')
            k[1] = shorten(k[1])
            k.extend([v])
            closest_dist_arr.append(k)

    return closest_dist_arr


def shorten(x):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    if len(x) % 3 != 0:
        raise ValueError('Input length should be a multiple of three')

    y = ''
    for i in range(int(len(x)/3)):
            y += d[x[3*i:3*i+3]]
    return y