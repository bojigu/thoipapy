import math
import os
import re
import sys
from pathlib import Path

import pandas as pd
import thoipapy

def calc_closedist_from_PREDDIMER_TMDOCK_best_model(s, df_set, logging):
    """Calculate the closest heavy atom distances from PREDDIMER or TMDOCK structures for list of TMDs.

    See closedist_calculate_from_dimer function for details.
        This function simply runs closedist_calculate_from_dimer on each TMD in a list.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    """
    logging.info("starting calc_closedist_from_PREDDIMER_TMDOCK_best_model for {}".format(s["setname"]))

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        sys.stdout.write("{}, ".format(acc)), sys.stdout.flush()
        database = df_set.loc[i, "database"]
        protein = acc
        other_predictors_dir = Path(s["thoipapy_data_folder"]) / "Predictions/other_predictors"
        pdb_file_preddimer = os.path.join(other_predictors_dir,database,"{}.preddimer.pdb".format(protein))
        pdb_file_tmdock = os.path.join(other_predictors_dir, database, "{}.tmdock.pdb".format(protein))
        preddimer_closedist_file = os.path.join(other_predictors_dir, database,"{}.preddimer.closedist.csv".format(protein))
        tmdock_closedist_file = os.path.join(other_predictors_dir, database,"{}.tmdock.closedist.csv".format(protein))

        #closedist_calculate_from_dimer(s,pdb_file_preddimer,preddimer_closedist_file)
        if os.path.isfile(pdb_file_tmdock):
            closedist_calculate_from_dimer(acc, s, logging, pdb_file_tmdock, tmdock_closedist_file)
        else:
            logging.warning("{} : TMDOCK pdb file not found".format(protein))
            continue

        if os.path.isfile(pdb_file_preddimer):
            closedist_calculate_from_dimer(acc, s, logging, pdb_file_preddimer, preddimer_closedist_file)
        else:
            logging.warning("{} : PREDDIMER pdb file not found".format(protein))
            continue

        if i % 10 == 0:
            sys.stdout.write("\n")

    sys.stdout.write("\n")

def closedist_calculate_from_dimer(acc, s, logging, pdb_file, closedist_out_csv):
    """Calculate the closest heavy atom distance between TMD residues from top-ranked PREDDIMER or TMDOCK structure.

    Distances are calculated between all heavy atoms.

    For each residue in the TMD, finds the closest heavy atom distance to any other residue in the opposing chain.

    This is used to define "interface" residues according to the PREDDIMER/TMDOCK prediction, in exactly the same way
    as the experimental crystal or NMR structures.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    pdb_file : str
        Path to pdb file that is the output from PREDDIMER or TMDOCK.
    closedist_out_csv : str
        Path to output file
    """
    if not os.path.isfile(pdb_file):
        logging.warning("the pdb file: {} is not existed, and skipped\n, this could be the tmdock prediction for this protein didn't exists, please try to run TMDOCK"
                         " server again with different sequence TMD length input".format(pdb_file))
        return None
    hash1arrayx = {}
    hash1arrayy = {}
    hash1arrayz = {}
    hash2arrayx = {}
    hash2arrayy = {}
    hash2arrayz = {}
    hashclosedist = {}

    pdb = pdb_file.strip().split('\\')[-1][0:6]
    chain1 = 'A'
    chain2 = 'B'

    with open(pdb_file, "r") as f:
        for line in f:
            if re.search(r"^MODEL\s+2", line):
                break
            if re.search(r"^ATOM", line):
                atom = line[12:16]
                #sys.stdout.write("{},".format(atom))
                # skip any hydrogens
                if re.search(r"^\s*H", atom):  # non-H atom distance
                    continue
                index = line[6:11]
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                chain = line[21:22]
                residue_num = line[22:26]
                residue_name = line[17:20]
                # skip any "OCT" atoms, which are C-terminal atoms that are often erroneously labelled by TMDOCK as
                # part of a new residue (e.g. ATALFVLVPSVFLIILYV of 2axt chain M, TMD1, where residue 18 is actually a Val
                # but  "275  OCT TYR B  18" suggests a TYR at positions 17 and 18
                if atom == " OCT" and residue_name != "CEN":
                    sys.stdout.write("\n")
                    logging.warning("{}  possible mislabelling of terminal residue (OCT for residue {})\nFull line:\n{}".format(acc, residue_name, line))
                    # skip this line
                    continue
                if chain == chain1:
                    kk = ':'.join([pdb, chain, residue_num, residue_name])
                    if kk not in hash1arrayx:
                        hash1arrayx[kk] = x
                        hash1arrayy[kk] = y
                        hash1arrayz[kk] = z
                    else:
                        hash1arrayx[kk] = hash1arrayx[kk] + ':' + x
                        hash1arrayy[kk] = hash1arrayy[kk] + ':' + y
                        hash1arrayz[kk] = hash1arrayz[kk] + ':' + z

                if chain == chain2:
                    kk = ':'.join([pdb, chain, residue_num, residue_name])
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
    closest_dist_arr = get_closedist_between_chainA_and_chainB(hashclosedist)
    closest_dist_df = pd.DataFrame.from_records(closest_dist_arr, columns = ["residue_num","residue_name","closedist"])
    closest_dist_df.residue_num = closest_dist_df.residue_num.astype(int)
    closest_dist_df=closest_dist_df.sort_values(by=["residue_num"])
    closest_dist_df.reset_index(inplace=True,drop=True)
    closest_dist_df.to_csv(closedist_out_csv)

def get_closedist_between_chainA_and_chainB(hashclosedist):
    i = 0
    j = 0
    hashA = {}
    hashB = {}
    closest_dist_arr = []

    for k, v in sorted(hashclosedist.items()):
        if re.search(r'NEN', k) or re.search(r'CEN', k):
            continue
        k = k.split(':')
        chain = k[1]
        if chain == 'A':
            i = i + 1
            k1 = ':'.join([str(i), k[3]])
            hashA[k1] = v
        if chain == 'B':
            j = j + 1
            k2 = ':'.join([str(j), k[3]])
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