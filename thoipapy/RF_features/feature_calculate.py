from Bio import AlignIO
from Bio.Align import AlignInfo
import csv
import os
import re
import pandas as pd
import scipy as sc
from pandas import Series
import scipy.stats
import thoipapy.mtutiles as mtutils
import subprocess,threading
import glob
import re
import math
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from difflib import SequenceMatcher

def calc_lipophilicity(seq, method = "mean"):
    """ Calculates the average hydrophobicity of a sequence according to the Hessa biological scale.
    """
    # hydrophobicity scale
    hessa_scale = np.array([0.11, -0.13, 3.49, 2.68, -0.32, 0.74, 2.06, -0.6, 2.71,
                            -0.55, -0.1, 2.05, 2.23, 2.36, 2.58, 0.84, 0.52, -0.31,
                            0.3, 0.68])
    # convert to biopython analysis object
    analysed_seq = ProteinAnalysis(seq)
    # biopython count_amino_acids returns a dictionary.
    aa_counts_dict = analysed_seq.count_amino_acids()
    # convert dictionary to array, sorted by aa
    aa_counts_arr = np.array([value for (key, value) in sorted(aa_counts_dict.items())])
    multiplied = aa_counts_arr * hessa_scale
    if method == "mean":
        return multiplied.mean()
    if method == "sum":
        return multiplied.sum()


def mem_a3m_homologous_filter(tmp_lists,test_protein,pathdict,set_,logging):
    ######for single protein filtering
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein= set_["tm_protein_name"]
        tm_start = int(set_["tm_start"])
        tm_end = int(set_["tm_end"])
        homo_a3m_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein
        homo_a3m_file_handle = open(homo_a3m_file, "r")
        homo_filter_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m.mem.uniq.2gaps") % tm_protein
        homo_filter_file_handle = open(homo_filter_file, "w")
        homo_mem_lips_input_file = os.path.join(set_["output_oa3m_homologous"], "%s.mem.lips.input") % tm_protein
        homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
        if os.path.isfile(homo_a3m_file):
            i = 0
            for line in homo_a3m_file_handle:
                tm_str = line[(tm_start - 1):tm_end]
                if i == 0:
                    tm_query = tm_str
                    i = i + 1
                mean_hydrophobicity = calc_lipophilicity(tm_str)
                ratio = SequenceMatcher(None, tm_query, tm_str).ratio()
                if not re.search("-", tm_str) and not re.search("X",
                                                                tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    print("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                gap_num = tm_str.count("-")
                if (gap_num <= 3 and not re.search("X",
                                                   tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    print("{}".format(tm_str), file=homo_filter_file_handle)
                    # homo_filter_file_handle.write(line)
                    # homo_mem_lips_input_file_handle.write(tm_str)
            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()

    # ########for multiple proteins filtering
    else:
        tmp_list_loc = set_["list_of_tmd_start_end"]
        tmp_file_handle = open(tmp_list_loc, 'r')
        for row in tmp_file_handle:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #if tm_protein=="P00533":
            tm_start=int(row.strip().split()[3])
            tm_end = int(row.strip().split()[4])
            #return tm_end
            homo_a3m_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein[0:6]
            if os.path.isfile(homo_a3m_file):
                homo_a3m_file_handle=open(homo_a3m_file,"r")
                homo_filter_file=os.path.join(set_["output_oa3m_homologous"], "%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"],tm_protein)
                homo_filter_file_handle = open(homo_filter_file,"w")
                homo_mem_lips_input_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.mem.lips.input") %(set_["Datatype"],tm_protein)
                homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
                i = 0
                for line in homo_a3m_file_handle:
                    tm_str = line[(tm_start - 1):tm_end]
                    if i == 0:
                        tm_query = tm_str
                        i = i + 1
                    mean_hydrophobicity = calc_lipophilicity(tm_str)
                    ratio = SequenceMatcher(None, tm_query, tm_str).ratio()
                    if not re.search("-", tm_str) and not re.search("X",
                                                                    tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                        print("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                    gap_num = tm_str.count("-")
                    if (gap_num <= 3 and not re.search("X",
                                                       tm_str) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                        print("{}".format(tm_str), file=homo_filter_file_handle)
                        # homo_filter_file_handle.write(line)
                        # homo_mem_lips_input_file_handle.write(tm_str
                homo_filter_file_handle.close()
                homo_mem_lips_input_file_handle.close()




def pssm_calculation(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('start pssm calculation')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        pssm_file = os.path.join(set_["feature_pssm"], "%s.mem.2gap.pssm.csv") % tm_protein
        pssm_file_handle = open(pssm_file, 'w')
        # homo_filter_fasta_file=os.path.join(set_["homologous_folder"],"%s_homo_filtered_aln.fasta") %tm_protein
        #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein
        #set_["homologous_folder"], "a3m/%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"], tm_protein)
        homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.a3m.mem.uniq.2gaps") % (set_["Datatype"],tm_protein)
        mat = []
        writer = csv.writer(pssm_file_handle, delimiter=',', quotechar='"',
                            lineterminator='\n',
                            quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
        writer.writerow(
            ['residue_num', 'residue_name', 'A', 'I', 'L', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'M', 'S', 'T', 'D', 'E', 'R',
             'H', 'K', 'G', 'P'])
        if os.path.isfile(homo_filter_fasta_file):
            with open(homo_filter_fasta_file) as f:
                for line in f.readlines():
            #for line in open(homo_filter_fasta_file).readlines():
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
                aa_num = [column.count('A'), column.count('I'), column.count('L'), column.count('V'), column.count('F'),
                          column.count('W'), column.count('Y'), column.count('N'), column.count('C'), column.count('Q'),
                          column.count('M'), column.count('S'), column.count('T'), column.count('D'), column.count('E'),
                          column.count('R'), column.count('H'), column.count('K'), column.count('G'), column.count('P')]
                # pssmrow = " ".join(str(item) for item in aa_num)
                # pssm_file_handle.write(mat[0][j] +' '+ str(pssmrow) + "\n")
                aa_num.insert(0, mat[0][j])  ###adding the residue name to the second column
                aa_num.insert(0, j + 1)  ##adding the residue number to the first column
                writer.writerow(aa_num)
                column = []
            pssm_file_handle.close()
            


    #####multiple calculation#########
    else:
        for row in tmp_lists:
            #if tm_protein==test_protein:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #pssm_file = os.path.join(set_["feature_pssm"],"NoRedundPro/%s.mem.2gap.pssm.csv") % tm_protein
            pssm_file = os.path.join(set_["feature_pssm"], "%s/%s.mem.2gap.pssm.csv") %(set_["Datatype"],tm_protein)
            if os.path.isfile(pssm_file):
                continue
            else:
                #homo_filter_fasta_file=os.path.join(set_["homologous_folder"],"%s_homo_filtered_aln.fasta") %tm_protein
                #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein
                #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
                homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"],"%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"],tm_protein)
                if os.path.isfile(homo_filter_fasta_file):
                    pssm_file_handle = open(pssm_file, 'w')
                    mat=[]
                    writer = csv.writer(pssm_file_handle, delimiter=',', quotechar='"',
                                        lineterminator='\n',
                                        quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                    writer.writerow(['residue_num','residue_name','A','I','L','V','F','W','Y','N','C','Q','M','S','T' ,'D','E','R','H','K','G','P'])
                    #if os.path.isfile(homo_filter_fasta_file):
                    #for line in open(homo_filter_fasta_file).readlines():
                    with open(homo_filter_fasta_file) as f:
                        for line in f.readlines():
                            if not re.search("^>", line):
                                mat.append(list(line))
                        rowlen = len(mat[0])
                        collen = len(mat)
                        column = []
                    #write 20 amino acids as the header of pssm output file
                    #pssm_file_handle.write(
                        #'residue'+' '+'A' + ' ' + 'I' + ' ' + 'L' + ' ' + 'V' + ' ' + 'F' + ' ' + 'W' + ' ' + 'Y' + ' ' + 'N' + ' ' + 'C' + ' ' + 'Q' + ' ' + 'M' + ' ' + 'S' + ' ' + 'T' + ' ' + 'D' + ' ' + 'E' + ' ' + 'R' + ' ' + 'H' + ' ' + 'K' + ' ' + 'G' + ' ' + 'P' + '\n')
                    for j in range(0, rowlen - 1):
                        for i in range(0, collen):
                            column.append(mat[i][j])
                        aa_num = [column.count('A'), column.count('I'), column.count('L'), column.count('V'), column.count('F'),
                                  column.count('W'), column.count('Y'), column.count('N'), column.count('C'), column.count('Q'),
                                  column.count('M'), column.count('S'), column.count('T'), column.count('D'), column.count('E'),
                                  column.count('R'), column.count('H'), column.count('K'), column.count('G'), column.count('P')]
                        #pssmrow = " ".join(str(item) for item in aa_num)
                        #pssm_file_handle.write(mat[0][j] +' '+ str(pssmrow) + "\n")
                        aa_num.insert(0,mat[0][j])     ###adding the residue name to the second column
                        aa_num.insert(0,j+1)           ##adding the residue number to the first column
                        writer.writerow(aa_num)
                        column = []
                    pssm_file_handle.close()
                    

def entropy_calculation(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('start entropy calculation')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        #for tm_protein in tmp_lists:
        #if tm_protein==test_protein:
        #homo_filter_fasta_file=os.path.join(set_["homologous_folder"],"%s_homo_filtered_aln.fasta") %tm_protein
        #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein
        homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.a3m.mem.uniq.2gaps") % (set_["Datatype"],tm_protein)
        if os.path.isfile(homo_filter_fasta_file):
            entropy_file = os.path.join(set_["feature_entropy"], "%s.mem.2gap.entropy.csv") % tm_protein
            entropy_file_handle = open(entropy_file, 'w')
            mat = []
            with open(homo_filter_fasta_file) as f:
                for line in f.readlines():    
            #for line in open(homo_filter_fasta_file).readlines():
                    if not re.search("^>", line):
                        mat.append(list(line))
                rowlen = len(mat[0])
                collen = len(mat)
                column = []
            writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num","residue_name","Entropy"])
            for j in range(0, rowlen - 1):
                for i in range(0, collen):
                    column.append(mat[i][j])
                column_serie = Series(column)
                p_data = column_serie.value_counts() / len(column_serie)  # calculates the probabilities
                entropy = sc.stats.entropy(p_data)  # input probabilities to get the entropy
                csv_header_for_ncbi_homologous_file = [j+1,mat[0][j],entropy]
                writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_ncbi_homologous_file)
                #entropy_file_handle.write(mat[0][j]+' '+ str(entropy)+'\n')
                column = []
            entropy_file_handle.close()


    ###########multiple calculation###########################
    else:
        for row in tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #homo_filter_fasta_file=os.path.join(set_["homologous_folder"],"%s_homo_filtered_aln.fasta") %tm_protein
            #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tm_protein
            #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
            homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"],"%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"],tm_protein)
            if os.path.isfile(homo_filter_fasta_file):
                #entropy_file = os.path.join(set_["feature_entropy"], "NoRedundPro/%s.mem.2gap.entropy.csv") % tm_protein
                entropy_file = os.path.join(set_["feature_entropy"], "%s/%s.mem.2gap.entropy.csv") %(set_["Datatype"],tm_protein)
                entropy_file_handle = open(entropy_file, 'w')
                mat = []
                with open(homo_filter_fasta_file) as f:
                    for line in f.readlines():
                #for line in open(homo_filter_fasta_file).readlines():
                        if not re.search("^>", line):
                            mat.append(list(line))
                    rowlen = len(mat[0])
                    collen = len(mat)
                    column = []
                writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num","residue_name","Entropy"])
                for j in range(0, rowlen - 1):
                    for i in range(0, collen):
                        column.append(mat[i][j])
                    column_serie = Series(column)
                    p_data = column_serie.value_counts() / len(column_serie)  # calculates the probabilities
                    entropy = sc.stats.entropy(p_data)  # input probabilities to get the entropy
                    csv_header_for_ncbi_homologous_file = [j+1,mat[0][j],entropy]
                    writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                        quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                    writer.writerow(csv_header_for_ncbi_homologous_file)
                    #entropy_file_handle.write(mat[0][j]+' '+ str(entropy)+'\n')
                    column = []
                entropy_file_handle.close()

def coevoluton_calculation_with_freecontact(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('start coevolution calculation by using freecontact')
    freecontact_loc=set_["freecontact_dir"]
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        #for tm_protein in tmp_lists:
        #if tm_protein==test_protein:
        #homo_filter_fasta_file=os.path.join(set_["homologous_folder"],"%s_homo_filtered_aln.fasta") %tm_protein
        homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.a3m.mem.uniq.2gaps") % (set_["Datatype"],tm_protein)
        if os.path.isfile(homo_filter_fasta_file):
            freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], "zpro/%s.mem.2gap.freecontact") % tm_protein
            freecontact_file_handle = open(freecontact_file, 'w')
            exect_str="grep -v '^>' {aln_file} |sed 's/[a-z]//g'|{freecontact} ".format( aln_file=homo_filter_fasta_file,freecontact=freecontact_loc)
            #return exect_str
            class Command(object):
                def __init__(self, cmd):
                    self.cmd = cmd
                    self.process = None

                def run(self, timeout):
                    def target():
                        print('Thread started')
                        #self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
                        subprocess.call(self.cmd, shell=True, stdout=freecontact_file_handle)
                        #print(self.process.communicate())
                        print('Thread finished')

                    thread = threading.Thread(target=target)
                    thread.start()

                    thread.join(timeout)
                    if thread.is_alive():
                        print('Terminating process')
                        self.process.terminate()
                        thread.join()
                    #print(self.process.returncode)

            command = Command(exect_str)
            #command=Command(exect_str)
            command.run(timeout=2000)
            # command.run(timeout=1)
            #command=mtutils.Command(exect_str)
            #command.run(timeout=120)
            logging.info("Output file: %s\n" %freecontact_file)
            freecontact_file_handle.close()



    #######multiple calculation#############
    else:
        for row in tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"], "NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
            homo_filter_fasta_file = os.path.join(set_["output_oa3m_homologous"],"%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"],tm_protein)
            if os.path.isfile(homo_filter_fasta_file):
                #freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], "zpro/NoRedundPro/%s.mem.2gap.freecontact") % tm_protein
                freecontact_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s/%s.mem.2gap.freecontact") %(set_["Datatype"],tm_protein)
                freecontact_file_handle = open(freecontact_file, 'w')
                exect_str = "grep -v '^>' {aln_file} |sed 's/[a-z]//g'|{freecontact} ".format(
                    aln_file=homo_filter_fasta_file, freecontact=freecontact_loc)

                # return exect_str
                class Command(object):
                    def __init__(self, cmd):
                        self.cmd = cmd
                        self.process = None

                    def run(self, timeout):
                        def target():
                            print('Thread started')
                            #self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
                            subprocess.call(self.cmd, shell=True, stdout=freecontact_file_handle)
                            # print(self.process.communicate())
                            print('Thread finished')

                        thread = threading.Thread(target=target)
                        thread.start()

                        thread.join(timeout)
                        if thread.is_alive():
                            print('Terminating process')
                            self.process.terminate()
                            thread.join()
                        #print(self.process.returncode)

                command = Command(exect_str)
                # command=Command(exect_str)
                command.run(timeout=400)
                # command.run(timeout=1)
                # command=mtutils.Command(exect_str)
                # command.run(timeout=120)
                logging.info("Output file: %s\n" % freecontact_file)





def cumulative_co_evolutionary_strength_parser(tmp_lists,test_protein,pathdict,set_,logging):
    """

    Parameters
    ----------
    tmp_lists : list

    pathdict : dict

    set_ : dict
        ......
    logging : logging.Logger
        ...

    Returns
    -------

    """
    logging.info('cumulative co-evolutionary strength parsing')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        #for tm_protein in tmp_lists:
            #if tm_protein == test_protein:
        freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], "zpro/%s.mem.2gap.freecontact") % tm_protein
        #freecontact_file = "/scratch2/zeng/thoipapy/Features/cumulative_coevolution/Q6ZRP7.freecontact"
        dict_di = {}
        dict_mi = {}
        dict_di_new = {}
        dict_mi_new = {}
        dict_di_list = {}
        dict_mi_list = {}
        s = "-"
        if os.path.isfile(freecontact_file):
            cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],
                                                        "zpro/%s.mem.2gap.cumul.coevostr.csv") % tm_protein
            #cumulative_coevlution_strength_file = "/scratch2/zeng/Q6ZRP7.csv"
            cumulative_coevlution_strength_file_handle = open(cumulative_coevlution_strength_file, 'w')
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
                    dict_di_list[residue_pairs[0]].append(residue_pairs[4])
                else:
                    dict_di_list[residue_pairs[0]].append(residue_pairs[4])
                if not residue_pairs[2] in dict_di_list:
                    dict_di_list[residue_pairs[2]] = []
                    dict_di_list[residue_pairs[2]].append(residue_pairs[4])
                else:
                    dict_di_list[residue_pairs[2]].append(residue_pairs[4])
                if not residue_pairs[0] in dict_mi_list:
                    dict_mi_list[residue_pairs[0]] = []
                    dict_mi_list[residue_pairs[0]].append(residue_pairs[5])
                else:
                    dict_mi_list[residue_pairs[0]].append(residue_pairs[5])
                if not residue_pairs[2] in dict_mi_list:
                    dict_mi_list[residue_pairs[2]] = []
                    dict_mi_list[residue_pairs[2]].append(residue_pairs[5])
                else:
                    dict_mi_list[residue_pairs[2]].append(residue_pairs[5])

            tmd_length = int(row.strip().split()[2])
            residue_di = [0] * tmd_length            #the sum of the tmd_length
            residue_mi = [0] * tmd_length
            residue_di_new = [0] * tmd_length         #the sum of the top 8
            residue_mi_new = [0] * tmd_length
            residue_di_max = [0] * tmd_length  # the sum of the top 8
            residue_mi_max = [0] * tmd_length
            for key in dict_di_list:
                residue_di_new[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:8]))
                residue_mi_new[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:8]))
                residue_di_max[int(key) - 1] = sorted(dict_di_list[key], reverse=True)[0]
                residue_mi_max[int(key) - 1] = sorted(dict_mi_list[key], reverse=True)[0]
                # print(str(key)+"corresponding to"+str(dict_di_list[key]))
            # print(residue_di_new)
            dict_di_value_sort = sorted(dict_di.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]
            dict_mi_value_sort = sorted(dict_mi.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]

            for i in range(0, tmd_length):
                res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
                res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
                residue_di[res0] = residue_di[res0] + float(dict_di_value_sort[i][1])
                residue_di[res1] = residue_di[res1] + float(dict_di_value_sort[i][1])
                residue_mi[res0] = residue_mi[res0] + float(dict_mi_value_sort[i][1])
                residue_mi[res1] = residue_mi[res1] + float(dict_mi_value_sort[i][1])

            writer = csv.writer(cumulative_coevlution_strength_file_handle, delimiter=',', quotechar='"',
                                lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "cumu_di", "cumu_mi", "cumu_di_new", "cumu_mi_new","di_max","mi_max"])
            for index in range(len(residue_di)):
                csv_header_for_cumulative_strength_file = [(index + 1), dict_residuenum_residuename[(index + 1)],
                                                           residue_di[index], residue_mi[index], residue_di_new[index],
                                                           residue_mi_new[index],residue_di_max[index],residue_mi_max[index]]
                # writer = csv.writer(cumulative_coevlution_strength_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                # quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_cumulative_strength_file)
            cumulative_coevlution_strength_file_handle.close()
                #cumulative_coevlution_strength_file_handle.write(str(residue_di[index])+"\t"+str(residue_mi[index])+"\n")
        freecontact_file_handle.close()


    # #######multiple calculation#############
    else:
        for row in tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #freecontact_file = os.path.join(set_["feature_cumulative_coevolution"], "zpro/NoRedundPro/%s.mem.2gap.freecontact") % tm_protein
            freecontact_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s/%s.mem.2gap.freecontact") %(set_["Datatype"],tm_protein)
            # freecontact_file = "/scratch2/zeng/thoipapy/Features/cumulative_coevolution/Q6ZRP7.freecontact"
            dict_di = {}
            dict_mi = {}
            dict_di_new = {}
            dict_mi_new = {}
            dict_di_list = {}
            dict_mi_list = {}
            s = "-";
            if os.path.isfile(freecontact_file):
                #cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/NoRedundPro/%s.mem.2gap.cumul.coevostr.csv") % tm_protein
                cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s/%s.mem.2gap.cumul.coevostr.csv") %(set_["Datatype"],tm_protein)
                # cumulative_coevlution_strength_file = "/scratch2/zeng/Q6ZRP7.csv"
                cumulative_coevlution_strength_file_handle = open(cumulative_coevlution_strength_file, 'w')
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
                        dict_di_list[residue_pairs[0]].append(residue_pairs[4])
                    else:
                        dict_di_list[residue_pairs[0]].append(residue_pairs[4])
                    if not residue_pairs[2] in dict_di_list:
                        dict_di_list[residue_pairs[2]] = []
                        dict_di_list[residue_pairs[2]].append(residue_pairs[4])
                    else:
                        dict_di_list[residue_pairs[2]].append(residue_pairs[4])
                    if not residue_pairs[0] in dict_mi_list:
                        dict_mi_list[residue_pairs[0]] = []
                        dict_mi_list[residue_pairs[0]].append(residue_pairs[5])
                    else:
                        dict_mi_list[residue_pairs[0]].append(residue_pairs[5])
                    if not residue_pairs[2] in dict_mi_list:
                        dict_mi_list[residue_pairs[2]] = []
                        dict_mi_list[residue_pairs[2]].append(residue_pairs[5])
                    else:
                        dict_mi_list[residue_pairs[2]].append(residue_pairs[5])

                tmd_length = int(row.strip().split()[2])
                residue_di = [0] * tmd_length
                residue_mi = [0] * tmd_length
                residue_di_new = [0] * tmd_length
                residue_mi_new = [0] * tmd_length
                residue_di_max = [0] * tmd_length  # the sum of the top 8
                residue_mi_max = [0] * tmd_length
                for key in dict_di_list:
                    residue_di_new[int(key) - 1] = sum(map(float, sorted(dict_di_list[key], reverse=True)[0:8]))
                    residue_mi_new[int(key) - 1] = sum(map(float, sorted(dict_mi_list[key], reverse=True)[0:8]))
                    residue_di_max[int(key) - 1] = sorted(dict_di_list[key], reverse=True)[0]
                    residue_mi_max[int(key) - 1] = sorted(dict_mi_list[key], reverse=True)[0]
                    # print(str(key)+"corresponding to"+str(dict_di_list[key]))
                # print(residue_di_new)
                dict_di_value_sort = sorted(dict_di.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]
                dict_mi_value_sort = sorted(dict_mi.items(), key=lambda x: x[1], reverse=True)[0:tmd_length]

                for i in range(0, tmd_length):
                    res0 = int(dict_di_value_sort[i][0].strip().split('-')[0]) - 1
                    res1 = int(dict_di_value_sort[i][0].strip().split('-')[1]) - 1
                    residue_di[res0] = residue_di[res0] + float(dict_di_value_sort[i][1])
                    residue_di[res1] = residue_di[res1] + float(dict_di_value_sort[i][1])
                    residue_mi[res0] = residue_mi[res0] + float(dict_mi_value_sort[i][1])
                    residue_mi[res1] = residue_mi[res1] + float(dict_mi_value_sort[i][1])

                writer = csv.writer(cumulative_coevlution_strength_file_handle, delimiter=',', quotechar='"',
                                    lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num", "residue_name", "cumu_di", "cumu_mi", "cumu_di_new", "cumu_mi_new","di_max","mi_max"])
                for index in range(len(residue_di)):
                    csv_header_for_cumulative_strength_file = [(index + 1), dict_residuenum_residuename[(index + 1)],
                                                               residue_di[index], residue_mi[index], residue_di_new[index],
                                                               residue_mi_new[index],residue_di_max[index],residue_mi_max[index]]
                    # writer = csv.writer(cumulative_coevlution_strength_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                    # quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                    writer.writerow(csv_header_for_cumulative_strength_file)
                cumulative_coevlution_strength_file_handle.close()
                # cumulative_coevlution_strength_file_handle.write(str(residue_di[index])+"\t"+str(residue_mi[index])+"\n")



def Lips_score_calculation(tmp_lists,tm_protein_name, pathdict, set_, logging):
    logging.info('start lips score calculation')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        sequence = ""
        #Lips_input_file = r"/scratch2/zeng/Exp_Pred_161020/a3m/Q9Y286.MEM.lips.input"  ##give the imput alignment without gaps
        Lips_input_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.mem.lips.input") % (set_["Datatype"],tm_protein)
    # for row in tmp_lists:
    #     tmp_protein_acc = row.strip().split("\t")[0]
    #     # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
    #
    #     if len(tmp_protein_acc) == 6:  # uniprot acc
    #         tm_protein = tmp_protein_acc
    #     else:  # PDB proteins
    #         tm_protein = tmp_protein_acc[0:7]
    #     Lips_input_file =os.path.join(set_["output_oa3m_homologous"], "TMPAD/%s.mem.lips.input") % tm_protein
        if os.path.isfile(Lips_input_file):
            file = open(Lips_input_file, "r")
            Lips_output_file = os.path.join(set_["feature_lips_score"], "zpro/%s.mem.lips.output") % tm_protein
            Lips_output_file_handle=open(Lips_output_file,"w")
            sequence = ' '.join(line.strip() for line in file)

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
                print("SURFACE", "%s" % i, ":",file=Lips_output_file_handle)
                #Lips_output_file_handle.write("SURFACE", "%s" % i, ":")
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
                          "%6.3f" % exp_entropy[j],file=Lips_output_file_handle)  # print residue information which is in surface i
                    #Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prop,"%6.3f" % exp_entropy[j])
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
                            print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],file=Lips_output_file_handle)
                            #Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                        k = k + 1
                    j = j + 7
            for i in range(4, 7):  # for surfaces from 4 to 6
                print("SURFACE", "%s" % i, ":",file=Lips_output_file_handle)
                #Lips_output_file_handle.write("SURFACE", "%s" % i, ":")
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
                    print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j],file=Lips_output_file_handle)
                    #Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j])
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
                            print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],file=Lips_output_file_handle)
                            #Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
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
                        print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],file=Lips_output_file_handle)
                        #Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
            print("SURFACE LIPOPHILICITY ENTROPY   LIPS",file=Lips_output_file_handle)
            #Lips_output_file_handle.write("SURFACE LIPOPHILICITY ENTROPY   LIPS")
            for i in sumim.keys():
                avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
                avpim = avpim * 2
                ave = sume[i] / aanum[i]  # average entropy for surface i
                peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
                print("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
                      "%8.3f" % peim,file=Lips_output_file_handle)  # print seven surfaces and see which surface with lowewst LIPS score
                #Lips_output_file_handle.write("%s" % i, "%10.3f" % avpim, "%8.3f" % ave, "%8.3f" % peim)
            Lips_output_file_handle.close()
            file.close()

    else:
        for row in tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #Lips_input_file =os.path.join(set_["output_oa3m_homologous"], "NoRedundPro/%s.mem.lips.input") % tm_protein
            Lips_input_file = os.path.join(set_["output_oa3m_homologous"], "%s/%s.mem.lips.input") %(set_["Datatype"],tm_protein)
            if os.path.isfile(Lips_input_file):
                file = open(Lips_input_file, "r")
                #Lips_output_file = os.path.join(set_["feature_lips_score"], "zpro/NoRedundPro/%s.mem.lips.output") % tm_protein
                Lips_output_file = os.path.join(set_["feature_lips_score"],"zpro/%s/%s.mem.lips.output") %(set_["Datatype"],tm_protein)
                Lips_output_file_handle = open(Lips_output_file, "w")
                sequence = ' '.join(line.strip() for line in file)

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
                    print("SURFACE", "%s" % i, ":", file=Lips_output_file_handle)
                    # Lips_output_file_handle.write("SURFACE", "%s" % i, ":")
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
                              file=Lips_output_file_handle)  # print residue information which is in surface i
                        # Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prop,"%6.3f" % exp_entropy[j])
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
                                      file=Lips_output_file_handle)
                                # Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                            k = k + 1
                        j = j + 7
                for i in range(4, 7):  # for surfaces from 4 to 6
                    print("SURFACE", "%s" % i, ":", file=Lips_output_file_handle)
                    # Lips_output_file_handle.write("SURFACE", "%s" % i, ":")
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
                        print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j], file=Lips_output_file_handle)
                        # Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j])
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
                                      file=Lips_output_file_handle)
                                # Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
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
                            print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k], file=Lips_output_file_handle)
                            # Lips_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                        k = k + 1
                print("SURFACE LIPOPHILICITY ENTROPY   LIPS", file=Lips_output_file_handle)
                # Lips_output_file_handle.write("SURFACE LIPOPHILICITY ENTROPY   LIPS")
                for i in sumim.keys():
                    avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
                    avpim = avpim * 2
                    ave = sume[i] / aanum[i]  # average entropy for surface i
                    peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
                    print("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
                          "%8.3f" % peim,
                          file=Lips_output_file_handle)  # print seven surfaces and see which surface with lowewst LIPS score
                    # Lips_output_file_handle.write("%s" % i, "%10.3f" % avpim, "%8.3f" % ave, "%8.3f" % peim)
                Lips_output_file_handle.close()


def Lips_score_parsing(tmp_lists,tm_protein_name, pathdict, set_, logging):
    logging.info('start parsing lips output to cons and lips scores')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        surface_num = 0
        surface_lips = 100  ##100 is an initialized big number assuming lips score will not bigger than this number
        Lips_output_file = os.path.join(set_["feature_lips_score"], "zpro/%s.mem.lips.output") % tm_protein
        if os.path.isfile(Lips_output_file):
            Lips_outnput_file_handle = open(Lips_output_file, "r")
            conslips_output_file = os.path.join(set_["feature_lips_score"], "zpro/%s.mem.lips.output.conslips.csv") % tm_protein
            conslips_output_file_handle = open(conslips_output_file, "w")
            i = 0
            array = []
            dict = {}
            for row in Lips_outnput_file_handle:
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
            Lips_outnput_file_handle.close()

            surface_find = 0
            dict1 = {}
            Lips_outnput_file_handle = open(Lips_output_file, "r")
            for row in Lips_outnput_file_handle:
                if re.search("^SURFACE\s" + surface_num, row):
                    surface_find = 1
                    continue
                if surface_find == 1 and re.search("^\s+\d+\s+[A-Z]", row):
                    array = row.split()
                    if not int(array[0]) in dict1:
                        dict1[int(array[0])] = " ".join([array[1], array[2], array[3]])
                    #print(row)
                else:
                    surface_find = 0
                    # print(row)
            Lips_outnput_file_handle.close()

            writer = csv.writer(conslips_output_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "lipophilicity", "entropy2","Lips_surface"])
            for k, v in sorted(dict.items()):
                v1 = v.split()
                v1.insert(0, k)
                if k in dict1:
                    v1.insert(4, 1)
                else:
                    v1.insert(4, 0)
                csv_header_for_cons_lips_score_file = v1
                writer.writerow(csv_header_for_cons_lips_score_file)
            conslips_output_file_handle.close()
            

    else:
        for row in tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
            if len(tmp_protein_acc) == 6:  # uniprot acc
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                tm_protein = tmp_protein_acc[0:7]
            #Lips_output_file = os.path.join(set_["feature_lips_score"], "zpro/NoRedundPro/%s.mem.lips.output") % tm_protein
            Lips_output_file = os.path.join(set_["feature_lips_score"],"zpro/%s/%s.mem.lips.output") %(set_["Datatype"],tm_protein)
            if os.path.isfile(Lips_output_file):
                surface_num = 0
                surface_lips = 100  ##100 is an initialized big number assuming lips score will not bigger than this number
                Lips_outnput_file_handle = open(Lips_output_file, "r")
                #conslips_output_file = os.path.join(set_["feature_lips_score"], "zpro/NoRedundPro/%s.mem.lips.output.conslips.csv") % tm_protein
                conslips_output_file = os.path.join(set_["feature_lips_score"],"zpro/%s/%s.mem.lips.output.conslips.csv") %(set_["Datatype"],tm_protein)
                conslips_output_file_handle = open(conslips_output_file, "w")
            i = 0
            array = []
            dict = {}
            for row in Lips_outnput_file_handle:
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

            surface_find = 0
            dict1 = {}
            Lips_outnput_file_handle = open(Lips_output_file, "r")
            for row in Lips_outnput_file_handle:
                if re.search("^SURFACE\s" + surface_num, row):
                    surface_find = 1
                    continue
                if surface_find == 1 and re.search("^\s+\d+\s+[A-Z]", row):
                    array = row.split()
                    if not int(array[0]) in dict1:
                        dict1[int(array[0])] = " ".join([array[1], array[2], array[3]])
                    print(row)
                else:
                    surface_find = 0


            #Lips_outnput_file_handle.close()

            writer = csv.writer(conslips_output_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "lipophilicity", "entropy2","Lips_surface"])
            for k, v in sorted(dict.items()):
                v1 = v.split()
                v1.insert(0, k)
                if k in dict1:
                    v1.insert(4, 1)
                else:
                    v1.insert(4, 0)
                csv_header_for_cons_lips_score_file = v1
                writer.writerow(csv_header_for_cons_lips_score_file)
            Lips_outnput_file_handle.close()




def convert_bind_data_to_csv(tmp_lists,test_protein,pathdict,set_,logging):
    for row in tmp_lists:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

        if len(tmp_protein_acc) == 6:  # uniprot acc
            #continue         ###only considers structure proteins
            tm_protein = tmp_protein_acc
        else:  # PDB proteins
            tm_protein = tmp_protein_acc[0:7]
        #if tm_protein == "5hej_A2":
        #bind_file=os.path.join(set_["structure_bind"],"NoRedundPro/%s.bind") %tm_protein
        #bind_file = os.path.join(set_["structure_bind"], "NoRedundPro/%s.bind") % tm_protein
        bind_file = os.path.join(set_["structure_bind"], "%s/%s.bind.closedist") %(set_["Datatype"],tm_protein)
        if os.path.isfile(bind_file):
            bind_file_handle=open(bind_file,"r")
            #csv_output_file=os.path.join(set_["structure_bind"],"NoRedundPro/%s.csv") %tm_protein
            csv_output_file = os.path.join(set_["structure_bind"], "%s/%s.bind.closedist.csv") %(set_["Datatype"],tm_protein)
            csv_output_file_handle=open(csv_output_file,"w")
            writer = csv.writer(csv_output_file_handle, delimiter=',', lineterminator='\n')
            writer.writerow(["residue_num", "residue_name", "bind","closedist"])
            i=1
            for row in bind_file_handle:
                array=row.split()
                if len(array) ==4:
                 csv_header_for_bind_file=[i,array[0].strip('"'),array[2],array[3]]
                 writer.writerow(csv_header_for_bind_file)
                 i=i+1


                         ###################################################################################################
                         #                                                                                                 #
                         #            combining train data and add physical parameters                                     #
                         #                                                                                                 #
                         ###################################################################################################

def features_combine_to_traindata(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('Combining features into traindata')
    freecontact_loc = set_["freecontact_dir"]
    #for tm_protein in tmp_lists:
    for row in tmp_lists:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

        if len(tmp_protein_acc) == 6:  # uniprot acc
            #continue                   ##only consider structure proteins
            tm_protein = tmp_protein_acc
        else:  # PDB proteins
            tm_protein = tmp_protein_acc[0:7]
        #if tm_protein == "4wit_A1":
        # entropy_file = os.path.join(set_["feature_entropy"], "NoRedundPro/%s.mem.2gap.entropy.csv") % tm_protein
        # pssm_file = os.path.join(set_["feature_pssm"], "NoRedundPro/%s.mem.2gap.pssm.csv") % tm_protein
        # cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/NoRedundPro/%s.mem.2gap.cumul.coevostr.csv") % tm_protein
        # conslips_score_file=os.path.join(set_["feature_lips_score"],"zpro/NoRedundPro/%s.mem.lips.output.conslips.csv") % tm_protein
        # bind_status_file=os.path.join(set_["RF_features"],'Structure/NoRedundPro/%s.csv') %tm_protein
        entropy_file = os.path.join(set_["feature_entropy"], "%s/%s.mem.2gap.entropy.csv") %(set_["Datatype"],tm_protein)
        pssm_file = os.path.join(set_["feature_pssm"], "%s/%s.mem.2gap.pssm.csv") %(set_["Datatype"],tm_protein)
        cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s/%s.mem.2gap.cumul.coevostr.csv") %(set_["Datatype"],tm_protein)
        #cumulative_coevlution_strength_file = os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.cumul.coevostr.csv") %tm_protein
        conslips_score_file=os.path.join(set_["feature_lips_score"],"zpro/%s/%s.mem.lips.output.conslips.csv") %(set_["Datatype"],tm_protein)
        bind_status_file=os.path.join(set_["RF_features"],'Structure/%s/%s.bind.closedist.csv') %(set_["Datatype"],tm_protein)
        #print(bind_status_file)
        if (os.path.isfile(entropy_file) and os.path.isfile(pssm_file) and os.path.isfile(cumulative_coevlution_strength_file) and os.path.isfile(conslips_score_file) and os.path.isfile(bind_status_file)):
            #print(bind_status_file)
            entropy_file_handle = pd.read_csv(entropy_file)
            pssm_file_handle = pd.read_csv(pssm_file)
            cumulative_coevlution_strength_file_handle = pd.read_csv(cumulative_coevlution_strength_file)
            conslips_score_file_handle=pd.read_csv(conslips_score_file)
            bind_status_file_handle=pd.read_csv(bind_status_file)
            #feature_combined_file = os.path.join(set_["RF_features"], "NoRedundPro/%s.mem.2gap.traindata1.csv") % tm_protein
            feature_combined_file=os.path.join(set_["RF_features"],"%s/%s.mem.2gap.traindata1.csv") %(set_["Datatype"],tm_protein)
            #feature_combined_file=os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.traindata1.csv") %tm_protein
            feature_combined_file_handle=open(feature_combined_file,'w')
            merge1=entropy_file_handle.merge(pssm_file_handle,on=['residue_num','residue_name'])
            merge2=merge1.merge(cumulative_coevlution_strength_file_handle,on=["residue_num","residue_name"])
            merge3=merge2.merge(conslips_score_file_handle,on=["residue_num","residue_name"])
            merge4 = merge3.merge(bind_status_file_handle, on=["residue_num", "residue_name"])
            merge4.to_csv(feature_combined_file_handle)
            #writer = csv.writer(feature_combined_file_handle, delimiter=',', quotechar='"',
                                #lineterminator='\n',
                                #quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            #writer.wrterow
        feature_combined_file_handle.close()




def adding_physical_parameters_to_train_data(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('adding physical parameters into traindata')
    for row in tmp_lists:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

        if len(tmp_protein_acc) == 6:  # uniprot acc
            #continue                   ##only consider structure proteins
            tm_protein = tmp_protein_acc
        else:  # PDB proteins
            tm_protein = tmp_protein_acc[0:7]
        #if tm_protein == "5hej_A2":
        dict = {}
        physical_parameter_file=os.path.join(set_["feature_physical_parameters"],"PhysicalProperty.txt")
        physical_parameter_file_handle=open(physical_parameter_file,"r")
        #file_handle = open(r"/scratch2/zeng/Exp_Pred_161020/a3m/PhysicalProperty.txt", "r")
        #train_data_file=os.path.join(set_["RF_features"],"NoRedundP/%s.mem.2gap.traindata1.csv") %tm_protein
        #train_data_file = os.path.join(set_["RF_features"], "%s/%s.mem.2gap.traindata1.csv") %(set_["Datatype"],tm_protein)
        train_data_file = os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact", "%s.traindata1.csv") %tm_protein
        #train_data_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.traindata.csv", "r")
        if os.path.isfile(train_data_file):
            train_data_file_handle=open(train_data_file,"r")
            #train_data_add_physical_parameter_file=os.path.join(set_["RF_features"],"NoRedundPro/%s.mem.2gap.physipara.traindata1.csv") %tm_protein
            train_data_add_physical_parameter_file = os.path.join(set_["RF_features"],"%s/%s.mem.2gap.physipara.traindata1.csv") %(set_["Datatype"],tm_protein)
            #train_data_add_physical_parameter_file = os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.physipara.traindata1.csv") %tm_protein
            train_data_add_physical_parameter_file_handle=open(train_data_add_physical_parameter_file,"w")
            #train_data_physical_parameter_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
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
                    array1[33:14]=["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn"]
                    #array2 = array1[0:31]
                    #array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                    writer.writerow(array1)
                else:
                    array3 = row1.rstrip().split(",")
                    array3[33:14] = dict[array3[2]]
                    writer.writerow(array3)
            train_data_file_handle.close()
        physical_parameter_file_handle.close()



                                             ###################################################################################################
                                             #                                                                                                 #
                                             #            combining test data and add physical parameters                                      #
                                             #                                                                                                 #
                                             ###################################################################################################


def features_combine_to_testdata(test_tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('Combining features into testdata')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        entropy_file = os.path.join(set_["feature_entropy"], "%s.mem.2gap.entropy.csv") % tm_protein
        pssm_file = os.path.join(set_["feature_pssm"], "%s.mem.2gap.pssm.csv") % tm_protein
        cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s.mem.2gap.cumul.coevostr.csv") % tm_protein
        #cumulative_coevlution_strength_file = os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.cumul.coevostr.csv") % tm_protein
        conslips_score_file=os.path.join(set_["feature_lips_score"],"zpro/%s.mem.lips.output.conslips.csv") % tm_protein
        #bind_status_file=os.path.join(set_["RF_features"],'Structure/%s.csv') %tm_protein
        #print(bind_status_file)
        if (os.path.isfile(entropy_file) and os.path.isfile(pssm_file) and os.path.isfile(cumulative_coevlution_strength_file) and os.path.isfile(conslips_score_file)):
            entropy_file_handle = pd.read_csv(entropy_file)
            pssm_file_handle = pd.read_csv(pssm_file)
            cumulative_coevlution_strength_file_handle = pd.read_csv(cumulative_coevlution_strength_file)
            conslips_score_file_handle=pd.read_csv(conslips_score_file)
            #bind_status_file_handle=pd.read_csv(bind_status_file)
            feature_combined_file=os.path.join(set_["RF_loc"],"TestData/%s/%s.mem.2gap.testdata.csv") %(set_["Datatype"],tm_protein)
            #feature_combined_file=os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.testdata.csv") %tm_protein
            feature_combined_file_handle=open(feature_combined_file,'w')
            merge1=entropy_file_handle.merge(pssm_file_handle,on=['residue_num','residue_name'])
            merge2=merge1.merge(cumulative_coevlution_strength_file_handle,on=["residue_num","residue_name"])
            merge3=merge2.merge(conslips_score_file_handle,on=["residue_num","residue_name"])
            #merge4 = merge3.merge(bind_status_file_handle, on=["residue_num", "residue_name"])
            merge3.to_csv(feature_combined_file_handle)
            #writer = csv.writer(feature_combined_file_handle, delimiter=',', quotechar='"',
                                #lineterminator='\n',
                                #quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
        feature_combined_file_handle.close()
            #writer.wrterow
    else:
        for row in test_tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                #continue                   ##only consider structure proteins
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                continue                #only considers test data
                #tm_protein = tmp_protein_acc[0:7]
            #if tm_protein == "Q7L4S7":
            entropy_file = os.path.join(set_["feature_entropy"], "%s/%s.mem.2gap.entropy.csv") %(set_["Datatype"],tm_protein)
            pssm_file = os.path.join(set_["feature_pssm"], "%s/%s.mem.2gap.pssm.csv") %(set_["Datatype"],tm_protein)
            #cumulative_coevlution_strength_file = os.path.join(set_["feature_cumulative_coevolution"],"zpro/%s/%s.mem.2gap.cumul.coevostr.csv") %(set_["Datatype"],tm_protein)
            cumulative_coevlution_strength_file = os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.cumul.coevostr.csv") % tm_protein
            conslips_score_file=os.path.join(set_["feature_lips_score"],"zpro/%s/%s.mem.lips.output.conslips.csv") %(set_["Datatype"],tm_protein)
            #bind_status_file=os.path.join(set_["RF_features"],'Structure/%s.csv') %tm_protein
            #print(bind_status_file)
            if (os.path.isfile(entropy_file) and os.path.isfile(pssm_file) and os.path.isfile(cumulative_coevlution_strength_file) and os.path.isfile(conslips_score_file)):
                entropy_file_handle = pd.read_csv(entropy_file)
                pssm_file_handle = pd.read_csv(pssm_file)
                cumulative_coevlution_strength_file_handle = pd.read_csv(cumulative_coevlution_strength_file)
                conslips_score_file_handle=pd.read_csv(conslips_score_file)
                #bind_status_file_handle=pd.read_csv(bind_status_file)
                #feature_combined_file=os.path.join(set_["RF_loc"],"TestData/%s/%s.mem.2gap.testdata.csv") %(set_["Datatype"],tm_protein)
                feature_combined_file=os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.testdata.csv") %tm_protein
                feature_combined_file_handle=open(feature_combined_file,'w')
                merge1=entropy_file_handle.merge(pssm_file_handle,on=['residue_num','residue_name'])
                merge2=merge1.merge(cumulative_coevlution_strength_file_handle,on=["residue_num","residue_name"])
                merge3=merge2.merge(conslips_score_file_handle,on=["residue_num","residue_name"])
                #merge4 = merge3.merge(bind_status_file_handle, on=["residue_num", "residue_name"])
                merge3.to_csv(feature_combined_file_handle)
                #writer = csv.writer(feature_combined_file_handle, delimiter=',', quotechar='"',
                                    #lineterminator='\n',
                                    #quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                #writer.wrterow




def adding_physical_parameters_to_test_data(test_tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('adding physical parameters into testdata')
    if not set_["multiple_tmp_simultaneous"]:
        tm_protein = set_["tm_protein_name"]
        dict = {}
        physical_parameter_file = os.path.join(set_["feature_physical_parameters"], "PhysicalProperty.txt")
        physical_parameter_file_handle = open(physical_parameter_file, "r")
        # file_handle = open(r"/scratch2/zeng/Exp_Pred_161020/a3m/PhysicalProperty.txt", "r")
        train_data_file = os.path.join(set_["RF_loc"], "TestData/%s/%s.mem.2gap.testdata.csv") %(set_["Datatype"], tm_protein)
        # train_data_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.traindata.csv", "r")
        if os.path.isfile(train_data_file):
            train_data_file_handle = open(train_data_file, "r")
            train_data_add_physical_parameter_file = os.path.join(set_["RF_loc"],
                                                                  "TestData/%s/%s.mem.2gap.physipara.testdata.csv") %(set_["Datatype"], tm_protein)
            train_data_add_physical_parameter_file_handle = open(train_data_add_physical_parameter_file, "w")
            # train_data_physical_parameter_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
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
                    array1[33:14] = ["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic",
                                     "Aromatic", "Polar", "Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn"]
                    # array2 = array1[0:31]
                    # array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                    writer.writerow(array1)
                else:
                    array3 = row1.rstrip().split(",")
                    array3[33:14] = dict[array3[2]]
                    writer.writerow(array3)
            train_data_file_handle.close()
            train_data_add_physical_parameter_file_handle.close()
        physical_parameter_file_handle.close()


    else:
        for row in test_tmp_lists:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                #continue                   ##only consider structure proteins
                tm_protein = tmp_protein_acc
            else:  # PDB proteins
                continue            ##onlys considers test proteins
                tm_protein = tmp_protein_acc[0:7]
            #if tm_protein == "Q7L4S7":
            dict = {}
            physical_parameter_file=os.path.join(set_["feature_physical_parameters"],"PhysicalProperty.txt")
            physical_parameter_file_handle=open(physical_parameter_file,"r")
            #file_handle = open(r"/scratch2/zeng/Exp_Pred_161020/a3m/PhysicalProperty.txt", "r")
            #train_data_file=os.path.join(set_["RF_loc"],"TestData/%s/%s.mem.2gap.testdata.csv") %(set_["Datatype"], tm_protein)
            train_data_file=os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact","%s.testdata.csv") %tm_protein
            #train_data_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.traindata.csv", "r")
            if os.path.isfile(train_data_file):
                train_data_file_handle=open(train_data_file,"r")
                #train_data_add_physical_parameter_file=os.path.join(set_["RF_loc"],"TestData/%s/%s.mem.2gap.physipara.testdata.csv") %(set_["Datatype"], tm_protein)
                train_data_add_physical_parameter_file=os.path.join("/scratch2/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact/TestData","%s.physipara.testdata.csv") %tm_protein
                train_data_add_physical_parameter_file_handle=open(train_data_add_physical_parameter_file,"w")
                #train_data_physical_parameter_file_handle = open(r"/scratch2/zeng/thoipapy/Features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
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
                        array1[33:14]=["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn"]
                        #array2 = array1[0:31]
                        #array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                        writer.writerow(array1)
                    else:
                        array3 = row1.rstrip().split(",")
                        array3[33:14] = dict[array3[2]]
                        writer.writerow(array3)



def combine_all_train_data_for_random_forest(tmp_lists,test_protein,pathdict,set_,logging):
    logging.info('creating train and test data for random forest')
    # path = r"/scratch/zeng/thoipapy/Features/"
    # interesting_files = glob.glob("/scratch/zeng/thoipapy/RandomForest/TrainData/*.csv")
    #train_data_add_physical_parameter_files=glob.glob(os.path.join(set_["RF_features"],"TMPAD/????-??.mem.2gap.physipara.traindata.csv"))
    #train_data_add_physical_parameter_files = glob.glob(os.path.join(set_["RF_features"], "NoRedundPro/*.mem.2gap.physipara.traindata1.csv"))
    train_data_add_physical_parameter_files = glob.glob(os.path.join(set_["RF_features"], "%s/*.mem.2gap.physipara.traindata1.csv") %set_["Datatype"])
    #train_data_add_physical_parameter_files = glob.glob(os.path.join("/scratch/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact/Train68", "*.physipara.traindata2.csv") )
    #train_data_add_physical_parameter_files = glob.glob(os.path.join(set_["RF_features"], "%s/TrainData1_68/7Resd/*.mem.2gap.physipara.traindata7.csv") %set_["Datatype"])
    #interesting_files = glob.glob("/scratch2/zeng/thoipapy/RandomForest/TestData/*.csv")
    header_saved = False
    # with open('/scratch/zeng/thoipapy/RandomForest/traindata.csv','w') as fout:
    #with open('/scratch/zeng/thoipapy/RandomForest/pb.traindata.csv', 'w') as fout:
    #with open('/scratch/zeng/thoipapy/RandomForest/NoRedundPro/TRAINDATA_orig1.csv', 'w') as fout:
    with open('/scratch/zeng/thoipapy/RandomForest/%s/TRAINDATA_orig1.csv' %set_["Datatype"], 'w')   as fout :
    #with open('/scratch/zeng/thoipapy/Features/cumulative_coevolution/zfullfreecontact/Train68/TRAINDATA_orig2.csv' , 'w')   as fout :
        for filename in train_data_add_physical_parameter_files:
            with open(filename) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)