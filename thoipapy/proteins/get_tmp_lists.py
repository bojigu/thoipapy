import re
import sys
from Bio import SeqIO


def calculate_fasta_file_length(set_):
    """
    :param set_:
    :return: tmp_list
    """
    seqLen=int()
    FastaFile = open(set_["nput_fasta_file"], 'r')
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
    FastaFile.close()
    return seqLen


def extract_tmps_from_input_file(set_):
    """
    :param set_:
    :return: tmp_list
    """
    tmp_lists={}
    tmp_list_loc = set_["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    for row in tmp_file_handle:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
        if len(tmp_protein_acc) == 6:  # uniprot acc
            acc = tmp_protein_acc
            tmp_lists[acc]=1
        else:  # PDB proteins
            acc = tmp_protein_acc[0:7]
            tmp_lists[acc] = 1
    tmp_file_handle.close()
    return tmp_lists

def extract_test_tmps_from_input_file(set_):
    """
    :param set_:
    :return: test tmp lists
    """
    tmp_lists={}
    tmp_list_loc = set_["list_of_test_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    for row in tmp_file_handle:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
        if len(tmp_protein_acc) == 6:  # uniprot acc
            acc = tmp_protein_acc
            tmp_lists[acc]=1
        else:  # PDB proteins
            acc = tmp_protein_acc[0:7]
            tmp_lists[acc] = 1
    tmp_file_handle.close()
    return tmp_lists