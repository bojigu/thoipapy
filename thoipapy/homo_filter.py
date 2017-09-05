###for all files in /homologous/*.txt, do such following filtering###
###filter 1: any homologous have less than 3 gaps###
###filter 2: sequence identity is between: 40%-1####
###:homologous have hydrophobicity cutoff :12####

import glob
import os
import sys
import csv
import pandas as pd
from time import strftime
import signal
import unicodedata
import re
f=open('zTmpPosCat.txt', 'r')
for row in f.readlines():
    tmp_protein_num=row.strip().split("\t")[0]
    if len(tmp_protein_num)==7:
            tmp_protein=tmp_protein_num[0:6]
    else:
        tmp_protein=tmp_protein_num
    #if tmp_protein=='Q6ZRP7':
    homo_file=r"homologous/%s_homo.txt" %tmp_protein
    homo_file_handle = open(homo_file, 'r')
    homo_filter_file=r'homologous/filter/%s_filter.txt' %tmp_protein
    homo_filter_file_handle=open(homo_filter_file,'w')
    i=0
    query_seq=''
    sbjt_seq=''
    match_seq=''
    alignment_dict={}
    for line in homo_file_handle.readlines():
        if re.search(r'Alignment',line):
            i=0
        if i==2:
            query_seq=line
            if not query_seq in alignment_dict:
                homo_filter_file_handle.write(query_seq)
                alignment_dict[query_seq]=1

        if i==3:
            sbjt_seq=line
        if i==4:
            match_seq=line
        i=i+1
        if i==5:
            query_seq_len=len(query_seq)
            sbjt_seq_gap_num=sbjt_seq.count('-')
            query_sbjt_identity_num=[x == y for (x, y) in zip(query_seq, sbjt_seq)].count(True)
            if query_seq_len*0.4 <= query_sbjt_identity_num < query_seq_len:
                if sbjt_seq_gap_num <=2:
                    if not sbjt_seq in alignment_dict:
                        alignment_dict[sbjt_seq]=1
    alignment_dict.pop(query_seq, None)
    for key in alignment_dict:
        homo_filter_file_handle.write(key)
        #if sbjt_seq_gap_num<=2:

    #print(query_seq)