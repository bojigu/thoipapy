from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import glob
import os
import sys
import csv
import pandas as pd
from time import strftime
import signal
import unicodedata
import re

proteins=[]
f=open('zTmpPosCat.txt', 'r')
for row in f.readlines():
    tmp_protein_num=row.strip().split("\t")[0]
    if len(tmp_protein_num)==7:
            tmp_protein=tmp_protein_num[0:6]
            tm_start=int(row.strip().split("\t")[3])
            tm_end=int(row.strip().split("\t")[4])
    else:
        tmp_protein=tmp_protein_num
        tm_start=int(row.strip().split("\t")[3])
        tm_end=int(row.strip().split("\t")[4])
    #print tmp_protein,"\t",tm_start,"\t",tm_end
    if tmp_protein=="Q6ZRP7":
        homo_out_file=r"homologous/homo1/%s_homo.txt" %tmp_protein
        homo_out_file_handel=open(homo_out_file, "w")
        #with open(homo_out_file, "w") as homo_out_file_handel:
        try:
            xml_file=r"xml/xml1/%s.xml" %tmp_protein
            xml_result_handle = open(xml_file)
            blast_record = NCBIXML.read(xml_result_handle)
            E_VALUE_THRESH = 0.01
            i=0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        i=i+1
                        query_start=hsp.query_start
                        query_end=hsp.query_end
                        query_seq=hsp.query
                        sbjt_seq=hsp.sbjct
                        query_seq_no_gap=re.sub('-','',query_seq)
                        match_seq=hsp.match
                        tmstr=''
                        sbjstr=''
                        if query_start<=tm_start and query_end>=tm_end:
                            homo_out_file_handel.write("****Alignment****"+str(i)+"\n")
                            print('****Alignment****',i)
                            #print('>', alignment.title)
                            #print('length:', alignment.length)
                            #print('identity:', hsp.identities)
                            print('query_start:', hsp.query_start)
                            print('query_end:', hsp.query_end)
                            print('subject_start:', hsp.sbjct_start)
                            print('subject_end:', hsp.sbjct_end)
                            print('tm_start:', tm_start)
                            print("tm_end",tm_end)
                            homo_out_file_handel.write("e value:" + str(hsp.expect) + "\n")
                            print('e value:', hsp.expect)
                            print(hsp.query)
                            print(hsp.match)
                            print(hsp.sbjct)
                            tm_str_start=tm_start-query_start
                            tm_str_end=tm_end-query_start+1
                            k=0
                            j=0
                            tmstr=""
                            sbjtstr=""
                            matchstr=""
                            for char in query_seq:
                                if char != '-':
                                    if j >= tm_str_start and j <tm_str_end:
                                        tmstr+=query_seq_no_gap[j]
                                        sbjtstr+=sbjt_seq[k]
                                        matchstr+=match_seq[k]
                                    j=j+1
                                k=k+1
                            #homo_out_file_handel.write(str(tmstr) + "\n")
                            #homo_out_file_handel.write(str(sbjtstr) + "\n")
                            #homo_out_file_handel.write(str(matchstr) + "\n")
                            print(tmstr)
                            print(sbjtstr)
                            print(matchstr)
        except:
            print("tmp file not existed")