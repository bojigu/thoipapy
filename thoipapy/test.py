import os
import re
import csv
from pandas import Series
import scipy as sc
import scipy.stats

tm_protein="P0A6S5"
homo_filter_fasta_file="give your alignment file path"
if os.path.isfile(homo_filter_fasta_file):
    entropy_file = "give the entropy file path that you will save the entropy result" 
    entropy_file_handle = open(entropy_file, 'w')
    mat = []
    for line in open(homo_filter_fasta_file).readlines():
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
