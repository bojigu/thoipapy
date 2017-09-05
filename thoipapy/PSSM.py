from Bio import AlignIO
from Bio.Align import AlignInfo

f=open('zTmpPosCat.txt', 'r')
for row in f.readlines():
    tmp_protein_num=row.strip().split("\t")[0]
    if len(tmp_protein_num)==7:
            tmp_protein=tmp_protein_num[0:6]
    else:
        tmp_protein=tmp_protein_num
    #if tmp_protein=='Q6ZRP7':
    #print(row)
    homo_filter_file=r"homologous/filter/%s_filter.txt" %tmp_protein
    homo_filter_file_handle = open(homo_filter_file, 'r')
    homo_filter_fasta_file=r"homologous/filter/%s_filter.fasta" %tmp_protein
    homo_filter_fasta_file_handle = open(homo_filter_fasta_file, 'w')
    pssm_file=r'PSSM/%s_pssm.txt' %tmp_protein
    pssm_file_handle=open(pssm_file,'w')
    i=0
    for line in homo_filter_file_handle.readlines():
        homo_filter_fasta_file_handle.write(">"+str(i)+"\n"+line)
        i=i+1
    print(homo_filter_fasta_file)
    try:
        alignment = AlignIO.read(homo_filter_fasta_file, "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        my_pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore = ['N', '-'])
        print(my_pssm)
    except:
        print("pssm calculation error")