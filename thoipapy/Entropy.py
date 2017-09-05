
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
    #homo_filter_file_handle = open(homo_filter_file, 'r')
    entropy_file=r'Entropy/%s_entropy.txt' %tmp_protein
    entropy_file_handle=open(entropy_file,'w')
    array = []
    #if tmp_protein=="Q6ZRP7":
    try:
        with open(homo_filter_file, "r") as homo_filter_file_handle:
            #array = []
            for line in homo_filter_file_handle:
                array.append(line)
                #print(line)
        position_freq=[]
        for collen in range(0,len(array[0])):
            col_array=[]
            for rowlen in range(0,len(array)):
                col_array.append(array[rowlen][collen])
            #print(col_array)
            position_freq.append(col_array.count(array[0][collen])/len(col_array))
        #print(position_freq)
        for entropy in position_freq:
            entropy_file_handle.write(str(entropy)+"\n")
        #entropy(position_freq))
            #print(rowlen)
        #print(len(array))
        #print(len(array[0]))
    except:
        print("filter number is empty")
