def extract_tmps_from_input_file(pathdict,set_,logging):
    tmp_lists={}
    tmp_list_loc = set_["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    for row in tmp_file_handle:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
        if len(tmp_protein_acc) == 6:  # uniprot acc
            tmp_protein_name = tmp_protein_acc
            tmp_lists[tmp_protein_name]=1
        else:  # PDB proteins
            tmp_protein_name = tmp_protein_acc[0:7]
            tmp_lists[tmp_protein_name] = 1
    tmp_file_handle.close()
    return tmp_lists

def extract_test_tmps_from_input_file(pathdict,set_,logging):
    tmp_lists={}
    tmp_list_loc = set_["list_of_test_tmps_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    for row in tmp_file_handle:
        tmp_protein_acc = row.strip().split("\t")[0]
        # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
        if len(tmp_protein_acc) == 6:  # uniprot acc
            tmp_protein_name = tmp_protein_acc
            tmp_lists[tmp_protein_name]=1
        else:  # PDB proteins
            tmp_protein_name = tmp_protein_acc[0:7]
            tmp_lists[tmp_protein_name] = 1
    tmp_file_handle.close()
    return tmp_lists