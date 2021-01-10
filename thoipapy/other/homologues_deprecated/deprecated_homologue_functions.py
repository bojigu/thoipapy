import os
import re
import sys
from difflib import SequenceMatcher
import pandas as pd
import thoipapy

# replace print function, so a search for "print (" can find areas being debugged
p_r_i_n_t = print

def extract_filtered_csv_homologues_to_alignments_orig_handle_method(s,logging):

    logging.info('start extract filtered csv homologues to alignments')
    min_identity_of_TMD_seq=s["min_identity_of_TMD_seq"]
    max_identity_of_TMD_seq=s["max_identity_of_TMD_seq"]
    max_n_gaps_in_TMD_subject_seq=s["max_n_gaps_in_TMD_subject_seq"]
    tmp_list_loc = s["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for line in tmp_file_handle:
        acc = line.strip().split(",")[0]
        database = line.strip().split(",")[6]
        #if acc == "2j58_C":
        alignment_dict={}
        alignment_dict1 = {}
        alignment_dict2 = {}
        #homo_csv_file_loc = os.path.join(s["data_dir"], "homologues", "NoRedundPro/%s_homo.csv") % acc
        homo_csv_file_loc = os.path.join(s["data_dir"], "homologues","a3m",database, "%s_homo%s.csv") % (acc,s["surres"])
        if (os.path.isfile(homo_csv_file_loc) and os.path.getsize(homo_csv_file_loc) > 0):  # whether protein homologues csv file exists
            df = pd.read_csv(homo_csv_file_loc)
            #homo_filtered_out_csv_file = os.path.join(s["data_dir"], "homologues", "NoRedundPro/%s_homo_filtered_aln.fasta") % acc
            homo_filtered_out_csv_file = os.path.join(s["data_dir"], "homologues","a3m",database,
                                                      "%s_homo_filtered_aln%s.fasta") % (acc,s["surres"])
            homo_filtered_out_csv_file_handle=open(homo_filtered_out_csv_file, 'w')
            #homo_filter_file = os.path.join(s["data_dir"], "homologues", "a3m/NoRedundPro/%s.a3m.mem.uniq.2gaps") % acc
            homo_filter_file = os.path.join(s["data_dir"], "homologues","a3m",database, "%s.a3m.mem.uniq.2gaps%s") % (acc,s["surres"])
            homo_filter_file_handle = open(homo_filter_file, "w")
            #homo_mem_lips_input_file = os.path.join(s["data_dir"], "homologues","a3m/NoRedundPro/%s.mem.lips.input") % acc
            homo_mem_lips_input_file = os.path.join(s["data_dir"], "homologues","a3m",database,"%s.mem.lips.input%s") % (acc,s["surres"])
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            #with open(homo_filtered_out_csv_file, 'w') as homo_filtered_out_csv_file_handle:

            i=0
            try:
                """WARNING : CURRENTLY Q12983 (BNIP3) for ETRA gives an error here
                 - need to re-organise script so the sequences are simply saved based on a filtered dataframe
                
                KeyError: 'tm_query_seq'
                Exception ignored in: <_io.FileIO name='D:\\\\data_thoipapy\\\\Input_data\\set12_BNIP3_ETRA_processed.csv' mode='rb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\\\data_thoipapy\\\\Input_data\\set12_BNIP3_ETRA_processed.csv' mode='r' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983_homo_filtered_aln_surr0.fasta' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983_homo_filtered_aln_surr0.fasta' mode='w' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983.a3m.mem.uniq.2gaps_surr0' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983.a3m.mem.uniq.2gaps_surr0' mode='w' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983.mem.lips.input_surr0' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologues\\a3m\\ETRA\\Q12983.mem.lips.input_surr0' mode='w' encoding='cp1252'>
                """


                for index, row in df.iterrows():
                    if i==0:
                        homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                        homo_filtered_out_csv_file_handle.write(row["tm_query_seq"] + "\n")
                        p_r_i_n_t("{}".format(row["tm_query_seq"]), file=homo_mem_lips_input_file_handle)
                        p_r_i_n_t("{}".format(row["tm_query_seq"]), file=homo_filter_file_handle)
                        alignment_dict[row["tm_query_seq"]] = 1
                        i=1
                        continue
                    tm_query_seq_len = len(row["tm_query_seq"])
                    tm_sbjt_seq_gap_num = row["tm_sbjt_seq"].count('-')
                    query_sbjt_identity_num = [x == y for (x, y) in zip(row["tm_query_seq"], row["tm_sbjt_seq"])].count(True)
                    mean_hydrophobicity= thoipapy.common.calc_lipophilicity(row["tm_sbjt_seq"])
                    ratio = SequenceMatcher(None, row["tm_query_seq"], row["tm_sbjt_seq"]).ratio()
                    gap_num = row["tm_sbjt_seq"].count("-")
                    if not re.search(r"-", row["tm_sbjt_seq"]) and not re.search(r"X",row["tm_sbjt_seq"]) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"] :  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict:
                            alignment_dict[row["tm_sbjt_seq"]]=1
                            homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                            p_r_i_n_t("{}".format(row["tm_sbjt_seq"]), file=homo_filtered_out_csv_file_handle)
                    if not re.search(r"-", row["tm_sbjt_seq"]) and not re.search(r"X",row["tm_sbjt_seq"]) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"] :  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict1:
                            alignment_dict1[row["tm_sbjt_seq"]]=1
                            p_r_i_n_t("{}".format(row["tm_sbjt_seq"]), file=homo_mem_lips_input_file_handle)
                    if (gap_num <= 2 and not re.search(r"X",row["tm_sbjt_seq"]) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"] ):  # gap number le 3 and no X in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict2:
                            alignment_dict2[row["tm_sbjt_seq"]]=1
                            p_r_i_n_t("{}".format(row["tm_sbjt_seq"]), file=homo_filter_file_handle)
            except:
                logging.warning("\n------------------\n{} parsing error in extract_filtered_csv_homologues_to_alignments_mult_prot\n{} could not be created. Files will be deleted.\n------------------\n".format(acc, homo_filter_file))
                homo_filter_file_handle.close()
                homo_mem_lips_input_file_handle.close()
                homo_filtered_out_csv_file_handle.close()
                if os.path.isfile(homo_filtered_out_csv_file):
                    os.remove(homo_filtered_out_csv_file)
                    os.remove(homo_filter_file)
                    os.remove(homo_mem_lips_input_file)

            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()
            homo_filtered_out_csv_file_handle.close()
            logging.info("{} extract_filtered_csv_homologues_to_alignments_mult_prot finished ({})".format(acc, homo_filter_file))
        else:
            sys.stdout.write("{} not found".format(homo_csv_file_loc))
    tmp_file_handle.close()