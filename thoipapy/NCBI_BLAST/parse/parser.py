import csv
import os
from Bio.Blast import NCBIXML
import re
import pandas as pd
from difflib import SequenceMatcher
import thoipapy

def parse_NCBI_xml_to_csv(set_, df_set, logging):
    """
    extract multiple sequence alignment and useful informations from xml homologous file
    Parameters
    ----------
    set_
    logging

    Returns
    -------
    csv file of query protein alignment information

    """
    logging.info('~~~~~~~~~~~~                 starting parse_NCBI_xml_to_csv              ~~~~~~~~~~~~')

    ##############################################################################################
    #                                                                                            #
    #     parse multiple csv files simultaneously                                                #
    #                                                                                            #
    ##############################################################################################


    tmp_list_loc = set_["list_of_tmd_start_end"]

    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0]
        database = row.strip().split(",")[6]
        #tm_start=int(row.strip().split(",")[4])+1
        #tm_start=int(row.strip().split(",")[4])+1-5 ###for surr20
        # set_["surres"] set how many residues on each sides of tmd sequence
        if set_["surres"] == "_surr0":
            tm_start = int(row.strip().split(",")[2])
            tm_end = int(row.strip().split(",")[3])
        elif set_["surres"] == "_surr5":
            tm_start = int(row.strip().split(",")[2]) -5 ###for fullseq
            if(tm_start<=0):
                tm_start=1
            #tm_end = int(row.strip().split(",")[4])+(int(row.strip().split(",")[3])-int(row.strip().split(",")[2])+1)
            #tm_end = int(row.strip().split(",")[4])+(int(row.strip().split(",")[3])-int(row.strip().split(",")[2])+1)+5 ##for surr20
            tm_end = int(row.strip().split(",")[3])  + 5  ###for fullseq
            if tm_end>int(row.strip().split(",")[1]):
                tm_end=int(row.strip().split(",")[1]) # quals to the full sequence length
        #xml_file = os.path.join(set_["xml_file_folder"], database, "%s.xml") % acc
        blast_xml_file = os.path.join(set_["xml_file_folder"], database, "{}.surr{}.xml".format(acc,set_["num_of_sur_residues"]))
        match_details_dict = {}

        if os.path.isfile(blast_xml_file):
            #homo_out_csv_file=os.path.join(set_["homologous_folder"],"NoRedundPro/%s_homo.csv") %acc
            homo_out_csv_file = os.path.join(set_["homologous_folder"], "a3m",database,
                                             "%s_homo%s.csv") % (acc,set_["surres"])
            with open(homo_out_csv_file,'w') as homo_out_csv_file_handle:
                xml_result_handle=open(blast_xml_file)
                xml_record=NCBIXML.read(xml_result_handle)
                E_VALUE_THRESH=set_["e_value_cutoff"]
                hit_num=0
                for alignment in xml_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect<=E_VALUE_THRESH:  #set homologous evalue cutoff
                            match_details_dict['hit_num']=hit_num
                            query_seq_no_gap=re.sub('-','',hsp.query)
                            if hsp.query_start <= tm_start and hsp.query_end >= tm_end:
                                tm_str_start = tm_start - hsp.query_start
                                tm_str_end = tm_end - hsp.query_start + 1
                                k = 0
                                j = 0
                                tm_query_str = ''
                                tm_sbjt_str = ''
                                tm_match_str = ""
                                for char in hsp.query:
                                    if char != '-':
                                        if j >= tm_str_start and j < tm_str_end:
                                            tm_query_str += query_seq_no_gap[j]
                                            tm_sbjt_str += hsp.sbjct[k]
                                            tm_match_str += hsp.match[k]
                                        j = j + 1
                                    k = k + 1
                                match_details_dict["tm_query_seq"]=tm_query_str
                                match_details_dict["tm_match_seq"] = tm_match_str
                                match_details_dict["tm_sbjt_seq"] = tm_sbjt_str
                            if(hit_num)==0:
                                description="%s_NCBI_query_sequence" % acc
                            else:
                                description=alignment.title
                            match_details_dict["description"]=description
                            taxonomy = re.search('\[(.*?)\]',alignment.title)
                            if taxonomy:
                                taxonomyNode=taxonomy.group(1)
                                match_details_dict["organism"]=taxonomyNode
                            #sequence has no organism in the database
                            match_details_dict["organism"]="no_organism"
                            #e_value for hit
                            match_details_dict["FASTA_expectation"]=hsp.expect
                            #convert identity from e.g. 80 (80%) to 0.8
                            match_details_dict["FASTA_identity"]=hsp.identities/100
                            match_details_dict["query_align_seq"]=hsp.query
                            match_details_dict["subject_align_seq"] = hsp.sbjct
                            match_details_dict["match_markup_seq"] = hsp.match
                            match_details_dict["query_start"] = hsp.query_start
                            match_details_dict["query_end"] = hsp.query_end
                            match_details_dict["subject_start"] = hsp.sbjct_start
                            match_details_dict["subject_end"] = hsp.sbjct_end

                            #write the header to the header of the csv file
                            if hit_num==0:
                                csv_header_for_ncbi_homologous_file=sorted(list(match_details_dict.keys()))
                                writer=csv.writer(homo_out_csv_file_handle, delimiter=',',quotechar='"',lineterminator='\n',quoting=csv.QUOTE_NONNUMERIC,doublequote=True)
                                writer.writerow(csv_header_for_ncbi_homologous_file)
                            #save the math_details_dict into the csv file
                            writer=csv.DictWriter(homo_out_csv_file_handle, fieldnames=csv_header_for_ncbi_homologous_file,extrasaction='ignore',delimiter=',',quotechar='"',lineterminator='\n',quoting=csv.QUOTE_MINIMAL,doublequote=True)
                            writer.writerow(match_details_dict)
                            hit_num+=1
                xml_result_handle.close()
            homo_out_csv_file_handle.close()
            logging.info("{} parse_NCBI_xml_to_csv finished ({})".format(acc, homo_out_csv_file))
    tmp_file_handle.close()
    logging.info('~~~~~~~~~~~~                 finished parse_NCBI_xml_to_csv              ~~~~~~~~~~~~')


def extract_filtered_csv_homologous_to_alignments(set_,logging):

    logging.info('start extract filtered csv homologues to alignments')
    min_identity_of_TMD_seq=set_["min_identity_of_TMD_seq"]
    max_identity_of_TMD_seq=set_["max_identity_of_TMD_seq"]
    max_n_gaps_in_TMD_seq=set_["max_n_gaps_in_TMD_seq"]
    max_n_gaps_in_TMD_seq = 1
    tmp_list_loc = set_["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for line in tmp_file_handle:
        tm_protein = line.strip().split(",")[0]
        database = line.strip().split(",")[6]
        #if tm_protein == "2j58_C":
        alignment_dict={}
        alignment_dict1 = {}
        alignment_dict2 = {}
        #homo_csv_file_loc = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo.csv") % tm_protein
        homo_csv_file_loc = os.path.join(set_["homologous_folder"],"a3m",database, "%s_homo%s.csv") % (tm_protein,set_["surres"])
        if (os.path.isfile(homo_csv_file_loc) and os.path.getsize(homo_csv_file_loc) > 0):  # whether protein homologous csv file exists
            #homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo_filtered_aln.fasta") % tm_protein
            homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"],"a3m",database,
                                                      "%s_homo_filtered_aln%s.fasta") % (tm_protein,set_["surres"])
            homo_filtered_out_csv_file_handle=open(homo_filtered_out_csv_file, 'w')
            #homo_filter_file = os.path.join(set_["homologous_folder"], "a3m/NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
            homo_filter_file = os.path.join(set_["homologous_folder"],"a3m",database, "%s.a3m.mem.uniq.2gaps%s") % (tm_protein,set_["surres"])
            homo_filter_file_handle = open(homo_filter_file, "w")
            #homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m/NoRedundPro/%s.mem.lips.input") % tm_protein
            homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m",database,"%s.mem.lips.input%s") % (tm_protein,set_["surres"])
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            #with open(homo_filtered_out_csv_file, 'w') as homo_filtered_out_csv_file_handle:
            df=pd.read_csv(homo_csv_file_loc)
            i=0
            try:
                """WARNING : CURRENTLY Q12983 (BNIP3) for ETRA gives an error here
                 - need to re-organise script so the sequences are simply saved based on a filtered dataframe
                
                KeyError: 'tm_query_seq'
                Exception ignored in: <_io.FileIO name='D:\\\\data_thoipapy\\\\Input_data\\set12_BNIP3_ETRA_processed.csv' mode='rb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\\\data_thoipapy\\\\Input_data\\set12_BNIP3_ETRA_processed.csv' mode='r' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983_homo_filtered_aln_surr0.fasta' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983_homo_filtered_aln_surr0.fasta' mode='w' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983.a3m.mem.uniq.2gaps_surr0' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983.a3m.mem.uniq.2gaps_surr0' mode='w' encoding='cp1252'>
                Exception ignored in: <_io.FileIO name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983.mem.lips.input_surr0' mode='wb' closefd=True>
                ResourceWarning: unclosed file <_io.TextIOWrapper name='D:\\data_thoipapy\\homologous\\a3m\\ETRA\\Q12983.mem.lips.input_surr0' mode='w' encoding='cp1252'>
                """

                for index, row in df.iterrows():
                    if i==0:
                        homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                        homo_filtered_out_csv_file_handle.write(row["tm_query_seq"] + "\n")
                        print("{}".format(row["tm_query_seq"]), file=homo_mem_lips_input_file_handle)
                        print("{}".format(row["tm_query_seq"]), file=homo_filter_file_handle)
                        alignment_dict[row["tm_query_seq"]] = 1
                        i=1
                        continue
                    tm_query_seq_len = len(row["tm_query_seq"])
                    tm_sbjt_seq_gap_num = row["tm_sbjt_seq"].count('-')
                    query_sbjt_identity_num = [x == y for (x, y) in zip(row["tm_query_seq"], row["tm_sbjt_seq"])].count(True)
                    mean_hydrophobicity=thoipapy.common.calc_lipophilicity(row["tm_sbjt_seq"])
                    ratio = SequenceMatcher(None, row["tm_query_seq"], row["tm_sbjt_seq"]).ratio()
                    gap_num = row["tm_sbjt_seq"].count("-")
                    if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and ratio < set_["max_identity_of_TMD_seq"] :  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict:
                            alignment_dict[row["tm_sbjt_seq"]]=1
                            homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_filtered_out_csv_file_handle)
                    if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and ratio < set_["max_identity_of_TMD_seq"] :  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict1:
                            alignment_dict1[row["tm_sbjt_seq"]]=1
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_mem_lips_input_file_handle)
                    if (gap_num <= 2 and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and ratio < set_["max_identity_of_TMD_seq"] ):  # gap number le 3 and no X in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict2:
                            alignment_dict2[row["tm_sbjt_seq"]]=1
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_filter_file_handle)
            except:
                logging.warning("\n------------------\n{} parsing error in extract_filtered_csv_homologous_to_alignments\n{} could not be created. Files will be deleted.\n------------------\n".format(tm_protein, homo_filter_file))
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
            logging.info("{} extract_filtered_csv_homologous_to_alignments finished ({})".format(tm_protein, homo_filter_file))
        else:
            print("{} not found".format(homo_csv_file_loc))
    tmp_file_handle.close()

