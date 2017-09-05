import csv
import os
from Bio.Blast import NCBIXML
import re
import pandas as pd
import pickle
import tarfile
import xml.etree.ElementTree as ET
import thoipapy.mtutiles as utils
import zipfile
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from difflib import SequenceMatcher

def parse_NCBI_xml_to_csv(pathdict,test_protein, set_,logging):
    counter_xml_to_csv=0
    logging.info('~~~~~~~~~~starting parse_NCBI_xml_to_csv')


    ##############################################################################################
    #                                                                                            #
    #     parse multiple csv files simultaneously                                                #
    #                                                                                            #
    ##############################################################################################


    if (set_["multiple_tmp_simultaneous"]):
        tmp_list_loc=set_["list_of_tmd_start_end"]
        tmp_file_handle=open(tmp_list_loc,'r')
        match_details_dict = {}
        for row in tmp_file_handle:
            tmp_protein_acc=row.strip().split("\t")[0]
            #differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)
            if len(tmp_protein_acc) == 6:    #uniprot acc
                tmp_protein_name=tmp_protein_acc
                tm_start = int(row.strip().split("\t")[3])
                tm_end = int(row.strip().split("\t")[4])
                #xml_file = os.path.join(set_["xml_file_folder"], "%s.xml") % tmp_protein_name
                xml_file = os.path.join(set_["xml_file_folder"], "%s/%s.xml") %(set_["Datatype"],tmp_protein_name)
            else:                  #PDB proteins
                tmp_protein_name=tmp_protein_acc[0:7]
                tm_start = int(row.strip().split("\t")[3])
                tm_end = int(row.strip().split("\t")[4])
                #xml_file = os.path.join(set_["xml_file_folder"], "%s.xml") % tmp_protein_name[0:6]
                xml_file = os.path.join(set_["xml_file_folder"], "%s/%s.xml") %(set_["Datatype"],tmp_protein_name[0:6])
            #if tmp_protein_name == test_protein:
            #xml_file = os.path.join(set_["xml_file_folder"], "%s.xml") % tmp_protein_name

            if os.path.isfile(xml_file):
                #homo_out_csv_file=os.path.join(set_["homologous_folder"],"NoRedundPro/%s_homo.csv") %tmp_protein_name
                homo_out_csv_file = os.path.join(set_["homologous_folder"],
                                                 "%s/%s_homo.csv") %(set_["Datatype"], tmp_protein_name)
                with open(homo_out_csv_file,'w') as homo_out_csv_file_handle:
                    xml_result_handle=open(xml_file)
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
                                    description="%s_NCBI_query_sequence" %tmp_protein_name
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


    ##############################################################################################
    #                                                                                            #
    #     parse single tmp csv file                                                              #
    #                                                                                            #
    ##############################################################################################


    else:
            tmp_protein_name = set_["tm_protein_name"]
            xml_file = os.path.join(set_["xml_file_folder"], "%s.xml") % tmp_protein_name
            match_details_dict = {}
            tm_start=int(set_["tm_start"])
            tm_end=int(set_["tm_end"])
            if os.path.isfile(xml_file):
                homo_out_csv_file = os.path.join(set_["homologous_folder"], "%s/%s_homo.csv") %(set_["Datatype"],tmp_protein_name)
                with open(homo_out_csv_file, 'w') as homo_out_csv_file_handle:
                    xml_result_handle = open(xml_file)
                    xml_record = NCBIXML.read(xml_result_handle)
                    E_VALUE_THRESH = set_["e_value_cutoff"]
                    hit_num = 0
                    for alignment in xml_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect <= E_VALUE_THRESH:  # set homologous evalue cutoff
                                match_details_dict['hit_num'] = hit_num
                                query_seq_no_gap = re.sub('-', '', hsp.query)
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
                                    match_details_dict["tm_query_seq"] = tm_query_str
                                    match_details_dict["tm_match_seq"] = tm_match_str
                                    match_details_dict["tm_sbjt_seq"] = tm_sbjt_str
                                if (hit_num) == 0:
                                    description = "%s_NCBI_query_sequence" % tmp_protein_name
                                else:
                                    description = alignment.title
                                match_details_dict["description"] = description
                                taxonomy = re.search('\[(.*?)\]', alignment.title)
                                if taxonomy:
                                    taxonomyNode = taxonomy.group(1)
                                    match_details_dict["organism"] = taxonomyNode
                                # sequence has no organism in the database
                                match_details_dict["organism"] = "no_organism"
                                # e_value for hit
                                match_details_dict["FASTA_expectation"] = hsp.expect
                                # convert identity from e.g. 80 (80%) to 0.8
                                match_details_dict["FASTA_identity"] = hsp.identities / 100
                                match_details_dict["query_align_seq"] = hsp.query
                                match_details_dict["subject_align_seq"] = hsp.sbjct
                                match_details_dict["match_markup_seq"] = hsp.match
                                match_details_dict["query_start"] = hsp.query_start
                                match_details_dict["query_end"] = hsp.query_end
                                match_details_dict["subject_start"] = hsp.sbjct_start
                                match_details_dict["subject_end"] = hsp.sbjct_end

                                # write the header to the header of the csv file
                                if hit_num == 0:
                                    csv_header_for_ncbi_homologous_file = sorted(list(match_details_dict.keys()))
                                    writer = csv.writer(homo_out_csv_file_handle, delimiter=',', quotechar='"',
                                                        lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC,
                                                        doublequote=True)
                                    writer.writerow(csv_header_for_ncbi_homologous_file)
                                # save the math_details_dict into the csv file
                                writer = csv.DictWriter(homo_out_csv_file_handle,
                                                        fieldnames=csv_header_for_ncbi_homologous_file,
                                                        extrasaction='ignore', delimiter=',', quotechar='"',
                                                        lineterminator='\n', quoting=csv.QUOTE_MINIMAL,
                                                        doublequote=True)
                                writer.writerow(match_details_dict)
                                hit_num += 1
                xml_result_handle.close()
            


    # ##############################################################################################
    # #                                                                                            #
    # #     create csv file for 3 xml files of mark                                                #
    # #                                                                                            #
    # ##############################################################################################
    # else:
    #         tmp_protein_name = set_["tm_protein_name"]
    #         xml_file = os.path.join(set_["xml_file_folder"], "mark/%s.xml") % tmp_protein_name
    #         match_details_dict = {}
    #         if os.path.isfile(xml_file):
    #             homo_out_csv_file = os.path.join(set_["homologous_folder"], "%s_homo.csv") % tmp_protein_name
    #             with open(homo_out_csv_file, 'w') as homo_out_csv_file_handle:
    #                 xml_result_handle = open(xml_file)
    #                 xml_record = NCBIXML.read(xml_result_handle)
    #                 E_VALUE_THRESH = set_["e_value_cutoff"]
    #                 hit_num = 0
    #                 for alignment in xml_record.alignments:
    #                     for hsp in alignment.hsps:
    #                         if hsp.expect <= E_VALUE_THRESH:  # set homologous evalue cutoff
    #                             match_details_dict['hit_num'] = hit_num
    #                             # query_seq_no_gap = re.sub('-', '', hsp.query)
    #                             # if hsp.query_start <= tm_start and hsp.query_end >= tm_end:
    #                             #     tm_str_start = tm_start - hsp.query_start
    #                             #     tm_str_end = tm_end - hsp.query_start + 1
    #                             #     k = 0
    #                             #     j = 0
    #                             #     tm_query_str = ''
    #                             #     tm_sbjt_str = ''
    #                             #     tm_match_str = ""
    #                             #     for char in hsp.query:
    #                             #         if char != '-':
    #                             #             if j >= tm_str_start and j < tm_str_end:
    #                             #                 tm_query_str += query_seq_no_gap[j]
    #                             #                 tm_sbjt_str += hsp.sbjct[k]
    #                             #                 tm_match_str += hsp.match[k]
    #                             #             j = j + 1
    #                             #         k = k + 1
    #                             #     match_details_dict["tm_query_seq"] = tm_query_str
    #                             #     match_details_dict["tm_match_seq"] = tm_match_str
    #                             #     match_details_dict["tm_sbjt_seq"] = tm_sbjt_str
    #                             if (hit_num) == 0:
    #                                 description = "%s_NCBI_query_sequence" % tmp_protein_name
    #                             else:
    #                                 description = alignment.title
    #                             match_details_dict["description"] = description
    #                             taxonomy = re.search('\[(.*?)\]', alignment.title)
    #                             if taxonomy:
    #                                 taxonomyNode = taxonomy.group(1)
    #                                 match_details_dict["organism"] = taxonomyNode
    #                             # sequence has no organism in the database
    #                             match_details_dict["organism"] = "no_organism"
    #                             # e_value for hit
    #                             match_details_dict["FASTA_expectation"] = hsp.expect
    #                             # convert identity from e.g. 80 (80%) to 0.8
    #                             match_details_dict["FASTA_identity"] = hsp.identities / 100
    #                             match_details_dict["query_align_seq"] = hsp.query
    #                             match_details_dict["subject_align_seq"] = hsp.sbjct
    #                             match_details_dict["match_markup_seq"] = hsp.match
    #                             match_details_dict["query_start"] = hsp.query_start
    #                             match_details_dict["query_end"] = hsp.query_end
    #                             match_details_dict["subject_start"] = hsp.sbjct_start
    #                             match_details_dict["subject_end"] = hsp.sbjct_end
    #
    #                             # write the header to the header of the csv file
    #                             if hit_num == 0:
    #                                 csv_header_for_ncbi_homologous_file = sorted(list(match_details_dict.keys()))
    #                                 writer = csv.writer(homo_out_csv_file_handle, delimiter=',', quotechar='"',
    #                                                     lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC,
    #                                                     doublequote=True)
    #                                 writer.writerow(csv_header_for_ncbi_homologous_file)
    #                             # save the math_details_dict into the csv file
    #                             writer = csv.DictWriter(homo_out_csv_file_handle,
    #                                                     fieldnames=csv_header_for_ncbi_homologous_file,
    #                                                     extrasaction='ignore', delimiter=',', quotechar='"',
    #                                                     lineterminator='\n', quoting=csv.QUOTE_MINIMAL,
    #                                                     doublequote=True)
    #                             writer.writerow(match_details_dict)
    #                             hit_num += 1


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



def extract_filtered_csv_homologous_to_alignments(tmp_lists, test_protein,pathdict, set_,logging):
    if (set_["multiple_tmp_simultaneous"]):
        logging.info('start extract filterd csv homologous')
        tmp_lists=tmp_lists
        min_identity_of_TMD_seq=set_["min_identity_of_TMD_seq"]
        max_identity_of_TMD_seq=set_["max_identity_of_TMD_seq"]
        max_n_gaps_in_TMD_seq=set_["max_n_gaps_in_TMD_seq"]
        max_n_gaps_in_TMD_seq = 1
        for tm_protein in tmp_lists:
            #if tm_protein==test_protein:
            alignment_dict={}
            alignment_dict1 = {}
            alignment_dict2 = {}
            #homo_csv_file_loc = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo.csv") % tm_protein
            homo_csv_file_loc = os.path.join(set_["homologous_folder"], "%s/%s_homo.csv") %(set_["Datatype"], tm_protein)
            if (os.path.isfile(homo_csv_file_loc)):  # whether protein homologous csv file exists
                #homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo_filtered_aln.fasta") % tm_protein
                homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"],
                                                          "%s/%s_homo_filtered_aln.fasta") %(set_["Datatype"], tm_protein)
                homo_filtered_out_csv_file_handle=open(homo_filtered_out_csv_file, 'w')
                #homo_filter_file = os.path.join(set_["homologous_folder"], "a3m/NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
                homo_filter_file = os.path.join(set_["homologous_folder"], "a3m/%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"], tm_protein)
                homo_filter_file_handle = open(homo_filter_file, "w")
                #homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m/NoRedundPro/%s.mem.lips.input") % tm_protein
                homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m/%s/%s.mem.lips.input") %(set_["Datatype"],tm_protein)
                homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
                #with open(homo_filtered_out_csv_file, 'w') as homo_filtered_out_csv_file_handle:
                df=pd.read_csv(homo_csv_file_loc,sep=",",quoting=csv.QUOTE_NONNUMERIC, index_col=0)
                i=0
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
                    mean_hydrophobicity=calc_lipophilicity(row["tm_sbjt_seq"])
                    ratio = SequenceMatcher(None, row["tm_query_seq"], row["tm_sbjt_seq"]).ratio()
                    gap_num = row["tm_sbjt_seq"].count("-")
                    if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict:
                            alignment_dict[row["tm_sbjt_seq"]]=1
                            homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_filtered_out_csv_file_handle)
                    if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict1:
                            alignment_dict1[row["tm_sbjt_seq"]]=1
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_mem_lips_input_file_handle)
                    if (gap_num <= 2 and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                        if not row["tm_sbjt_seq"] in alignment_dict2:
                            alignment_dict2[row["tm_sbjt_seq"]]=1
                            print("{}".format(row["tm_sbjt_seq"]), file=homo_filter_file_handle)
                homo_filter_file_handle.close()
                homo_mem_lips_input_file_handle.close()
                homo_filtered_out_csv_file_handle.close()


    else:
        logging.info('start extract filterd csv homologous')
        tmp_lists=tmp_lists
        min_identity_of_TMD_seq=set_["min_identity_of_TMD_seq"]
        max_identity_of_TMD_seq=set_["max_identity_of_TMD_seq"]
        max_n_gaps_in_TMD_seq=set_["max_n_gaps_in_TMD_seq"]
        max_n_gaps_in_TMD_seq = 1
        #for tm_protein in tmp_lists:
        tm_protein=set_["tm_protein_name"]
        #if tm_protein==test_protein:
        alignment_dict={}
        alignment_dict1 = {}
        alignment_dict2 = {}
        #homo_csv_file_loc = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo.csv") % tm_protein
        homo_csv_file_loc = os.path.join(set_["homologous_folder"], "%s/%s_homo.csv") %(set_["Datatype"], tm_protein)
        if (os.path.isfile(homo_csv_file_loc)):  # whether protein homologous csv file exists
            #homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"], "NoRedundPro/%s_homo_filtered_aln.fasta") % tm_protein
            homo_filtered_out_csv_file = os.path.join(set_["homologous_folder"],
                                                        "%s/%s_homo_filtered_aln.fasta") %(set_["Datatype"], tm_protein)
            homo_filtered_out_csv_file_handle=open(homo_filtered_out_csv_file, 'w')
            #homo_filter_file = os.path.join(set_["homologous_folder"], "a3m/NoRedundPro/%s.a3m.mem.uniq.2gaps") % tm_protein
            homo_filter_file = os.path.join(set_["homologous_folder"], "a3m/%s/%s.a3m.mem.uniq.2gaps") %(set_["Datatype"], tm_protein)
            homo_filter_file_handle = open(homo_filter_file, "w")
            #homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m/NoRedundPro/%s.mem.lips.input") % tm_protein
            homo_mem_lips_input_file = os.path.join(set_["homologous_folder"],"a3m/%s/%s.mem.lips.input") %(set_["Datatype"],tm_protein)
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            #with open(homo_filtered_out_csv_file, 'w') as homo_filtered_out_csv_file_handle:
            df=pd.read_csv(homo_csv_file_loc,sep=",",quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            i=0
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
                mean_hydrophobicity=calc_lipophilicity(row["tm_sbjt_seq"])
                ratio = SequenceMatcher(None, row["tm_query_seq"], row["tm_sbjt_seq"]).ratio()
                gap_num = row["tm_sbjt_seq"].count("-")
                if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    if not row["tm_sbjt_seq"] in alignment_dict:
                        alignment_dict[row["tm_sbjt_seq"]]=1
                        homo_filtered_out_csv_file_handle.write(">" + row["description"] + "\n")
                        print("{}".format(row["tm_sbjt_seq"]), file=homo_filtered_out_csv_file_handle)
                if not re.search("-", row["tm_sbjt_seq"]) and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    if not row["tm_sbjt_seq"] in alignment_dict1:
                        alignment_dict1[row["tm_sbjt_seq"]]=1
                        print("{}".format(row["tm_sbjt_seq"]), file=homo_mem_lips_input_file_handle)
                if (gap_num <= 2 and not re.search("X",row["tm_sbjt_seq"]) and ratio >= set_["min_identity_of_TMD_seq"] and mean_hydrophobicity < set_["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    if not row["tm_sbjt_seq"] in alignment_dict2:
                        alignment_dict2[row["tm_sbjt_seq"]]=1
                        print("{}".format(row["tm_sbjt_seq"]), file=homo_filter_file_handle)
            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()
            homo_filtered_out_csv_file_handle.close()
        

                    # if (gap_num <= 3 and not re.search("X",
                    #                                    row["tm_sbjt_seq"]) and ratio >= set_[
                    #     "min_identity_of_TMD_seq"] and mean_hydrophobicity < set_[
                    #     "max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    #     print("{}".format(row["tm_sbjt_seq"]), file=homo_filter_file_handle)

                    # if tm_query_seq_len * min_identity_of_TMD_seq <= query_sbjt_identity_num < tm_query_seq_len*max_identity_of_TMD_seq and not re.search("X",row["tm_sbjt_seq"]) and mean_hydrophobicity<set_["max_hydrophilicity_Hessa"]:
                    #     #if tm_sbjt_seq_gap_num <= max_n_gaps_in_TMD_seq:
                    #     if tm_sbjt_seq_gap_num < max_n_gaps_in_TMD_seq:
                    #         if not row["tm_sbjt_seq"] in alignment_dict:
                    #             alignment_dict[row["tm_sbjt_seq"]] = 1
                    #             #alignment_dict.pop(row["tm_query_seq"], None)
                    #             homo_filtered_out_csv_file_handle.write(">" + row["description"] +"\n")
                    #             homo_filtered_out_csv_file_handle.write(row["tm_sbjt_seq"] + "\n")
