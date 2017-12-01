import csv
import os
from Bio.Blast import NCBIXML
import re
import pandas as pd
from difflib import SequenceMatcher
import thoipapy
import numpy as np
import sys
import tarfile

def parse_NCBI_xml_to_csv_mult_prot(set_, df_set, logging):
    """
    extract multiple sequence alignment and useful informations from xml homologues file
    Parameters
    ----------
    set_
    logging

    Returns
    -------
    csv file of query protein alignment information

    """
    logging.info('~~~~~~~~~~~~                 starting parse_NCBI_xml_to_csv_mult_prot              ~~~~~~~~~~~~')

    ##############################################################################################
    #                                                                                            #
    #     parse multiple csv files simultaneously                                                #
    #                                                                                            #
    ##############################################################################################

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_start = df_set.loc[i, "TMD_start"]
        TMD_end = df_set.loc[i, "TMD_end"]
        seqlen = df_set.loc[i, "seqlen"]
        if set_["surres"] == "_surr0":
            pass
        elif set_["surres"] == "_surr5":
            # start 5 residues earlier
            TMD_start = TMD_start - 5 ###for fullseq
            if TMD_start <= 0:
                TMD_start = 1
            # end 5 residues later
            TMD_end = TMD_end  + 5  ###for fullseq
            if TMD_end > seqlen:
                TMD_end = seqlen # quals to the full sequence length
        else:
            raise ValueError('set_["surres"] does not seem to be correct')

        BLAST_xml_tar = os.path.join(set_["xml_file_folder"], database, "{}.surr{}.BLAST.xml.tar.gz".format(acc, set_["num_of_sur_residues"]))
        homo_out_dir = os.path.join(set_["homologues_folder"], "ncbi", database)
        if not os.path.isdir(homo_out_dir):
            os.makedirs(homo_out_dir)
        BLAST_csv_tar = os.path.join(homo_out_dir, "{}.surr{}.BLAST.csv.tar.gz".format(acc, set_["num_of_sur_residues"]))

        parse_NCBI_xml_to_csv(set_, acc, BLAST_xml_tar, BLAST_csv_tar, TMD_start, TMD_end, logging)

    logging.info('~~~~~~~~~~~~                 finished parse_NCBI_xml_to_csv_mult_prot              ~~~~~~~~~~~~')

def parse_NCBI_xml_to_csv(set_, acc, blast_xml_tar, BLAST_csv_tar, TMD_start, TMD_end, logging):
    # remove the final ".tar.gz" to get the xml and csv filename
    BLAST_xml_file = blast_xml_tar[:-7]
    BLAST_csv_file = BLAST_csv_tar[:-7]

    match_details_dict = {}

    if not os.path.isfile(blast_xml_tar):
        warning = "{} parse_NCBI_xml_to_csv_mult_prot failed, blast_xml_tar not found = {}".format(acc, blast_xml_tar)
        logging.warning(warning)
        return acc, False, warning

    tar_size = os.path.getsize(blast_xml_tar)
    if tar_size < 100:
        warning = "{} parse_NCBI_xml_to_csv_mult_prot failed, blast_xml_tar seems to be empty, and will be removed".format(acc)
        logging.warning(warning)
        os.remove(blast_xml_tar)
        return acc, False, warning

    # unpack all files in the tarball to the same folder
    # opening as a handle currently doesn't work, as not recognised by NCBIXML.read
    with tarfile.open(blast_xml_tar, 'r:gz') as tar:
        tar.extractall(os.path.dirname(blast_xml_tar))

    with open(BLAST_csv_file, 'w') as homo_out_csv_file_handle:
        with open(BLAST_xml_file) as xml_result_handle:
            xml_record = NCBIXML.read(xml_result_handle)
            E_VALUE_THRESH = set_["e_value_cutoff"]
            hit_num = 0
            n_hsps_excluded_due_to_e_value_cutoff = 0
            for alignment in xml_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= E_VALUE_THRESH:  # set homologues evalue cutoff
                        match_details_dict['hit_num'] = hit_num
                        query_seq_no_gap = re.sub('-', '', hsp.query)
                        if hsp.query_start <= TMD_start and hsp.query_end >= TMD_end:
                            tm_str_start = TMD_start - hsp.query_start
                            tm_str_end = TMD_end - hsp.query_start + 1
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
                            description = "%s_NCBI_query_sequence" % acc
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
                            csv_header_for_ncbi_homologues_file = sorted(list(match_details_dict.keys()))
                            writer = csv.writer(homo_out_csv_file_handle, delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                            writer.writerow(csv_header_for_ncbi_homologues_file)
                        # save the math_details_dict into the csv file
                        writer = csv.DictWriter(homo_out_csv_file_handle, fieldnames=csv_header_for_ncbi_homologues_file, extrasaction='ignore', delimiter=',', quotechar='"', lineterminator='\n',
                                                quoting=csv.QUOTE_MINIMAL, doublequote=True)
                        writer.writerow(match_details_dict)
                        hit_num += 1
                    else:
                        n_hsps_excluded_due_to_e_value_cutoff += 1
                        #sys.stdout.write("|")

    logging.info("n_hsps_excluded_due_to_e_value_cutoff = {}".format(n_hsps_excluded_due_to_e_value_cutoff))

    # delete the extracted xml file
    thoipapy.utils.delete_BLAST_xml(BLAST_xml_file)

    with tarfile.open(BLAST_csv_tar, mode='w:gz') as tar:
        # add the files to the compressed tarfile
        tar.add(BLAST_csv_file, arcname=os.path.basename(BLAST_csv_file))

    # delete the original homologue csv files
    try:
        os.remove(BLAST_csv_file)
    except:
        logging.warning("{} could not be deleted".format(BLAST_csv_file))
    logging.info("{} parse_NCBI_xml_to_csv finished ({})".format(acc, BLAST_csv_tar))

    return acc, True, "no errors"




def get_start_and_end_of_TMD_in_query(x, TMD_regex_ss):
    '''
    define function to obtain regex output (start, stop, etc) as a tuple
    function taken from the korbinian module of Mark Teese
    '''
    m = re.search(TMD_regex_ss, x)
    if m:
        # if the tmd is in the query, return True, start, stop
        return [bool(m), m.start(), m.end()]
    else:
        # if the tmd is not in the query, return False, 0, 0
        return np.nan

def slice_query_TMD_seq(x):
    return x['query_align_seq'][int(x["start"]):int(x["end"])]
def slice_markup_TMD_seq(x):
    return x['match_markup_seq'][int(x["start"]):int(x["end"])]
def slice_match_TMD_seq(x):
    return x['subject_align_seq'][int(x["start"]):int(x["end"])]
def slice_match_TMD_seq_surr5(x):
    return x['subject_align_seq'][int(x["start_min_5"]) : int(x["end"]) + 5]

def save_fasta(df, col_with_seqs, filepath, acc, query_TMD_seq):
    with open(filepath, "w") as f:
        # write orig seq
        f.write(">{}_orig_seq\n{}\n".format(acc, query_TMD_seq))
        for n in df.index:
            seq = df.loc[n, col_with_seqs]
            if seq == query_TMD_seq:
                # do not write the orig seq again. continue to next seq.
                continue
            description = df.loc[n, "description"]
            f.write(">{}\n{}\n".format(description, seq))

def save_fasta_from_array(array_of_seqs, filepath, acc, query_TMD_seq):
    with open(filepath, "w") as f:
        if query_TMD_seq is not None:
            f.write(">{}_orig_seq\n{}\n".format(acc, query_TMD_seq))
            for n, seq in enumerate(array_of_seqs):
                if seq == query_TMD_seq:
                    # do not write the orig seq again. continue to next seq.
                    continue
                f.write(">{}\n{}\n".format(n, seq))
        else:
            for n, seq in enumerate(array_of_seqs):
                f.write(">{}\n{}\n".format(n, seq))

def save_seqs(array_of_seqs, filepath, query_TMD_seq):
    with open(filepath, "w") as f:
        if query_TMD_seq is not None:
            f.write("{}\n".format(query_TMD_seq))
            for seq in array_of_seqs:
                if seq == query_TMD_seq:
                    # do not write the orig seq again. continue to next seq.
                    continue
                f.write("{}\n".format(seq))
        else:
            for seq in array_of_seqs:
                f.write("{}\n".format(seq))

def extract_filtered_csv_homologues_to_alignments_mult_prot(set_, df_set, logging):

    logging.info('start extract filtered csv homologues to alignments')
    out_dict = {}

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        query_TMD_seq = df_set.loc[i, "TMD_seq"]
        query_full_seq = df_set.loc[i, "full_seq"]
        TMD_len = df_set.loc[i, "TMD_len"]
        query_TMD_seq_surr5 = df_set.loc[i, "TMD_seq_pl_surr5"]

        homo_out_dir = os.path.join(set_["homologues_folder"], "ncbi", database)
        BLAST_csv_tar = os.path.join(homo_out_dir, "{}.surr{}.BLAST.csv.tar.gz".format(acc, set_["num_of_sur_residues"]))

        alignments_dir = os.path.join(set_["homologues_folder"], "alignments", database)
        if not os.path.isdir(alignments_dir):
            os.makedirs(alignments_dir)

        fasta_all_TMD_seqs = os.path.join(alignments_dir, "{}.surr{}.gaps{}.redundant.fas".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, set_["num_of_sur_residues"], set_["max_n_gaps_in_TMD_subject_seq"]))
        path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, set_["num_of_sur_residues"]))
        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.txt".format(acc, set_["max_n_gaps_in_TMD_subject_seq"]))


        single_prot_dict = extract_filtered_csv_homologues_to_alignments(set_, acc, TMD_len, fasta_all_TMD_seqs, path_uniq_TMD_seqs_for_PSSM_FREECONTACT,
                                                                         path_uniq_TMD_seqs_no_gaps_for_LIPS, path_uniq_TMD_seqs_surr5_for_LIPO, BLAST_csv_tar,
                                                                         query_TMD_seq, query_TMD_seq_surr5, logging)
        out_dict[acc] = single_prot_dict

    df_align_results = pd.DataFrame(out_dict).T
    df_align_results.index.name = "acc"
    align_results_csv = os.path.join(set_["set_results_folder"], "{}_alignment_summary.csv".format(set_["setname"]))
    df_align_results.to_csv(align_results_csv)

    logging.info('finished extract filtered csv homologues to alignments for {} proteins. Output = {}'.format(df_align_results.shape[0], align_results_csv))


def extract_filtered_csv_homologues_to_alignments(set_, acc, TMD_len, fasta_all_TMD_seqs, path_uniq_TMD_seqs_for_PSSM_FREECONTACT,
                                                  path_uniq_TMD_seqs_no_gaps_for_LIPS, path_uniq_TMD_seqs_surr5_for_LIPO, BLAST_csv_tar,
                                                  query_TMD_seq, query_TMD_seq_surr5, logging):
    fasta_uniq_TMD_seqs_for_PSSM_FREECONTACT = path_uniq_TMD_seqs_for_PSSM_FREECONTACT[:-4] + ".fas"
    fasta_uniq_TMD_seqs_no_gaps_for_LIPS = path_uniq_TMD_seqs_no_gaps_for_LIPS[:-4] + ".fas"
    fasta_uniq_TMD_seqs_surr5_for_LIPO = path_uniq_TMD_seqs_surr5_for_LIPO[:-4] + ".fas"
    alignment_summary_csv = fasta_all_TMD_seqs[:-13] + "alignment_summary.csv"

    single_prot_dict = {}

    # remove the final ".tar.gz" to get the csv filename
    BLAST_csv_file = BLAST_csv_tar[:-7]
    print(BLAST_csv_tar)

    if os.path.isfile(BLAST_csv_tar):
        with tarfile.open(BLAST_csv_tar, 'r:gz') as tar:
            BLAST_csv_file_basename = os.path.basename(BLAST_csv_file)
            with tar.extractfile(BLAST_csv_file_basename) as BLAST_csv_extracted:
                df = pd.read_csv(BLAST_csv_extracted)
                n_total_BLAST_hits = df.shape[0]

                # create regex search string, e.g. V-*L-*L-*G-*A-*V-*G-*G-*A-*G-*A-*T-*A-*L-*V-*F-*L-*S-*F-*C
                # (finds sequence regardless of gaps)
                TMD_regex_ss = thoipapy.utils.create_regex_string(query_TMD_seq)

                # obtain the bool, start, end of TMD seqs in the match sequences. Add to the new TMD-specific dataframe.
                start_end_list_in_alignment_ser = df.query_align_seq.apply(get_start_and_end_of_TMD_in_query, args=(TMD_regex_ss,)).dropna()
                df['start'] = start_end_list_in_alignment_ser.apply(lambda x: x[1])
                df['end'] = start_end_list_in_alignment_ser.apply(lambda x: x[2])
                # slice TMDs out of dfs_sel, and save them in the new df_TMD
                df.dropna(subset=["start", "end"], inplace=True)
                df["start"] = df["start"].astype(int)
                df["end"] = df["end"].astype(int)

                if df.empty:
                    BLAST_xml_tar = "{}.surr{}.BLAST.xml.tar.gz".format(acc, set_["num_of_sur_residues"])
                    raise ValueError("{} extract_filtered_csv_homologues_to_alignments failed\n. None of the homologues contained the original TMD sequence.\n"
                                     "This can occur when the homologue xml file is old, and the sequences from BLAST and the protein set don't match. Delete the"
                                     "homologue xml file ({}), re-run NCBI blast, and try-again".format(acc, BLAST_xml_tar))

                df['query_TMD_align_seq'] = df.apply(slice_query_TMD_seq, axis=1)
                df['markup_TMD_align_seq'] = df.apply(slice_markup_TMD_seq, axis=1)
                df['subject_TMD_align_seq'] = df.apply(slice_match_TMD_seq, axis=1)

                # slice out the TMD plus surrounding 5 residues
                # the start index needs to be adjusted as necessary
                # (as far as I understand, indexing "12345"[0:10] is allowed
                start_min_5 = df["start"] - 5
                start_min_5[start_min_5 < 0] = 0
                df['start_min_5'] = start_min_5
                df['subject_TMD_align_seq_surr5'] = df.apply(slice_match_TMD_seq_surr5, axis=1)

                df['X_in_subject_TMD_align_seq_surr5'] = df['subject_TMD_align_seq_surr5'].str.contains("X")

                # calculate fraction identity of TMD region only
                df["frac_ident_TMD"] = (df.markup_TMD_align_seq.str.len() - df.markup_TMD_align_seq.str.replace("+", " ").str.count(" ")) / TMD_len

                # count gaps
                df["n_gaps_query_align_seq"] = df.query_TMD_align_seq.str.count("-")
                df["n_gaps_subject_align_seq"] = df.subject_TMD_align_seq.str.count("-")

                # filter by gaps in query, gaps in subject, and fraction identity of TMD
                df.query("n_gaps_query_align_seq <= {} & n_gaps_subject_align_seq <= {} & X_in_subject_TMD_align_seq_surr5 == False & "
                         "frac_ident_TMD > {}".format(int(set_["max_n_gaps_in_TMD_query_seq"]), int(set_["max_n_gaps_in_TMD_subject_seq"]),
                                                      float(set_["min_identity_of_TMD_seq"])), inplace=True)

                # save all TMD sequences (NON-UNIQUE) for manual inspection of alignments only
                n_total_filtered_seqs = df.shape[0]

                save_fasta(df, "subject_TMD_align_seq", fasta_all_TMD_seqs, acc, query_TMD_seq)

                # save unique sequences WITH gaps (FOR COEVOLUTION WITH FREECONTACT, ETC)
                uniq_TMD_seqs_for_PSSM_FREECONTACT = df.subject_TMD_align_seq.unique()
                save_seqs(uniq_TMD_seqs_for_PSSM_FREECONTACT, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, query_TMD_seq)
                save_fasta_from_array(uniq_TMD_seqs_for_PSSM_FREECONTACT, fasta_uniq_TMD_seqs_for_PSSM_FREECONTACT, acc, query_TMD_seq)

                # save unique sequences WITHOUT gaps (FOR LIPS)
                uniq_TMD_seqs_no_gaps_for_LIPS = [seq for seq in uniq_TMD_seqs_for_PSSM_FREECONTACT if "-" not in seq]
                save_seqs(uniq_TMD_seqs_no_gaps_for_LIPS, path_uniq_TMD_seqs_no_gaps_for_LIPS, query_TMD_seq)
                save_fasta_from_array(uniq_TMD_seqs_no_gaps_for_LIPS, fasta_uniq_TMD_seqs_no_gaps_for_LIPS, acc, query_TMD_seq)

                # save unique sequences WITH gaps with 5 surrounding residues (FOR PSSM and Hessa LIPO)
                # delete any longer sequences, where the query had gaps
                # assume the first sequence has no gaps
                TMD_plus_5_len = len(df.iloc[0, :]["subject_TMD_align_seq_surr5"])

                # only keep the seqs that have the same length as the first one
                df_no_gaps_in_q_plus5 = df.loc[df['subject_TMD_align_seq_surr5'].str.len() == TMD_plus_5_len]
                uniq_TMD_seqs_surr5_for_LIPO = df_no_gaps_in_q_plus5['subject_TMD_align_seq_surr5'].unique()
                print("query_TMD_seq_surr5", query_TMD_seq_surr5)
                save_seqs(uniq_TMD_seqs_surr5_for_LIPO, path_uniq_TMD_seqs_surr5_for_LIPO, query_TMD_seq=query_TMD_seq_surr5)
                save_fasta_from_array(uniq_TMD_seqs_surr5_for_LIPO, fasta_uniq_TMD_seqs_surr5_for_LIPO, acc, query_TMD_seq=query_TMD_seq_surr5)

                single_prot_dict["n_total_BLAST_hits"] = n_total_BLAST_hits
                single_prot_dict["n_total_filtered_seqs"] = n_total_filtered_seqs
                single_prot_dict["n_uniq_TMD_seqs_for_PSSM_FREECONTACT"] = len(uniq_TMD_seqs_for_PSSM_FREECONTACT)
                single_prot_dict["n_uniq_TMD_seqs_no_gaps_for_LIPS"] = len(uniq_TMD_seqs_no_gaps_for_LIPS)
                single_prot_dict["n_uniq_TMD_seqs_surr5_for_LIPO"] = len(uniq_TMD_seqs_surr5_for_LIPO)

                single_prot_aln_result_ser = pd.Series(single_prot_dict)
                single_prot_aln_result_ser.to_csv(alignment_summary_csv)

                logging.info("{} extract_filtered_csv_homologues_to_alignments finished ({}). {}, {}, and {} valid seqs "
                             "from {} total".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, single_prot_dict["n_uniq_TMD_seqs_for_PSSM_FREECONTACT"], single_prot_dict["n_uniq_TMD_seqs_no_gaps_for_LIPS"],
                                                    single_prot_dict["n_uniq_TMD_seqs_surr5_for_LIPO"], n_total_filtered_seqs))

    else:
        print("{} not found".format(BLAST_csv_tar))

    return single_prot_dict

# DEPRECATED METHOD
def extract_filtered_csv_homologues_to_alignments_orig_handle_method(set_,logging):

    logging.info('start extract filtered csv homologues to alignments')
    min_identity_of_TMD_seq=set_["min_identity_of_TMD_seq"]
    max_identity_of_TMD_seq=set_["max_identity_of_TMD_seq"]
    max_n_gaps_in_TMD_subject_seq=set_["max_n_gaps_in_TMD_subject_seq"]
    tmp_list_loc = set_["list_of_tmd_start_end"]
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
        #homo_csv_file_loc = os.path.join(set_["homologues_folder"], "NoRedundPro/%s_homo.csv") % acc
        homo_csv_file_loc = os.path.join(set_["homologues_folder"],"a3m",database, "%s_homo%s.csv") % (acc,set_["surres"])
        print(homo_csv_file_loc)
        if (os.path.isfile(homo_csv_file_loc) and os.path.getsize(homo_csv_file_loc) > 0):  # whether protein homologues csv file exists
            df = pd.read_csv(homo_csv_file_loc)
            print(df["tm_sbjt_seq"])
            #homo_filtered_out_csv_file = os.path.join(set_["homologues_folder"], "NoRedundPro/%s_homo_filtered_aln.fasta") % acc
            homo_filtered_out_csv_file = os.path.join(set_["homologues_folder"],"a3m",database,
                                                      "%s_homo_filtered_aln%s.fasta") % (acc,set_["surres"])
            homo_filtered_out_csv_file_handle=open(homo_filtered_out_csv_file, 'w')
            #homo_filter_file = os.path.join(set_["homologues_folder"], "a3m/NoRedundPro/%s.a3m.mem.uniq.2gaps") % acc
            homo_filter_file = os.path.join(set_["homologues_folder"],"a3m",database, "%s.a3m.mem.uniq.2gaps%s") % (acc,set_["surres"])
            homo_filter_file_handle = open(homo_filter_file, "w")
            #homo_mem_lips_input_file = os.path.join(set_["homologues_folder"],"a3m/NoRedundPro/%s.mem.lips.input") % acc
            homo_mem_lips_input_file = os.path.join(set_["homologues_folder"],"a3m",database,"%s.mem.lips.input%s") % (acc,set_["surres"])
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
                        print("{}".format(row["tm_query_seq"]), file=homo_mem_lips_input_file_handle)
                        print("{}".format(row["tm_query_seq"]), file=homo_filter_file_handle)
                        alignment_dict[row["tm_query_seq"]] = 1
                        i=1
                        continue
                    tm_query_seq_len = len(row["tm_query_seq"])
                    tm_sbjt_seq_gap_num = row["tm_sbjt_seq"].count('-')
                    print('row["tm_sbjt_seq"]', row["tm_sbjt_seq"])
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
            print("{} not found".format(homo_csv_file_loc))
    tmp_file_handle.close()

