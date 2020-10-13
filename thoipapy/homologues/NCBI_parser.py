import csv
import os
import re
import sys
import tarfile
import numpy as np
import pandas as pd
from Bio.Blast import NCBIXML
from pathlib import Path
from thoipapy.utils import create_regex_string, delete_BLAST_xml


def parse_NCBI_xml_to_csv_mult_prot(s, df_set, logging):
    """
    extract multiple sequence alignment and useful informations from xml homologues file
    Parameters
    ----------
    s
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
        #if acc == "Q99IB8":
        database = df_set.loc[i, "database"]
        TMD_start = df_set.loc[i, "TMD_start"]
        TMD_end = df_set.loc[i, "TMD_end"]
        seqlen = df_set.loc[i, "seqlen"]
        if s["surres"] == "_surr0":
            pass
        elif s["surres"] == "_surr5":
            # start 5 residues earlier
            TMD_start = TMD_start - 5 ###for fullseq
            if TMD_start <= 0:
                TMD_start = 1
            # end 5 residues later
            TMD_end = TMD_end  + 5  ###for fullseq
            if TMD_end > seqlen:
                TMD_end = seqlen # quals to the full sequence length
        else:
            raise ValueError('s["surres"] does not seem to be correct')

        BLAST_xml_tar = os.path.join(s["thoipapy_data_folder"], "homologues", "xml", database, "{}.surr{}.BLAST.xml.tar.gz".format(acc, s["num_of_sur_residues"]))
        homo_out_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "ncbi", database)
        if not os.path.isdir(homo_out_dir):
            os.makedirs(homo_out_dir)
        BLAST_csv_tar = os.path.join(homo_out_dir, "{}.surr{}.BLAST.csv.tar.gz".format(acc, s["num_of_sur_residues"]))
        parse_NCBI_xml_to_csv(acc, BLAST_xml_tar, BLAST_csv_tar, TMD_start, TMD_end, s["e_value_cutoff"], logging)

    logging.info('~~~~~~~~~~~~                 finished parse_NCBI_xml_to_csv_mult_prot              ~~~~~~~~~~~~')

def parse_NCBI_xml_to_csv(acc, blast_xml_tar, BLAST_csv_tar, TMD_start, TMD_end, e_value_cutoff, logging):
    # remove the final ".tar.gz" to get the xml and csv filename
    BLAST_xml_file = str(blast_xml_tar)[:-7]
    BLAST_csv_file = str(BLAST_csv_tar)[:-7]

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
            hit_num = 0
            n_hsps_excluded_due_to_e_value_cutoff = 0
            for alignment in xml_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= e_value_cutoff:  # set homologues evalue cutoff
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
                        taxonomy = re.search(r'\[(.*?)\]', alignment.title)
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

    if n_hsps_excluded_due_to_e_value_cutoff > 0:
        logging.info("n_hsps_excluded_due_to_e_value_cutoff = {}".format(n_hsps_excluded_due_to_e_value_cutoff))

    # delete the extracted xml file
    delete_BLAST_xml(BLAST_xml_file)

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

def extract_filtered_csv_homologues_to_alignments_mult_prot(s, df_set, logging):

    logging.info('start extract filtered csv homologues to alignments')
    out_dict = {}

    num_of_sur_residues = s["num_of_sur_residues"]
    max_n_gaps_in_TMD_subject_seq = s["max_n_gaps_in_TMD_subject_seq"]

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        acc_db = acc + "-" + database
        query_TMD_seq = df_set.loc[i, "TMD_seq"]
        query_full_seq = df_set.loc[i, "full_seq"]
        TMD_len = df_set.loc[i, "TMD_len"]
        query_TMD_seq_surr5 = df_set.loc[i, "TMD_seq_pl_surr5"]

        homo_out_dir: Path = Path(s["thoipapy_data_folder"]) / "homologues/ncbi/{database}"
        BLAST_csv_tar: Path = homo_out_dir / f"{acc}.surr{num_of_sur_residues}.BLAST.csv.tar.gz"

        alignments_dir: Path = Path(s["thoipapy_data_folder"]) / f"homologues/alignments/{database}"
        if not alignments_dir.is_dir():
            alignments_dir.mkdir(parents=True)

        fasta_all_TMD_seqs: Path = alignments_dir / f"{acc}.surr{num_of_sur_residues}.gaps{max_n_gaps_in_TMD_subject_seq}.redundant.fas"
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT: Path = alignments_dir / f"{acc}.surr{num_of_sur_residues}.gaps{max_n_gaps_in_TMD_subject_seq}.uniq.for_PSSM_FREECONTACT.txt"
        path_uniq_TMD_seqs_no_gaps_for_LIPS: Path = alignments_dir / f"{acc}.surr{num_of_sur_residues}.gaps0.uniq.for_LIPS.txt"
        path_uniq_TMD_seqs_surr5_for_LIPO: Path = alignments_dir / f"{acc}.surr5.gaps{max_n_gaps_in_TMD_subject_seq}.uniq.for_LIPO.txt"

        single_prot_dict = extract_filtered_csv_homologues_to_alignments(s, acc, TMD_len, fasta_all_TMD_seqs, path_uniq_TMD_seqs_for_PSSM_FREECONTACT,
                                                                         path_uniq_TMD_seqs_no_gaps_for_LIPS, path_uniq_TMD_seqs_surr5_for_LIPO, BLAST_csv_tar,
                                                                         query_TMD_seq, query_TMD_seq_surr5, logging)
        out_dict[acc_db] = single_prot_dict

    df_align_results = pd.DataFrame(out_dict).T
    df_align_results.index.name = "acc"
    align_results_csv = os.path.join(s["thoipapy_data_folder"], "results", s["setname"], "{}_alignment_summary.csv".format(s["setname"]))
    df_align_results.to_csv(align_results_csv)

    logging.info('finished extract filtered csv homologues to alignments for {} proteins. Output = {}'.format(df_align_results.shape[0], align_results_csv))


def extract_filtered_csv_homologues_to_alignments(s: dict,
                                                  acc: str,
                                                  TMD_len: int,
                                                  fasta_all_TMD_seqs: Path,
                                                  path_uniq_TMD_seqs_for_PSSM_FREECONTACT: Path,
                                                  path_uniq_TMD_seqs_no_gaps_for_LIPS: Path,
                                                  path_uniq_TMD_seqs_surr5_for_LIPO: Path,
                                                  BLAST_csv_tar: Path,
                                                  query_TMD_seq: str,
                                                  query_TMD_seq_surr5: str,
                                                  logging):
    fasta_uniq_TMD_seqs_for_PSSM_FREECONTACT: Path = Path(str(path_uniq_TMD_seqs_for_PSSM_FREECONTACT)[:-4] + ".fas")
    fasta_uniq_TMD_seqs_no_gaps_for_LIPS: Path = Path(str(path_uniq_TMD_seqs_no_gaps_for_LIPS)[:-4] + ".fas")
    fasta_uniq_TMD_seqs_surr5_for_LIPO: Path = Path(str(path_uniq_TMD_seqs_surr5_for_LIPO)[:-4] + ".fas")
    alignment_summary_csv: Path = Path(str(fasta_all_TMD_seqs)[:-13] + "alignment_summary.csv")

    single_prot_dict = {}

    # remove the final ".tar.gz" to get the csv filename
    BLAST_csv_file = Path(str(BLAST_csv_tar)[:-7])

    if BLAST_csv_tar.is_file:
        with tarfile.open(BLAST_csv_tar, 'r:gz') as tar:
            BLAST_csv_file_basename = os.path.basename(BLAST_csv_file)
            with tar.extractfile(BLAST_csv_file_basename) as BLAST_csv_extracted:
                df = pd.read_csv(BLAST_csv_extracted)
                n_total_BLAST_hits = df.shape[0]

                # create regex search string, e.g. V-*L-*L-*G-*A-*V-*G-*G-*A-*G-*A-*T-*A-*L-*V-*F-*L-*S-*F-*C
                # (finds sequence regardless of gaps)
                TMD_regex_ss = create_regex_string(query_TMD_seq)

                # obtain the bool, start, end of TMD seqs in the match sequences. Add to the new TMD-specific dataframe.
                start_end_list_in_alignment_ser = df.query_align_seq.apply(get_start_and_end_of_TMD_in_query, args=(TMD_regex_ss,)).dropna()
                df['start'] = start_end_list_in_alignment_ser.apply(lambda x: x[1])
                df['end'] = start_end_list_in_alignment_ser.apply(lambda x: x[2])
                # slice TMDs out of dfs_sel, and save them in the new df_TMD
                df.dropna(subset=["start", "end"], inplace=True)
                df["start"] = df["start"].astype(int)
                df["end"] = df["end"].astype(int)

                if df.empty:
                    BLAST_xml_tar = "{}.surr{}.BLAST.xml.tar.gz".format(acc, s["num_of_sur_residues"])
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
                         "frac_ident_TMD > {}".format(int(s["max_n_gaps_in_TMD_query_seq"]), int(s["max_n_gaps_in_TMD_subject_seq"]),
                                                      float(s["min_identity_of_TMD_seq"])), inplace=True)

                # save all TMD sequences (NON-UNIQUE) for manual inspection of alignments only
                n_total_filtered_seqs = df.shape[0]

                save_fasta(df, "subject_TMD_align_seq", fasta_all_TMD_seqs, acc, query_TMD_seq)

                # save unique sequences WITH gaps (FOR COEVOLUTION WITH FREECONTACT, ETC)
                uniq_TMD_seqs_for_PSSM_FREECONTACT = df.subject_TMD_align_seq.unique()
                # remove seqs with O, U, or J that are not accepted by rate4site
                uniq_TMD_seqs_for_PSSM_FREECONTACT = [x for x in uniq_TMD_seqs_for_PSSM_FREECONTACT if not contains_unaccepted_letter(x)]

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
                uniq_TMD_seqs_surr5_for_LIPO = [x for x in uniq_TMD_seqs_surr5_for_LIPO if not contains_unaccepted_letter(x)]
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
        sys.stdout.write("{} not found".format(BLAST_csv_tar))

    return single_prot_dict

def contains_unaccepted_letter(seq):
    unaccepted_letters = ['O', 'U', 'J']
    return any([u in seq for u in unaccepted_letters])