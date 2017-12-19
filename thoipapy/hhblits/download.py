import sys
import os
import subprocess,threading
import time
import thoipapy.utils as utils


def parse_a3m_alignment(s, logging):
    logging.info('start parsing a3m file')

    tmp_list_loc = s["list_of_tmd_start_end"]
    #if (s["multiple_tmp_homo_download"]):
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0]

        output_oa3m_homologues_file = os.path.join(s["Protein_folder"], "%s.fasta.surr20.a3m") % acc
        output_oa3m_homologues_processed_file = os.path.join(s["output_oa3m_homologues"],"%s.surr20.parse.a3m") % acc
        exect_str1 = "grep -v '^>' {a3m_file} |sed 's/[a-z]//g' >{a3m_processed_file}".format(a3m_file=output_oa3m_homologues_file,
                                                                                              a3m_processed_file=output_oa3m_homologues_processed_file)

        if os.path.isfile(output_oa3m_homologues_file):
            command1 = utils.Command(exect_str1)
            command1.run(timeout=20)
            logging.info("a3m processed Output file: %s\n" % output_oa3m_homologues_processed_file)
    tmp_file_handle.close()





def download_homologues_with_hhblits(s,logging):
    logging.info('start downloading protein homologues with hhblits')
    hhblits_loc = s["hhblits_dir"]
    uniprot_database_loc = s["uniprot_database_dir"]

    tmp_list_loc = s["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    #skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0]
        database = row.strip().split(",")[6]
        tmp_protein_fasta = os.path.join(s["Protein_folder"], "hhblits", database, "%s.surr20.fasta") % acc
        if os.path.isfile(tmp_protein_fasta):
            output_oa3m_homologues_file = os.path.join(s["Protein_folder"], "hhblits", database,  "%s.fasta.surr20.a3m") % acc
            if acc=="QuePro" and os.path.isfile(output_oa3m_homologues_file):
                os.system("rm %s" % output_oa3m_homologues_file)          ###remove QuePro.fasta.a3m that is the previous queried protein

            exect_str = "qsub -b y -q all.q -N 'hhblitsjob' {hhblits} -i {fasta_input} -oa3m {a3m_output} -d {uniprot_database} -Z 999999999 " \
                        "-B 999999999 -maxfilt 999999999 -id 99 -diff inf".format(hhblits=hhblits_loc, fasta_input=tmp_protein_fasta,
                                                                                  a3m_output=output_oa3m_homologues_file,uniprot_database=uniprot_database_loc)


            command = utils.Command(exect_str)
            command.run(timeout=20)
            sys.stdout.write("hhblits is running ......................")
            sys.stdout.flush()
            while not os.path.exists(output_oa3m_homologues_file):
                time.sleep(1)
            logging.info("Output file: %s\n" % output_oa3m_homologues_file)
        else:
            logging.warning("hhblits failed, input file not found: %s\n" % tmp_protein_fasta)

    tmp_file_handle.close()
