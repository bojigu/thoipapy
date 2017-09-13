import os
import subprocess,threading
import time
import thoipapy.utils as utils


def parse_a3m_alignment(set_, logging):
    logging.info('start parsing a3m file')

    tmp_list_loc = set_["list_of_tmd_start_end"]
    #if (set_["multiple_tmp_homo_download"]):
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        tmp_protein_name = row.strip().split(",")[0]

        output_oa3m_homologous_file = os.path.join(set_["Protein_folder"], "%s.fasta.surr20.a3m") % tmp_protein_name
        output_oa3m_homologous_processed_file = os.path.join(set_["output_oa3m_homologous"],"%s.surr20.parse.a3m") % tmp_protein_name
        exect_str1 = "grep -v '^>' {a3m_file} |sed 's/[a-z]//g' >{a3m_processed_file}".format(a3m_file=output_oa3m_homologous_file,
                                                                                              a3m_processed_file=output_oa3m_homologous_processed_file)

        if os.path.isfile(output_oa3m_homologous_file):
            command1 = utils.Command(exect_str1)
            command1.run(timeout=20)
            logging.info("a3m processed Output file: %s\n" % output_oa3m_homologous_processed_file)
    tmp_file_handle.close()





def download_homologous_with_hhblits(set_,logging):
    logging.info('start downloading protein homologous with hhblits')
    hhblits_loc = set_["hhblits_dir"]
    uniprot_database_loc=set_["uniprot_database_dir"]

    tmp_list_loc = set_["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    #skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        tmp_protein_name = row.strip().split(",")[0]
        tmp_protein_fasta = os.path.join(set_["Protein_folder"], "%s.surr20.fasta") % tmp_protein_name
        if os.path.isfile(tmp_protein_fasta):
            output_oa3m_homologous_file = os.path.join(set_["Protein_folder"], "%s.fasta.surr20.a3m") % tmp_protein_name
            if tmp_protein_name=="QuePro" and os.path.isfile(output_oa3m_homologous_file):
                os.system("rm %s" % output_oa3m_homologous_file)          ###remove QuePro.fasta.a3m that is the previous queried protein

            exect_str = "qsub -b y -q all.q -N 'hhblitsjob' {hhblits} -i {fasta_input} -oa3m {a3m_output} -d {uniprot_database} -Z 999999999 " \
                        "-B 999999999 -maxfilt 999999999 -id 99 -diff inf".format(hhblits=hhblits_loc, fasta_input=tmp_protein_fasta,
                                                                                  a3m_output=output_oa3m_homologous_file,uniprot_database=uniprot_database_loc)


            command = utils.Command(exect_str)
            command.run(timeout=20)
            print("hhblits is running ......................")
            while not os.path.exists(output_oa3m_homologous_file):
                time.sleep(1)
            logging.info("Output file: %s\n" % output_oa3m_homologous_file)
    tmp_file_handle .close()
