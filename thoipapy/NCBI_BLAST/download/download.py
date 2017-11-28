from Bio.Blast import NCBIWWW
from thoipapy import mtutiles as utils
import os
import tarfile
import platform
import time
from time import strftime
from thoipapy.utils import delete_BLAST_xml


def download_homologues_from_ncbi(set_, df_set, logging):
    """Download homologues with NCBI qblast
    :param set_:
    :param logging:
    :return: save homologues into xml files
    """
    evalue_score = set_["expect_value"]
    hit_list_size = int(set_["hit_list_size"])

    OS_description = platform.system()

    # try to detect if there's not enough HD space for the download.
    # currently not working for
    if "Linux" in OS_description or "Windows" in OS_description:
        try:
            byteformat = "GB"
            thoipapy_data_folder = set_["thoipapy_data_folder"]
            size = utils.get_free_space(thoipapy_data_folder, byteformat)
            # logging.info('Hard disk remaining space = {}'.format(size))

            if size[0] < 5:
                raise utils.HardDriveSpaceException(
                    "Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))

        # log the exception, so you actually can see what goes on in the logfile
        except utils.HardDriveSpaceException as e:
            logging.critical(e)
            # now stop all processes and raise an error
            raise utils.HardDriveSpaceException("process stopped")
    else:
        logging.warning("Your system does not seem to be Linux or Windows. Harddrive space check not conducted.")

    ##############################################################################################
    #                                                                                            #
    #      download homologues from protein lists file                                           #
    #                                                                                            #
    ##############################################################################################

    # tmp_list_loc = set_["list_of_tmd_start_end"]
    # tmp_file_handle = open(tmp_list_loc, 'r')
    # # skip header
    # next(tmp_file_handle)
    # for row in tmp_file_handle:
    #     tmp_protein_name = row.strip().split(",")[0]
    #
    #     tmp_protein_fasta = os.path.join(set_["Protein_folder"], database, "%s.surr20.fasta") % tmp_protein_name
    #
    #     if os.path.isfile(tmp_protein_fasta):  # whether fasta file in folder
    #         with open(tmp_protein_fasta, 'r') as f:
    #             tmp_protein_fasta_string = f.read()
    #             print("tmp_protein_fasta_string")
    #             print(tmp_protein_fasta_string)

    for tmp_protein_name in df_set.index:
        TMD_seq_pl_surr = df_set.loc[tmp_protein_name, "TMD_seq_pl_surr"]
        database = df_set.loc[tmp_protein_name, "database"]
        tmp_protein_fasta_string = ">{} TMD add surround 20 residues\n{}".format(tmp_protein_name, TMD_seq_pl_surr)

        # run online server NCBI blastp with biopython module
        blast_xml_file = os.path.join(set_["xml_file_folder"], database, "{}.surr{}.BLAST.xml".format(tmp_protein_name, set_["num_of_sur_residues"]))
        xml_tar_gz = blast_xml_file[:-4] + ".xml.tar.gz"
        xml_txt = blast_xml_file[:-4] + ".BLAST_details.txt"


        if not os.path.isfile(xml_tar_gz):
            run_download = True
            logging.info('{} starting download_homologues_from_ncbi'.format(tmp_protein_name))
        else:
            if set_["rerun_existing_blast_results"]:
                run_download = True
                logging.info('{} starting download_homologues_from_ncbi (EXISTING xml.tar.gz FILE WILL BE OVERWRITTEN)'.format(tmp_protein_name))
            elif set_["rerun_existing_blast_results"] in [False, 0]:
                run_download = False
                logging.info('{} download_homologues_from_ncbi skipped (EXISTING xml.tar.gz FILE)'.format(tmp_protein_name))
                # skip protein
                continue
            else:
                raise ValueError('set_["rerun_existing_blast_results"] does not seem to be True or False')

        if run_download:
            try:
                start = time.clock()
                tmp_protein_homologues_xml_handle = NCBIWWW.qblast("blastp", "nr", tmp_protein_fasta_string,
                                                                   expect=evalue_score,
                                                                   hitlist_size=hit_list_size)
                with open(blast_xml_file, "w") as save_tmp_xml_file:
                    save_tmp_xml_file.write(tmp_protein_homologues_xml_handle.read())

                #save_tmp_xml_file = open(blast_xml_file, "w")
                # save_tmp_xml_file.write(tmp_protein_homologues_xml_handle.read())
                # save_tmp_xml_file.close()
                tmp_protein_homologues_xml_handle.close()
                duration = time.clock() - start

            except:
                logging.warning("{} Query string not found in the CGI context in qblast".format(tmp_protein_name))

            # create an empty text file with the download date
            date = strftime("%Y%m%d")
            with open(xml_txt, "w") as f:
                f.write("acc\t{}\ndownload_date\t{}\ndatabase\tncbi_nr\ne_value\t1\n".format(tmp_protein_name, date))

            with tarfile.open(xml_tar_gz, mode='w:gz') as tar:
                # add the files to the compressed tarfile
                tar.add(blast_xml_file, arcname=os.path.basename(blast_xml_file))
                tar.add(xml_txt, arcname=os.path.basename(xml_txt))

            delete_BLAST_xml(blast_xml_file)

            logging.info("Output file: {}. (time taken = {:0.3f} min)".format(xml_tar_gz, duration / 60))


    logging.info('homologues download was finished')
