from Bio.Blast import NCBIWWW
from thoipapy import mtutiles as utils
import os

def download_homologous_from_ncbi(set_, logging):
    """From the list of proteins in csv format or single protein in setting file, begins downloading homologues from the NCBI nr database.

    """
    evalue_score=set_["e_value"]
    hit_list_size=int(set_["hit_list_size"])
    enough_hard_drive_space = True

    try:
        byteformat = "GB"
        data_harddrive = set_["data_harddrive"]
        size = utils.get_free_space(data_harddrive, byteformat)
        #logging.info('Hard disk remaining space = {}'.format(size))

        if size[0] < 5:
            raise utils.HardDriveSpaceException("Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))

    except utils.HardDriveSpaceException as e:
        logging.warning(e)


    if enough_hard_drive_space:

        ##############################################################################################
        #                                                                                            #
        #      download multiple homologous from tmp lists                                           #
        #                                                                                            #
        ##############################################################################################

        tmp_list_loc = set_["list_of_tmd_start_end"]
        tmp_file_handle = open(tmp_list_loc, 'r')
        # skip header
        next(tmp_file_handle)
        for row in tmp_file_handle:
            tmp_protein_name = row.strip().split(",")[0]

            tmp_protein_fasta = os.path.join(set_["Protein_folder"], "%s.surr20.fasta") % tmp_protein_name

            if os.path.isfile(tmp_protein_fasta):  # whether fasta file in folder
                with open(tmp_protein_fasta, 'r') as f:
                    tmp_protein_fasta_string=f.read()
                    #tmp_protein_fasta_string = open(tmp_protein_fasta).read()
                    # run online server NCBI blastp with biopython module
                    #logging.info('~~~~~~~~~~starting downloading homologous for protein %s') %tmp_protein_name
                    tmp_protein_xml_file = os.path.join(set_["xml_file_folder"], "%s.xml") % tmp_protein_name
                    if not os.path.isfile(tmp_protein_xml_file):
                        try:
                            tmp_protein_homologous_xml_handle = NCBIWWW.qblast("blastp", "nr", tmp_protein_fasta_string,expect=evalue_score, hitlist_size=hit_list_size)
                            save_tmp_xml_file = open(tmp_protein_xml_file, "w")
                            save_tmp_xml_file.write(tmp_protein_homologous_xml_handle.read())
                            save_tmp_xml_file.close()
                            tmp_protein_homologous_xml_handle.close()
                        except:
                            print("Query string not found in the CGI context in qblast")
                    f.close()
        logging.info("Output file: %s\n" % tmp_protein_xml_file )
        tmp_file_handle.close()
    logging.info('Homologous download was finished')






