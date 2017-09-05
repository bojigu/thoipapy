from Bio.Blast import NCBIWWW
from thoipapy import mtutiles as utils
import os

def download_homologous_from_ncbi(pathdict, set_, logging):
    """From the list of proteins in csv format or single protein in setting file, begins downloading homologues from the NCBI nr database.

    """
    tmp_list_loc=set_["list_of_tmd_start_end"]
    xml_file_loc=set_["xml_file_folder"]
    evalue_score=set_["e_value"]
    evalue_cutoff=set_["e_value_cutoff"]
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


        if(set_["multiple_tmp_simultaneous"]):
            tmp_file_handle=open(tmp_list_loc,'r')
            for row in tmp_file_handle:
                tmp_protein_acc=row.strip().split("\t")[0]
                #differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

                if len(tmp_protein_acc) == 6:    #uniprot acc
                    tmp_protein_name=tmp_protein_acc
                else:                  #PDB proteins
                    tmp_protein_name=tmp_protein_acc[0:6]
                #tmp_protein_name = row.strip().split("\t")[0]
                #if tmp_protein_name == "1A24_HUMAN_1_30":
                tmp_protein_fasta = os.path.join(set_["Protein_folder"], "TmBind1/%s.fasta") % tmp_protein_name

                if os.path.isfile(tmp_protein_fasta):  # whether fasta file in folder
                    tmp_protein_fasta_string = open(tmp_protein_fasta).read()
                    # run online server NCBI blastp with biopython module
                    #logging.info('~~~~~~~~~~starting downloading homologous for protein %s') %tmp_protein_name
                    tmp_protein_xml_file = os.path.join(xml_file_loc, "TmBind1/%s.xml") % tmp_protein_name
                    if not os.path.isfile(tmp_protein_xml_file):
                        try:
                            tmp_protein_homologous_xml_handle = NCBIWWW.qblast("blastp", "nr", tmp_protein_fasta_string,expect=evalue_score, hitlist_size=hit_list_size)
                            save_tmp_xml_file = open(tmp_protein_xml_file, "w")
                            save_tmp_xml_file.write(tmp_protein_homologous_xml_handle.read())
                            save_tmp_xml_file.close()
                            tmp_protein_homologous_xml_handle.close()
                        except:
                            print("Query string not found in the CGI context in qblast")


        ##############################################################################################
        #                                                                                            #
        #      download single homologous of tmp_protein_name in the setting file                    #
        #                                                                                            #
        ##############################################################################################

        else:
            tmp_protein_name = set_["tm_protein_name"]
            tmp_protein_fasta=os.path.join(set_["Protein_folder"],"NoRedundPro/%s.fasta") % tmp_protein_name
            tmp_protein_fasta=set_["input_fasta_file"]
            if os.path.isfile(tmp_protein_fasta):  
                tmp_protein_fasta_file= open(tmp_protein_fasta)         #whether fasta file in folder
                tmp_protein_fasta_string=tmp_protein_fasta_file.read()
                #tmp_protein_fasta_file.close()
                #run online server NCBI blastp with biopython module
                logging.info('~~~~~~~~~~starting downloading homologous for protein {}'.format(tmp_protein_name))
                tmp_protein_homologous_xml_handle=NCBIWWW.qblast("blastp", "nr", tmp_protein_fasta_string, expect=evalue_score,hitlist_size=hit_list_size)
                tmp_protein_xml_file = os.path.join(xml_file_loc,"%s.xml") % tmp_protein_name
                save_tmp_xml_file = open(tmp_protein_xml_file, "w")
                save_tmp_xml_file.write(tmp_protein_homologous_xml_handle.read())
                save_tmp_xml_file.close()
                tmp_protein_homologous_xml_handle.close()
                tmp_protein_fasta_file.close()
                #tmp_protein_fasta_file.close()
                #lose(tmp_protein_fasta)
                #tmp_protein_fasta_string.close()
                logging.info('Homologous download was finished')







