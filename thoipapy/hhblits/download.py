import os
import subprocess,threading
import time


def parse_a3m_alignment(pathdict, set_, logging):
    logging.info('start parsing a3m file')
    if not set_["multiple_tmp_simultaneous"]:
        tmp_protein_name = set_["tm_protein_name"]
        output_oa3m_homologous_file = os.path.join(set_["output_oa3m_homologous"], "%s.fasta.a3m") % tmp_protein_name
        output_oa3m_homologous_processed_file = os.path.join(set_["output_oa3m_homologous"], "%s.a3m") % tmp_protein_name
        exect_str1 = "grep -v '^>' {a3m_file} |sed 's/[a-z]//g' >{a3m_processed_file}".format(a3m_file=output_oa3m_homologous_file, a3m_processed_file=output_oa3m_homologous_processed_file)
        print(exect_str1)

        class Command(object):
            def __init__(self, cmd):
                self.cmd = cmd
                self.process = None

            def run(self, timeout):
                def target():
                    print('Thread started')
                    #self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
                    subprocess.call(self.cmd, shell=True)
                    # self.process.communicate()
                    # print(self.process.communicate())
                    print('Thread finished')

                thread = threading.Thread(target=target)
                thread.start()

                thread.join(timeout)
                if thread.is_alive():
                    print('Terminating process')
                    self.process.terminate()
                    thread.join()
                #print(self.process.returncode)

        if os.path.isfile(output_oa3m_homologous_file):
            command1 = Command(exect_str1)
            command1.run(timeout=40)
            logging.info("a3m processed Output file: %s\n" % output_oa3m_homologous_processed_file)
            #os.system("rm %s" % output_oa3m_homologous_file)


    # # ######for multiple download
    else:
        tmp_list_loc = set_["list_of_tmd_start_end"]
        #if (set_["multiple_tmp_homo_download"]):
        tmp_file_handle = open(tmp_list_loc, 'r')
        for row in tmp_file_handle:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tmp_protein_name = tmp_protein_acc
            else:  # PDB proteins
                tmp_protein_name = tmp_protein_acc[0:6]
            output_oa3m_homologous_file = os.path.join(set_["output_oa3m_homologous"], "TMPAD/%s.fasta.a3m") % tmp_protein_name
            output_oa3m_homologous_processed_file = os.path.join(set_["output_oa3m_homologous"],"TMPAD/%s.a3m") % tmp_protein_name
            exect_str1 = "grep -v '^>' {a3m_file} |sed 's/[a-z]//g' >{a3m_processed_file}".format(a3m_file=output_oa3m_homologous_file, a3m_processed_file=output_oa3m_homologous_processed_file)

            class Command(object):
                def __init__(self, cmd):
                    self.cmd = cmd
                    self.process = None

                def run(self, timeout):
                    def target():
                        print('Thread started')
                        #self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
                        subprocess.call(self.cmd, shell=True)
                        # self.process.communicate()
                        # print(self.process.communicate())
                        print('Thread finished')

                    thread = threading.Thread(target=target)
                    thread.start()

                    thread.join(timeout)
                    if thread.is_alive():
                        print('Terminating process')
                        self.process.terminate()
                        thread.join()
                    #print(self.process.returncode)

            if os.path.isfile(output_oa3m_homologous_file):
                command1 = Command(exect_str1)
                command1.run(timeout=20)
                logging.info("a3m processed Output file: %s\n" % output_oa3m_homologous_processed_file)






def download_homologous_with_hhblits(pathdict, set_, logging):
    logging.info('start downloading protein homologous with hhblits')
    hhblits_loc = set_["hhblits_dir"]
    uniprot_database_loc=set_["uniprot_database_dir"]


    ##############################################################################################
    #                                                                                            #
    #      download single homologous of tmp_protein_name in the setting file                    #
    #                                                                                            #
    ##############################################################################################

    if not set_["multiple_tmp_simultaneous"]:
        tmp_protein_name = set_["tm_protein_name"]
        #tmp_protein_fasta = os.path.join(set_["Protein_folder"], "experimental_protein/%s.fasta") % tmp_protein_name
        tmp_protein_fasta=set_["input_fasta_file"]
        if os.path.isfile(tmp_protein_fasta):
            output_oa3m_homologous_file = os.path.join(set_["output_oa3m_homologous"],"%s.fasta.a3m") % tmp_protein_name
            if os.path.isfile(output_oa3m_homologous_file):
                os.system("rm %s" % output_oa3m_homologous_file)          ###remove QuePro.fasta.a3m that is the previous queried protein
            exect_str = "qsub -b y -q all.q {hhblits} -i {fasta_input} -oa3m {a3m_output} -d {uniprot_database} -Z 999999999 -B 999999999 -maxfilt 999999999 -id 99 -diff inf".format(hhblits=hhblits_loc, fasta_input=tmp_protein_fasta, a3m_output=output_oa3m_homologous_file,uniprot_database=uniprot_database_loc)
            print(exect_str)
            class Command(object):
                def __init__(self, cmd):
                    self.cmd = cmd
                    self.process = None

                def run(self, timeout):
                    def target():
                        print('Thread started')
                        #self.process = subprocess.Popen(self.cmd,shell=True,stdout=subprocess.PIPE)
                        subprocess.call(self.cmd,shell=True)
                        #self.process.communicate()
                        # print(self.process.communicate())
                        print('Thread finished')

                    thread = threading.Thread(target=target)
                    thread.start()

                    thread.join(timeout)
                    if thread.is_alive():
                        print('Terminating process')
                        self.process.terminate()
                        thread.join()
                    #print(self.process.returncode)

            command = Command(exect_str)
            # command=Command(exect_str)
            command.run(timeout=20)                       ###since hhblits requres more than 10 minutes to finish, maybe later we should consider using qsub to the server
            # command.run(timeout=1)
            # command=mtutils.Command(exect_str)
            # command.run(timeout=120)
            print("hhblits is running ......................")
            while not os.path.exists(output_oa3m_homologous_file):
                time.sleep(1)
            logging.info("Output file: %s\n" % output_oa3m_homologous_file)



    ##############################################################################################
    #                                                                                            #
    #      download multiple homologous from tmp lists with hhblits                              #
    #                                                                                            #
    ##############################################################################################


    else:
        tmp_list_loc = set_["list_of_tmd_start_end"]
        #if (set_["multiple_tmp_homo_download"]):
        tmp_file_handle = open(tmp_list_loc, 'r')
        for row in tmp_file_handle:
            tmp_protein_acc = row.strip().split("\t")[0]
            # differ uniprot acc from PDB proteins (contains PDB name,chain name and TMD order, like 1ft8_A1)

            if len(tmp_protein_acc) == 6:  # uniprot acc
                tmp_protein_name = tmp_protein_acc
            else:  # PDB proteins
                tmp_protein_name = tmp_protein_acc[0:6]
            tmp_protein_fasta = os.path.join(set_["Protein_folder"], "TMPAD/%s.fasta") % tmp_protein_name
            if os.path.isfile(tmp_protein_fasta):
                output_oa3m_homologous_file = os.path.join(set_["output_oa3m_homologous"], "TMPAD/%s.fasta.a3m") % tmp_protein_name
                #exect_str = "{hhblits} -i {fasta_input} -oa3m {a3m_output} -d {uniprot_database} -Z 999999999 -B 999999999 -maxfilt 999999999 -id 99 -diff inf".format(hhblits=hhblits_loc, fasta_input=tmp_protein_fasta, a3m_output=output_oa3m_homologous_file_handle,uniprot_database=uniprot_database_loc)
                exect_str = "qsub -b y -q all.q -N 'hhblitsjob' {hhblits} -i {fasta_input} -oa3m {a3m_output} -d {uniprot_database} -Z 999999999 -B 999999999 -maxfilt 999999999 -id 99 -diff inf".format(hhblits=hhblits_loc, fasta_input=tmp_protein_fasta, a3m_output=output_oa3m_homologous_file,uniprot_database=uniprot_database_loc)
                print(exect_str)
                class Command(object):
                    def __init__(self, cmd):
                        self.cmd = cmd
                        self.process = None

                    def run(self, timeout):
                        def target():
                            print('Thread started')
                            #self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
                            subprocess.call(self.cmd, shell=True)
                            # print(self.process.communicate())
                            print('Thread finished')

                        thread = threading.Thread(target=target)
                        thread.start()

                        thread.join(timeout)
                        if thread.is_alive():
                            print('Terminating process')
                            self.process.terminate()
                            thread.join()
                        #print(self.process.returncode)

                command = Command(exect_str)
                # command=Command(exect_str)
                command.run(timeout=2000)
                # command.run(timeout=1)
                # command=mtutils.Command(exect_str)
                # command.run(timeout=120)
                logging.info("Output file: %s\n" % output_oa3m_homologous_file)

