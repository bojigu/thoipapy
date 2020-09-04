import csv
import os

import thoipapy


def create_PSSM_from_MSA_mult_prot(s, df_set, logging):
    """ Runs create_PSSM_from_MSA for each protein in a list.

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Returns
    -------

    """
    logging.info('start pssm calculation')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        TMD_seq_pl_surr5 = df_set.loc[i, "TMD_seq_pl_surr5"]
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir,"{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"],s["max_n_gaps_in_TMD_subject_seq"]))
        pssm_csv = os.path.join(s["thoipapy_data_folder"], "features", "pssm", database, "{}.surr{}.gaps{}.pssm.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, TMD_seq, logging)

        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.txt".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))

        pssm_csv_surr5 = os.path.join(s["thoipapy_data_folder"], "features", "pssm", database, "{}.surr5.gaps{}.pssm.csv".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        create_PSSM_from_MSA(path_uniq_TMD_seqs_surr5_for_LIPO, pssm_csv_surr5, acc, TMD_seq_pl_surr5, logging)


def create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, TMD_seq, logging):
    """Creates a PSSM from a multiple sequence alignment in FASTA format.

    Parameters
    ----------
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    pssm_csv : str
        Path to csv file with the PSSM for the TMD region.
    acc : str
        Protein accession (e.g. UniProt, PDB)
    logging : logging.Logger
        Python object with settings for logging to console and file.

    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        thoipapy.utils.make_sure_path_exists(pssm_csv, isfile=True)
        with open(pssm_csv, 'w') as pssm_file_handle:
            mat = []
            writer = csv.writer(pssm_file_handle)
            writer.writerow(['residue_num', 'residue_name', 'A', 'I', 'L', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'M', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P'])
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, "r") as f:
                for line in f.readlines():
                    # strip removes \n. Input currently is not FASTA, no need to check for >
                    mat.append(list(line.strip()))
                    # if not re.search("^>", line):
                    #     mat.append(list(line))

            """ mat is a nested list of all sequences
            [['L', 'G', 'C', 'S', 'A', 'V', 'G', 'G'], ['L', 'A', 'A', 'S', 'A', 'V', 'G', 'P', 'G', 'I', 'G', 'E', 'G'], [.. and so on
            """
            # number of residues in each sequence
            n_residues = len(mat[0])
            if n_residues != len(TMD_seq):
                raise ValueError("Alignment length does not match TMD length. Check that appropriate TMD_seq (or TMD_seq_surr5) has been inserted."
                                 "TMD_seq : {}, first line of alignment = {}".format(TMD_seq, mat[0]))
            # number of sequences in alignment
            n_seqs = len(mat)
            column = []
            # write 20 amino acids as the header of pssm output file
            # pssm_file_handle.write(
            # 'residue'+' '+'A' + ' ' + 'I' + ' ' + 'L' + ' ' + 'V' + ' ' + 'F' + ' ' + 'W' + ' ' + 'Y' + ' ' + 'N' + ' ' + 'C' + ' ' + 'Q' + ' ' + 'M' + ' ' + 'S' + ' ' + 'T' + ' ' + 'D' + ' ' + 'E' + ' ' + 'R' + ' ' + 'H' + ' ' + 'K' + ' ' + 'G' + ' ' + 'P' + '\n')

            #for j in range(0, rowlen - 1):

            # no need to exclude final \n anymore. iterate through number of residues
            for j in range(0, n_residues):
                # iterate through sequences
                for i in range(0, n_seqs):
                    # append the residue in that column
                    # this results in a list of residues at that column position, for all sequences
                    column.append(mat[i][j])
                aa_num = [column.count('A') / n_seqs, column.count('I') / n_seqs, column.count('L') / n_seqs, column.count('V') / n_seqs, column.count('F') / n_seqs,
                          column.count('W') / n_seqs, column.count('Y') / n_seqs, column.count('N') / n_seqs, column.count('C') / n_seqs, column.count('Q') / n_seqs,
                          column.count('M') / n_seqs, column.count('S') / n_seqs, column.count('T') / n_seqs, column.count('D') / n_seqs, column.count('E') / n_seqs,
                          column.count('R') / n_seqs, column.count('H') / n_seqs, column.count('K') / n_seqs, column.count('G') / n_seqs, column.count('P') / n_seqs]
                # add the residue name to the second column
                aa_num.insert(0, TMD_seq[j])
                #aa_num.insert(0, mat[0][j])  #DEPRECATED: Assumes that first sequence is always the original sequence (will break as soon as this is not the case...)
                # add the residue number to the first column
                aa_num.insert(0, j + 1)
                """ the aa_num now looks like this:
                [4, 'S', 0.0, 0.014, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0074, 0.0, 0.0, 0.933, 0.037, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007]
                """
                writer.writerow(aa_num)
                column = []
            logging.info('{} pssm calculation finished ({})'.format(acc, pssm_csv))

    else:
        logging.warning("{} homo_filter_fasta_file does not exist({})".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))