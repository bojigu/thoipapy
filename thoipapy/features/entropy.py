import csv
import os

import scipy as sc
from pandas import Series


def entropy_calculation_mult_prot(s, df_set, logging):
    """ Runs entropy_calculation for a set of proteins

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
    """
    logging.info('start entropy calculation')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        # homo_filter_fasta_file = os.path.join(s["data_dir"], "homologues", "a3m",database,"%s.a3m.mem.uniq.2gaps%s") % (acc,s["surres"])
        alignments_dir = os.path.join(s["data_dir"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        entropy_file = os.path.join(s["data_dir"], "features", "entropy", database, "{}.surr{}.gaps{}.uniq.entropy.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, TMD_seq, entropy_file, logging)


def entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, TMD_seq, entropy_file, logging):
    """Calculates conservation of positions using entropy formula.

    S = - sum (Pi + log(Pi))
    See the scipy documentation for more information regarding the entropy formula.
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html

    Note that this code currently considers gaps to be the 21st amino acid.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    entropy_file : str
        Path to csv file with entropy (conservation) data
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):

        if not os.path.isdir(os.path.dirname(entropy_file)):
            os.makedirs(os.path.dirname(entropy_file))

        with open(entropy_file, 'w') as entropy_file_handle:
            mat = []
            with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
                # iterate through each sequence
                for line in f.readlines():
                    # append each seq as a list to nested list, mat
                    mat.append(list(line.strip()))
                    # if not re.search("^>", line):
                    #     mat.append(list(line))
                n_residues = len(mat[0])
                n_seqs = len(mat)
                column = []
            writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            writer.writerow(["residue_num", "residue_name", "conservation"])
            # iterate through each residue position
            for j in range(0, n_residues):
                # iterate through each sequence
                for i in range(0, n_seqs):
                    # add the residue at that column position to list
                    column.append(mat[i][j])
                """
                UP UNTIL NOW, THE CODE MIRRORS EXACTLY WHAT IS DONE FOR THE PSSM CALCULATION
                """

                # convert to series
                column_serie = Series(column)
                # calculates the probabilities
                p_data = column_serie.value_counts() / len(column_serie)
                # input probabilities to get the entropy
                # note that this currently includes gaps "-" as the 21st amino acid
                entropy_orig = sc.stats.entropy(p_data)
                # convert to positive values that are inverted
                # in figures, for less technical readers, "high conservation and polarity" is easier to understand than "low entropy and high polarity"
                conservation = (entropy_orig * -1) + 3
                csv_header_for_ncbi_homologues_file = [j + 1, TMD_seq[j], conservation]
                writer = csv.writer(entropy_file_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(csv_header_for_ncbi_homologues_file)
                # entropy_file_handle.write(mat[0][j]+' '+ str(entropy)+'\n')
                column = []
            # entropy_file_handle.close()
            logging.info('{} entropy_calculation finished ({})'.format(acc, entropy_file))
    else:
        logging.warning("{} entropy_calculation failed. {} input file not found".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))
