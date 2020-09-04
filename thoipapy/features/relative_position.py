import csv
import os
import re

import pandas as pd

import thoipapy
from thoipapy.utils import normalise_0_1


def calc_relative_position_mult_prot(s, df_set, logging):
    """calculate the residue relative position on the TMD

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
    logging.info('start to calculate the relative positions')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        #TMD_start = int(row.strip().split(",")[2])
        #seqlen = int(row.strip().split(",")[1])
        TMD_start = df_set.loc[i, "TMD_start"]
        tm_seq = df_set.loc[i, "full_seq"]
        seqlen = df_set.loc[i, "seqlen"]

        relative_position_file = os.path.join(s["thoipapy_data_folder"], "features", "relative_position", database, "%s.relative_position%s.csv") % (acc, s["surres"])
        thoipapy.utils.make_sure_path_exists(relative_position_file, isfile=True)
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_for_PSSM_FREECONTACT = os.path.join(alignments_dir, "{}.surr{}.gaps{}.uniq.for_PSSM_FREECONTACT.txt".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

        calc_relative_position(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, relative_position_file, TMD_start, seqlen, logging)


def calc_relative_position(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, relative_position_file, TMD_start, seqlen, logging):
    """Calculate the residue relative position on the TMD

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT : str
        Path to text file with list of TMD sequences
    relative_position_file : str
        Path to csv file with the relative position of each residue in the TMD and full protein sequence
    TMD_start : int
        Start of TMD in full sequence
    seqlen : int
        Length of full sequence
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
        relative_position_file_handle = open(relative_position_file, 'w')
        mat = []
        writer = csv.writer(relative_position_file_handle, delimiter=',', quotechar='"',
                            lineterminator='\n',
                            quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
        writer.writerow(['residue_num', 'residue_name', 'RelPos_TMD', 'RelPos_fullseq'])

        with open(path_uniq_TMD_seqs_for_PSSM_FREECONTACT) as f:
            mat = []
            for line in f.readlines():
                if not re.search("^>", line):
                    mat.append(list(line))
            tm_seq = mat[0]
            tm_len = len(tm_seq)
            for i in range(1, tm_len):
                RelPos_TMD = i / tm_len
                RelPos_fullseq = (i + TMD_start - 1) / seqlen
                residue_num = i
                residue_name = tm_seq[i - 1]
                writer.writerow([residue_num, residue_name, RelPos_TMD, RelPos_fullseq])
        relative_position_file_handle.close()
        logging.info('{} relative position calculation finished ({})'.format(acc, relative_position_file))
        dfrp = pd.read_csv(relative_position_file, index_col=0)
        dfrp["residue_depth"] = (normalise_0_1(1 - abs(dfrp.RelPos_TMD - 0.5))[0]).round(1)
        dfrp.to_csv(relative_position_file)

    else:
        logging.warning("{} calc_relative_position failed, file not found ({})".format(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT))