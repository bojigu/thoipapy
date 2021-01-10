import csv
import math
import os
import re

import thoipapy


def LIPS_score_calculation_mult_prot(s, df_set, logging):
    """Run LIPS_score_calculation for a list of proteins.

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
    logging.info('start lips score calculation')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        # LIPS_input_file = os.path.join(s["data_dir"], "homologues", "a3m",database, "%s.mem.lips.input%s") % (acc,s["surres"])
        alignments_dir = os.path.join(s["data_dir"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, s["num_of_sur_residues"]))

        if os.path.isfile(path_uniq_TMD_seqs_no_gaps_for_LIPS):
            # LIPS_output_file = os.path.join(s["data_dir"], "features", "lips_score", "zpro/NoRedundPro/%s.mem.lips.output") % acc
            # path_uniq_TMD_seqs_no_gaps_for_LIPS = os.path.join(alignments_dir, "{}.surr{}.gaps0.uniq.for_LIPS.txt".format(acc, s["num_of_sur_residues"]))

            # LIPS_output_file = os.path.join(s["data_dir"], "features", "lips_score", database, "%s.mem.lips.output%s") % (acc, s["surres"])

            LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, s["num_of_sur_residues"]))

            LIPS_score_calculation(path_uniq_TMD_seqs_no_gaps_for_LIPS, LIPS_output_file)
        else:
            logging.warning("{} path_uniq_TMD_seqs_no_gaps_for_LIPS not found".format(acc))


def LIPS_score_calculation(input_seq_file, LIPS_output_file):
    """Python version of the LIPS algorithm by Adamian and Liang (2016) Prediction of transmembrane helix orientation in polytopic membrane protenis.

    This script should give exactly the same output as the original perl algorithm.

    Parameters
    ----------
    input_seq_file : str
        Path to text file with a list of sequences
    LIPS_output_file : str
        Path to file with LIPS output result.
    """
    # use different reference to the print function, to aid in finding areas that are being debugged
    p_r_i_n_t = print

    thoipapy.utils.make_sure_path_exists(LIPS_output_file, isfile=True)

    with open(input_seq_file, "r") as file:
        sequence = ' '.join(line.strip() for line in file)

    with open(LIPS_output_file, "w") as LIPS_output_file_handle:
        n = 0
        sump = 0
        sumlip = 0
        sume = {}  # the sum of entropy for each one of the seven surfaces
        sumf = 0
        sumim = {}  # the sum of lipophilicity for each surfaces
        aanum = {}  # the number of residues for each one of seven surfaces
        resnum = 1
        amino = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        ###propi and propm are both from TMLIP scale, propi means the residue lipophilicity propensity in membrane headgroup region
        ##and propm is the residue lipophilicity propensity in hydrophobic core region
        ##in TMLIP scale paper, the membrane headgroup regions is defined as the first (tmlen/5) and (tmlen-tmlen/5) residues which
        ##more likely the membrane bilayer
        ##while the other residues are defined as hydrophobic core region

        propi = {
            'A': 0.71,
            'R': 1.47,
            'N': 0.96,
            'D': 1.20,
            'C': 1.16,
            'Q': 0.61,
            'E': 0.90,
            'G': 0.48,
            'H': 0.82,
            'I': 1.11,
            'L': 1.18,
            'K': 2.38,
            'M': 1.38,
            'F': 1.57,
            'P': 0.99,
            'S': 0.69,
            'T': 0.72,
            'W': 2.45,
            'Y': 1.23,
            'V': 0.98
        }

        propm = dict(
            A=0.82,
            R=0.18,
            N=0.19,
            D=0.29,
            C=1.01,
            Q=0.26,
            E=0.19,
            G=0.35,
            H=0.12,
            I=1.88,
            L=1.71,
            K=0.42,
            M=1.02,
            F=1.97,
            P=0.65,
            S=0.55,
            T=0.66,
            W=1.65,
            Y=0.94,
            V=1.77
        )

        tmp = sequence.split(' ')
        nrow = len(tmp)
        ncol = len(tmp[0])
        bnum = ncol / 5
        oc = {}
        prob = {}
        entropy = {}  ##residue entropy
        exp_entropy = {}  # the final exponential entropy
        lips = {}  ##the lipophilicity score

        for i in range(nrow):
            for j in range(ncol):
                residue = tmp[i][j]
                res_j = ' '.join((residue, str(j)))
                if (res_j in oc.keys()):
                    oc[res_j] = oc[res_j] + 1
                else:
                    oc[res_j] = 1

        for j in range(ncol):
            for res in amino:
                if (' '.join((res, str(j))) in oc):
                    prob[res] = oc[' '.join((res, str(j)))] / nrow
                    if (j in entropy.keys()):
                        entropy[j] = entropy[j] + prob[res] * math.log(prob[res])  # the entropy calculation
                    else:
                        entropy[j] = prob[res] * math.log(prob[res])
                    if ((j <= bnum) or (j > ncol - bnum)):  ###here is the membrane headgroup residues
                        if (j in lips.keys()):
                            lips[j] = lips[j] + prob[res] * propi[res]
                        else:
                            lips[j] = prob[res] * propi[res]
                    else:  ###here is the hydrophobic region residues
                        if (j in lips.keys()):
                            lips[j] = lips[j] + prob[res] * propm[res]
                        else:
                            lips[j] = prob[res] * propm[res]
            exp_entropy[j] = 2.718 ** ((-1) * entropy[j])  # expontional entropy

        for j in sorted(exp_entropy):
            res = tmp[0][j]
            m = resnum + j
            sump = sump + exp_entropy[j]
            sumlip = sumlip + lips[j]

        for i in range(4):  # for the first 4 surfaces
            p_r_i_n_t("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
            # LIPS_output_file_handle.write("SURFACE", "%s" % i, ":")
            j = i
            while j < ncol:
                res = tmp[0][j]
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[j]  # the sum of lipophilicity for surface i
                else:
                    sumim[i] = lips[j]
                prop = lips[j]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[j]  # the sum of entropy for surface i
                else:
                    sume[i] = exp_entropy[j]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1  # the sum of the residue numbers for surface i
                else:
                    aanum[i] = 1
                rn = j + resnum
                # r3=residuename123(res)
                p_r_i_n_t("%3s" % rn, res, "%6.3f" % prop,
                          "%6.3f" % exp_entropy[j],
                          file=LIPS_output_file_handle)  # print residue information which is in surface i
                # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prop,"%6.3f" % exp_entropy[j])
                k = j + 3
                while (k <= j + 4):  # here add the the residues of i+3 and i+4 into surface i to form heptad repeat
                    if (k < ncol):
                        res = tmp[0][k]
                        # r3=residuename123(res)
                        if (i in sumim.keys()):
                            sumim[i] = sumim[i] + lips[k]
                        else:
                            sumim[i] = lips[k]
                        prob = lips[k]
                        if (i in sume.keys()):
                            sume[i] = sume[i] + exp_entropy[k]
                        else:
                            sume[i] = exp_entropy[k]
                        if (i in aanum.keys()):
                            aanum[i] = aanum[i] + 1
                        else:
                            aanum[i] = 1
                        rn = k + resnum
                        p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
                                  file=LIPS_output_file_handle)
                        # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
                j = j + 7
        for i in range(4, 7):  # for surfaces from 4 to 6
            p_r_i_n_t("SURFACE", "%s" % i, ":", file=LIPS_output_file_handle)
            # LIPS_output_file_handle.write("SURFACE", "%s" % i, ":")
            j = i
            while j < ncol:
                res = tmp[0][j]
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[j]
                else:
                    sumim[i] = lips[j]
                prob = lips[j]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[j]
                else:
                    sume[i] = exp_entropy[j]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1
                else:
                    aanum[i] = 1
                rn = j + resnum
                # r3=residuename123(res)
                p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j], file=LIPS_output_file_handle)
                # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j])
                k = j + 3
                while (k <= j + 4):
                    if (k < ncol):
                        res = tmp[0][k]
                        # r3=residuename123(res)
                        if (i in sumim.keys()):
                            sumim[i] = sumim[i] + lips[k]
                        else:
                            sumim[i] = lips[k]
                        prob = lips[k]
                        if (i in sume.keys()):
                            sume[i] = sume[i] + exp_entropy[k]
                        else:
                            sume[i] = exp_entropy[k]
                        if (i in aanum.keys()):
                            aanum[i] = aanum[i] + 1
                        else:
                            aanum[i] = 1
                        rn = k + resnum
                        p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k],
                                  file=LIPS_output_file_handle)
                        # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                    k = k + 1
                j = j + 7
            k = i - 4
            while (k <= i - 3):  # here adding residues at the first 7 positions
                if (k < ncol):
                    res = tmp[0][k]
                    # r3=residuename123(res)
                    if (i in sumim.keys()):
                        sumim[i] = sumim[i] + lips[k]
                    else:
                        sumim[i] = lips[k]
                    prob = lips[k]
                    if (i in sume.keys()):
                        sume[i] = sume[i] + exp_entropy[k]
                    else:
                        sume[i] = exp_entropy[k]
                    if (i in aanum.keys()):
                        aanum[i] = aanum[i] + 1
                    else:
                        aanum[i] = 1
                    rn = k + resnum
                    p_r_i_n_t("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k], file=LIPS_output_file_handle)
                    # LIPS_output_file_handle.write("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
                k = k + 1
        p_r_i_n_t("SURFACE LIPOPHILICITY ENTROPY   LIPS", file=LIPS_output_file_handle)
        # LIPS_output_file_handle.write("SURFACE LIPOPHILICITY ENTROPY   LIPS")
        for i in sumim.keys():
            avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
            avpim = avpim * 2
            ave = sume[i] / aanum[i]  # average entropy for surface i
            peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
            p_r_i_n_t("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
                      "%8.3f" % peim,
                      file=LIPS_output_file_handle)  # print seven surfaces and see which surface with lowewst LIPS score
            # LIPS_output_file_handle.write("%s" % i, "%10.3f" % avpim, "%8.3f" % ave, "%8.3f" % peim)


def parse_LIPS_score_mult_prot(s, df_set, logging):
    """Runs parse_LIPS_score for a list of sequences.

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
    logging.info('start parsing lips output to cons and lips scores')

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = os.path.join(s["data_dir"], "homologues", "alignments", database)
        LIPS_output_file = os.path.join(alignments_dir, "{}.surr{}.LIPS_output.csv".format(acc, s["num_of_sur_residues"]))
        LIPS_parsed_csv = os.path.join(s["data_dir"], "features", "lips_score", database, "{}.surr{}.LIPS_score_parsed.csv".format(acc, s["num_of_sur_residues"]))
        parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging)


def parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging):
    """Parse the LIPS output to create a CSV with features for input in machine learning algorithm.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt, PDB)
    LIPS_output_file : str
        Path to file with LIPS output result.
    LIPS_parsed_csv : str
        Path to csv with LIPS output organised into features for machine learning.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """

    thoipapy.utils.make_sure_path_exists(LIPS_parsed_csv, isfile=True)

    if os.path.isfile(LIPS_output_file):
        # try:
        surface_num = 0
        surface_lips = 100  ##100 is an initialized big number assuming lips score will not bigger than this number
        with open(LIPS_output_file, "r") as LIPS_output_handle:
            with open(LIPS_parsed_csv, "w") as LIPS_parsed_csv_handle:
                i = 0
                array = []
                dict = {}
                for row in LIPS_output_handle:
                    if re.search(r"^\s+\d+\s+[A-Z]", row):
                        array = row.split()
                        if not int(array[0]) in dict:
                            dict[int(array[0])] = " ".join([array[1], array[2], array[3]])

                    if re.search("^\d{1}\s+", row):
                        surface_num1 = row.split()[0]
                        surface_lips1 = row.split()[3]
                        if (float(surface_lips1) < float(surface_lips)):
                            surface_lips = surface_lips1
                            surface_num = surface_num1
                LIPS_output_handle.close()

                surface_find = 0
                dict1 = {}
                LIPS_output_handle = open(LIPS_output_file, "r")
                for row in LIPS_output_handle:
                    if re.search(r"^SURFACE\s" + surface_num, row):
                        surface_find = 1
                        continue
                    if surface_find == 1 and re.search(r"^\s+\d+\s+[A-Z]", row):
                        array = row.split()
                        if not int(array[0]) in dict1:
                            dict1[int(array[0])] = " ".join([array[1], array[2], array[3]])
                    else:
                        surface_find = 0
                LIPS_output_handle.close()

                writer = csv.writer(LIPS_parsed_csv_handle, delimiter=',', quotechar='"', lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                writer.writerow(["residue_num", "residue_name", "LIPS_polarity", "LIPS_entropy", "LIPS_surface"])
                for k, v in sorted(dict.items()):
                    v1 = v.split()
                    v1.insert(0, k)
                    if k in dict1:
                        v1.insert(4, 1)
                    else:
                        v1.insert(4, 0)
                    csv_header_for_cons_lips_score_file = v1
                    writer.writerow(csv_header_for_cons_lips_score_file)
                LIPS_parsed_csv_handle.close()
                logging.info('{} lips score parse finished ({})'.format(acc, LIPS_parsed_csv))
    else:
        logging.warning("{} LIPS_output_file not found.")
