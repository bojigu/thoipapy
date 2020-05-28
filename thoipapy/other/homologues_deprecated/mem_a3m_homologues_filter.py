import os
import re
from difflib import SequenceMatcher

from thoipapy.features.lipophilicity import calc_lipophilicity


def mem_a3m_homologues_filter(s,logging):
    p_r_i_n_t = print
    tmp_list_loc = s["list_of_tmd_start_end"]
    tmp_file_handle = open(tmp_list_loc, 'r')
    # skip header
    next(tmp_file_handle)
    for row in tmp_file_handle:
        acc = row.strip().split(",")[0][0:5]
        TMD_start = int(row.strip().split(",")[4]) + 1
        TMD_end = int(row.strip().split(",")[4]) + (int(row.strip().split(",")[3]) - int(row.strip().split(",")[2]) + 1)
        #return TMD_end
        homo_a3m_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.parse.a3m") % acc
        if os.path.isfile(homo_a3m_file):

            homo_a3m_file_handle=open(homo_a3m_file,"r")
            homo_filter_file=os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.a3m.mem.uniq.2gaps") %acc
            homo_filter_file_handle = open(homo_filter_file,"w")
            homo_mem_lips_input_file = os.path.join(s["thoipapy_data_folder"], "homologues", "a3m", "SinglePassTmd/%s.surr20.mem.lips.input") %acc
            homo_mem_lips_input_file_handle = open(homo_mem_lips_input_file, "w")
            logging.info("starting parsing a3m file: %s\n" %homo_filter_file)
            i = 0
            tm_query=""
            for line in homo_a3m_file_handle:
                tm_str = line[(TMD_start - 1):TMD_end]
                if i == 0:
                    tm_query = tm_str
                    i = i + 1
                    p_r_i_n_t("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                    p_r_i_n_t("{}".format(tm_str), file=homo_filter_file_handle)
                    continue
                mean_hydrophobicity = calc_lipophilicity(tm_str)
                ratio = SequenceMatcher(None, tm_query, tm_str).ratio()
                if not re.search("-", tm_str) and not re.search("X",
                                                                tm_str) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"]and mean_hydrophobicity < s["max_hydrophilicity_Hessa"]:  ##No X and gap in each alignment
                    p_r_i_n_t("{}".format(tm_str), file=homo_mem_lips_input_file_handle)
                gap_num = tm_str.count("-")
                if (gap_num <= 3 and not re.search("X",
                                                   tm_str) and ratio >= s["min_identity_of_TMD_seq"] and ratio < s["max_identity_of_TMD_seq"] and mean_hydrophobicity < s["max_hydrophilicity_Hessa"]):  # gap number le 3 and no X in each alignment
                    p_r_i_n_t("{}".format(tm_str), file=homo_filter_file_handle)
                    # homo_filter_file_handle.write(line)
                    # homo_mem_lips_input_file_handle.write(tm_str
            homo_filter_file_handle.close()
            homo_mem_lips_input_file_handle.close()
        homo_a3m_file_handle.close()
    tmp_file_handle.close()