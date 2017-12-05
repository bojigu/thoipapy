# intersect function
def intersect(a, b):
    return list(set(a) & set(b))


import pandas as pd
import numpy as np
import scipy
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import os
import subprocess, threading
import sys
import glob
import numpy as np
import tlabtools as tools
import csv


def Bo_Curve_Create(acc,prob_pos,Lips_score, disrupt_or_closedist,database):
    """
    Create Bo Curve parameter for protein acc and return the output as a dataframe
    Parameters
    ----------
    acc : str, protein name
    prob_pos : the thoipa prediction score for each protein positions
    Lips_score : the LIPS score predicted by LIPS for each protein positions
    disrupt_or_closedist : for experimental ToxR etra data, it is disruption, but for CRYSTAL or NMR data, it is closedistance

    Returns
    -------
    odf: output dataframe

    """
    if database == "crystal" or database == "NMR":
        disrupt_or_closedist = -1 * disrupt_or_closedist  #(it is closest distance and low value means high propencity of interfacial)
    if database == "ETRA":
        disrupt_or_closedist = disrupt_or_closedist      #(it is closest experimental disruption and high value means high propencity of interfacial)
    odf = pd.DataFrame()
    ind = 0
    for i in range(1, 11):
        tm_len = len(disrupt_or_closedist)
        inter_num_total = i
        non_inter_num_total = tm_len - inter_num_total
        thoipa_inter_num_ober = len(
            intersect(np.argsort(np.array((-1) * prob_pos))[0:i], np.argsort(np.array((-1) * disrupt_or_closedist))[0:i]))
        Lips_inter_num_ober = len(
            intersect(np.argsort(np.array(Lips_score))[0:i], np.argsort(np.array((-1) * disrupt_or_closedist))[0:i]))
        pval = scipy.special.comb(inter_num_total, thoipa_inter_num_ober) * scipy.special.comb(non_inter_num_total,
                                                                                               inter_num_total - thoipa_inter_num_ober) / scipy.special.comb(
            tm_len, inter_num_total)
        inter_num_random = 0
        for j in range(1, i + 1):
            random_inter_num_ober = j
            random_non_inter_num = tm_len - inter_num_total
            inter_num_random = inter_num_random + scipy.special.comb(inter_num_total,
                                                                     random_inter_num_ober) * scipy.special.comb(
                random_non_inter_num, inter_num_total - random_inter_num_ober) / scipy.special.comb(tm_len,
                                                                                                    inter_num_total) * j
        odf.set_value(ind, acc, thoipa_inter_num_ober)
        odf.set_value(ind, "sample_size", "Top" + str(i))
        odf.set_value(ind, "parameters", "Ono")
        ind = ind + 1
        odf.set_value(ind, acc, inter_num_random)
        odf.set_value(ind, "sample_size", "Top" + str(i))
        odf.set_value(ind, "parameters", "Rno")
        ind = ind + 1
        odf.set_value(ind, acc, Lips_inter_num_ober)
        odf.set_value(ind, "sample_size", "Top" + str(i))
        odf.set_value(ind, "parameters", "Lpno")
        ind = ind + 1
        odf.set_value(ind, acc, pval)
        odf.set_value(ind, "sample_size", "Top" + str(i))
        odf.set_value(ind, "parameters", "Pono")
        ind = ind + 1
    odf.set_index(["sample_size", "parameters"], inplace=True, drop=True)
    return odf
