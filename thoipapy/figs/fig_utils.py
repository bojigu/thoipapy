# intersect function
def intersect(a, b):
    return list(set(a) & set(b))


import pandas as pd
import scipy
import os
import thoipapy

#
# def Bo_Curve_Create(acc, prob_pos, Lips_score, disrupt_or_closedist, database):
#     """
#     Create Bo Curve parameter for protein acc and return the output as a dataframe
#     Parameters
#     ----------
#     acc : str, protein name
#     prob_pos : the thoipa prediction score for each protein positions
#     Lips_score : the LIPS score predicted by LIPS for each protein positions
#     disrupt_or_closedist : for experimental ToxR etra data, it is disruption, but for CRYSTAL or NMR data, it is closedistance
#
#     Returns
#     -------
#     odf: output dataframe
#
#     """
#     if database == "crystal" or database == "NMR":
#         disrupt_or_closedist = -1 * disrupt_or_closedist  #(it is closest distance and low value means high propencity of interfacial)
#     if database == "ETRA":
#         disrupt_or_closedist = disrupt_or_closedist      #(it is closest experimental disruption and high value means high propencity of interfacial)
#     odf = pd.DataFrame()
#     ind = 0
#     for sample_size in range(1, 11):
#         # length of residues tested (usually full TMD length, unless some are excluded)
#         tm_len = len(disrupt_or_closedist)
#         non_sample_size = tm_len - sample_size
#
#         observed_overlap = len(intersect(np.argsort(np.array((-1) * prob_pos))[0:sample_size], np.argsort(np.array((-1) * disrupt_or_closedist))[0:sample_size]))
#         LIPS_observed_overlap = len(intersect(np.argsort(np.array(Lips_score))[0:sample_size], np.argsort(np.array((-1) * disrupt_or_closedist))[0:sample_size]))
#         pval = scipy.special.comb(sample_size, observed_overlap) * scipy.special.comb(non_sample_size, sample_size - observed_overlap) / scipy.special.comb(tm_len, sample_size)
#         inter_num_random = 0
#         for j in range(1, sample_size + 1):
#             random_inter_num_ober = j
#             random_non_inter_num = tm_len - sample_size
#             inter_num_random = inter_num_random + scipy.special.comb(sample_size,random_inter_num_ober) * scipy.special.comb(random_non_inter_num, sample_size - random_inter_num_ober) / scipy.special.comb(tm_len,
#                                                                                                     sample_size) * j
#         odf.set_value(ind, acc, observed_overlap)
#         odf.set_value(ind, "sample_size", "Top" + str(sample_size))
#         odf.set_value(ind, "parameters", "observed_overlap")
#         ind = ind + 1
#         odf.set_value(ind, acc, inter_num_random)
#         odf.set_value(ind, "sample_size", "Top" + str(sample_size))
#         odf.set_value(ind, "parameters", "random_overlap")
#         ind = ind + 1
#         odf.set_value(ind, acc, LIPS_observed_overlap)
#         odf.set_value(ind, "sample_size", "Top" + str(sample_size))
#         odf.set_value(ind, "parameters", "LIPS_observed_overlap")
#         ind = ind + 1
#         odf.set_value(ind, acc, pval)
#         odf.set_value(ind, "sample_size", "Top" + str(sample_size))
#         odf.set_value(ind, "parameters", "p_value_from_obs_overlap")
#         ind = ind + 1
#     odf.set_index(["sample_size", "parameters"], inplace=True, drop=True)
#     return odf



def calc_best_overlap(acc_db, df, experiment_col="interface_score", pred_col="THOIPA"):
    """
    Create Bo Curve parameter for protein acc_db and return the output as a dataframe
    Parameters
    ----------
    acc_db : str, protein name
    prob_pos : the thoipa prediction score for each protein positions
    Lips_score : the LIPS score predicted by LIPS for each protein positions
    interface_score_arr : for experimental ToxR etra data, it is disruption, but for CRYSTAL or NMR data, it is closedistance

    Returns
    -------
    odf: output dataframe

    """
    tm_len = df.shape[0]

    if experiment_col not in df.columns:
        raise IndexError("{} {} is not in the columns.\n column list = {}".format(acc_db, experiment_col, df.columns))
    if pred_col not in df.columns:
        raise IndexError("{} {} is not in the columns.\n column list = {}".format(acc_db, pred_col, df.columns))

    # Give the rank of the values
    # NOTE THAT IT IS NECESSARY TO RUN ARGSORT TWICE TO ACHIEVE THIS
    # https://stackoverflow.com/questions/31910407/numpy-argsort-cant-see-whats-wrong

    # rank of experimental values
    df["exp_argsort"] = df[experiment_col].argsort().argsort()
    # rank of prediction values
    df["pred_argsort"] = df[pred_col].argsort().argsort()
    # sort by experimental values, so that iloc can be used below to select increasing sample sizes
    df.sort_values("exp_argsort", inplace=True, ascending=False)

    """Dataframe now contains the rank of experiment and prediction data
    Sorted descending by experiment data.
    
    GpA Elazar 2016 example.
    
       residue_name  interface_score  THOIPA  exp_argsort  pred_argsort
    12            T         0.737715   0.515           14            10
    11            G         0.631558   0.495           13             4
    8             G         0.622738   0.485           12             2
    4             G         0.551554   0.515           11             8
    7             A         0.201012   0.520           10            11
    """

    odf = pd.DataFrame()
    ind = 0
    previous_overlap = 0
    for sample_size in range(1, 11):
        # length of residues tested (usually full TMD length, unless some are excluded)
        #tm_len = len(interface_score_arr)
        non_sample_size = tm_len - sample_size

        # DEPRECATED
        # AM NOT SURE HOW IT WORKED
        #observed_overlap = len(intersect(np.argsort(np.array((-1) * prob_pos))[0:sample_size], np.argsort(np.array((-1) * interface_score_arr))[0:sample_size]))

        """
        get the set of the indices
        
        for the experiental data, this will simply decrease from the length of the TMD
        {19} 
        {18, 19}
        {17, 18, 19} 
        {16, 17, 18, 19}
        
        For the prediction, this will show the rank of the positions highest for the experimental data, e.g.
        {0}
        {0, 3}
        {0, 3, 7}
        {0, 17, 3, 7}
        
        The intersection in this case yield empty sets until a sample size of 4, in which 17 was shared.
        
        NOTE THAT IT IS POSSIBLE TO ADD TWO SHARED POSITIONS BETWEEN EXP. AND PRED. EVEN IF SAMPLE SIZE IS INCREASED BY 1
        """
        exp_set = set(df.exp_argsort.iloc[:sample_size])
        pred_set = set(df.pred_argsort.iloc[:sample_size])
        intersection_result = exp_set.intersection(pred_set)
        observed_overlap = len(intersection_result)

        pval = scipy.special.comb(sample_size, observed_overlap) * scipy.special.comb(non_sample_size, sample_size - observed_overlap) / scipy.special.comb(tm_len, sample_size)
        inter_num_random = 0
        for j in range(1, sample_size + 1):
            random_inter_num_ober = j
            random_non_inter_num = tm_len - sample_size
            inter_num_random = inter_num_random + scipy.special.comb(sample_size, random_inter_num_ober) * scipy.special.comb(random_non_inter_num, sample_size - random_inter_num_ober) / scipy.special.comb(tm_len,
                                                                                                                                                                                                              sample_size) * j
        odf.set_value(ind, acc_db, observed_overlap)
        odf.set_value(ind, "sample_size", "Top" + str(sample_size))
        odf.set_value(ind, "parameters", "observed_overlap")
        ind = ind + 1
        odf.set_value(ind, acc_db, inter_num_random)
        odf.set_value(ind, "sample_size", "Top" + str(sample_size))
        odf.set_value(ind, "parameters", "random_overlap")
        ind = ind + 1
        odf.set_value(ind, acc_db, pval)
        odf.set_value(ind, "sample_size", "Top" + str(sample_size))
        odf.set_value(ind, "parameters", "p_value_from_obs_overlap")
        ind = ind + 1

        previous_overlap = observed_overlap
    odf.set_index(["sample_size", "parameters"], inplace=True, drop=True)
    return odf

def create_one_out_train_data(acc_db,set_path,s):
    df_train = pd.DataFrame()
    df_set04 = pd.read_excel(set_path, sheetname='proteins')
    for j in df_set04.index:
        acc1 = df_set04.loc[j, "acc_db"]
        if not acc1 == acc_db:
            database = df_set04.loc[j, "database"]
            feature_combined_file = os.path.join(s["features_folder"], "combined", database,
                                                 "{}.surr20.gaps5.combined_features.csv".format(acc1))

            df_features_new_protein1 = pd.read_csv(feature_combined_file, index_col=0)
            df_features_new_protein1["acc_db"] = "{}-{}".format(acc1, database)

            # reorder the columns
            df_features_new_protein1 = thoipapy.utils.reorder_dataframe_columns(df_features_new_protein1,
                                                                                ['acc_db', 'residue_num', 'residue_name',
                                                                           'n_homologues'])
            # for the first protein, replace the empty dataframe
            if df_train.empty:
                df_train = df_features_new_protein1
            else:
                # concatenate the growing dataframe of combined proteins and new dataframe
                df_train = pd.concat([df_train, df_features_new_protein1])

                # reset the index to be a range (0,...).
    df_train.index = range(df_train.shape[0])
    return df_train

def get_test_and_train_set_lists(s):

    if s["test_datasets"] == True:
        test_set_list = ["1"]
    elif s["test_datasets"] == False:
        test_set_list = ["0"]
    elif isinstance(s["test_datasets"], int):
        test_set_list = [str(s["test_datasets"])]
    elif isinstance(s["test_datasets"], str):
        test_set_list = s["test_datasets"].split(",")
    else:
        raise ValueError("test_datasets type is not correct {} ({})".format(s["test_datasets"], type(s["test_datasets"])))

    if s["train_datasets"] == True:
        train_set_list = ["1"]
    elif s["train_datasets"] == False:
        train_set_list = ["0"]
    elif isinstance(s["train_datasets"], int):
        train_set_list = [str(s["train_datasets"])]
    elif isinstance(s["train_datasets"], str):
        train_set_list = s["train_datasets"].split(",")
    else:
        raise ValueError("train_datasets type is not correct {} ({})".format(s["train_datasets"], type(s["train_datasets"])))

    return test_set_list, train_set_list

def get_set_lists(s):

    if s["set_list"] == True:
        set_list = ["1"]
    elif s["set_list"] == False:
        set_list = ["0"]
    elif isinstance(s["set_list"], int):
        set_list = [str(s["set_list"])]
    elif isinstance(s["set_list"], str):
        set_list = s["set_list"].split(",")
    else:
        raise ValueError("set_list type is not correct {} ({})".format(s["set_list"], type(s["set_list"])))

    return set_list