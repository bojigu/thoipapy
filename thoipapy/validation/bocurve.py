import numpy as np
import pandas as pd
from scipy.special import comb


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

    if experiment_col not in df.columns:
        raise IndexError("{} {} is not in the columns.\ntry re-running add_predictions_to_combined_files\ncolumn list = {}".format(acc_db, experiment_col, df.columns))
    if pred_col not in df.columns:
        raise IndexError("{} {} is not in the columns.\ntry re-running add_predictions_to_combined_files\ncolumn list = {}".format(acc_db, pred_col, df.columns))

    # Give the rank of the values
    # NOTE THAT IT IS NECESSARY TO RUN ARGSORT TWICE TO ACHIEVE THIS
    # https://stackoverflow.com/questions/31910407/numpy-argsort-cant-see-whats-wrong

    # drop any nan values in either the experimental column, or the prediction column
    df_sel = df.dropna(subset=[experiment_col, pred_col])

    # get the TMD length from the number of rows with data, after dropping nan
    # as PREDDIMER has longer TMDs, this will result in different statistics regarding the overlaps
    tm_len = df_sel.shape[0]

    # rank of experimental values
    df_sel["exp_argsort"] = df_sel[experiment_col].argsort().argsort()
    # rank of prediction values
    df_sel["pred_argsort"] = df_sel[pred_col].argsort().argsort()
    # sort by experimental values, so that iloc can be used below to select increasing sample sizes
    df_sel.sort_values("exp_argsort", inplace=True, ascending=False)

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
        exp_set = set(df_sel.exp_argsort.iloc[:sample_size])
        pred_set = set(df_sel.pred_argsort.iloc[:sample_size])
        intersection_result = exp_set.intersection(pred_set)
        observed_overlap = len(intersection_result)

        pval = comb(sample_size, observed_overlap) * comb(non_sample_size, sample_size - observed_overlap) / comb(tm_len, sample_size)

        #inter_num_random = calc_rand_overlap(tm_len, sample_size)
        inter_num_random = sample_size**2 / tm_len

        odf.at[ind, acc_db] = observed_overlap
        odf.at[ind, "sample_size"] = "Top" + str(sample_size)
        odf.at[ind, "parameters"] = "observed_overlap"
        ind = ind + 1
        odf.at[ind, acc_db] = inter_num_random
        odf.at[ind, "sample_size"] = "Top" + str(sample_size)
        odf.at[ind, "parameters"] = "random_overlap"
        ind = ind + 1
        odf.at[ind, acc_db] = pval
        odf.at[ind, "sample_size"] = "Top" + str(sample_size)
        odf.at[ind, "parameters"] = "p_value_from_obs_overlap"
        ind = ind + 1

    odf.set_index(["sample_size", "parameters"], inplace=True, drop=True)

    return odf


def parse_BO_data_csv_to_excel(bo_data_csv, bocurve_data_xlsx, logging, predictor_name=""):
    """

    """
    dfb = pd.read_csv(bo_data_csv, index_col=0)

    """ORIGINAL BO DATA CSV LOOKS LIKE THIS
    Top1 = sample size 1
    observed_overlap = overlap in data
    random_overlap = random overlap based on that sequence length and sample size
    p_value_from_obs_overlap = p-value for finding that overlap

                               parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    sample_size                                                                                                                                                                               ...                                                                                                                                       
    Top1                 observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1                   random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1         p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2                 observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2                   random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20
    """

    dfb.index = dfb.index + "_" + dfb.parameters

    """NOW INDICES ARE UNIQUE

                                             parameters  1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    parameters                                                                                                                                                                                                  ...                                                                                                                                       
    Top1_observed_overlap                  observed_overlap        0.000000            0.00            0.00            1.00            0.00        0.000000        0.000000        0.000000        0.000000     ...              1.00     1.000000         0.00     0.000000     0.000000         0.00     0.000000     0.000000     1.000000         0.00
    Top1_random_overlap                      random_overlap        0.035714            0.05            0.05            0.05            0.05        0.045455        0.043478        0.034483        0.043478     ...              0.05     0.041667         0.04     0.047619     0.047619         0.05     0.047619     0.041667     0.047619         0.05
    Top1_p_value_from_obs_overlap  p_value_from_obs_overlap        0.964286            0.95            0.95            0.05            0.95        0.954545        0.956522        0.965517        0.956522     ...              0.05     0.041667         0.96     0.952381     0.952381         0.95     0.952381     0.958333     0.047619         0.95
    Top2_observed_overlap                  observed_overlap        0.000000            0.00            0.00            2.00            0.00        1.000000        0.000000        0.000000        1.000000     ...              1.00     1.000000         0.00     1.000000     0.000000         2.00     2.000000     1.000000     1.000000         2.00
    Top2_random_overlap                      random_overlap        0.142857            0.20            0.20            0.20            0.20        0.181818        0.173913        0.137931        0.173913     ...              0.20     0.166667         0.16     0.190476     0.190476         0.20     0.190476     0.166667     0.190476         0.20

    """
    cols_to_drop = set(["parameters", "Ratio (Average(Ono)/Average(Rno))"]).intersection(set(dfb.columns))
    dfb.drop(cols_to_drop, axis=1, inplace=True)

    # split into separate dataframes. Relabel index to match sample size.
    # dataframe of observed overlaps
    dfobs = dfb[dfb.index.str.contains("observed_overlap")].astype(int)
    dfobs.index = range(1, dfobs.shape[0] + 1)
    # dataframe of random calculated overlaps
    dfrand = dfb[dfb.index.str.contains("random_overlap")]
    dfrand.index = range(1, dfrand.shape[0] + 1)
    # dataframe of p-values
    dfp = dfb[dfb.index.str.contains("p_value_from_obs_overlap")]
    dfp.index = range(1, dfp.shape[0] + 1)

    """FOR EXAMPLE df_obs NOW LOOKS LIKE THIS, WITH A ROW FOR EACH SAMPLE SIZE:

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1               0               0               0               1               0               0               0               0               0               0     ...                 1            1            0            0            0            0            0            0            1            0
    2               0               0               0               2               0               1               0               0               1               0     ...                 1            1            0            1            0            2            2            1            1            2
    3               1               0               0               3               1               1               0               0               1               0     ...                 3            3            0            2            0            2            2            1            1            2
    4               1               0               0               3               1               1               0               0               3               0     ...                 3            4            1            4            0            2            3            2            2            2
    5               1               0               0               4               2               2               2               0               4               0     ...                 3            4            1            4            1            2            3            3            3            2
    """

    # convert the observed overlap (e.g. 4 residues from sample size 10) to a fraction (0.4 of 1)
    dfobs_frac = dfobs.divide(np.array(dfobs.index), axis=0)
    dfrand_frac = dfrand.divide(np.array(dfrand.index), axis=0)

    """If dfobs looks like this:
        ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1               0           0            0
    2               0           0            0
    3               0           0            1
    4               1           0            1
    5               1           1            2
    
    
    Then the dfobs_frac is simply the same data as a fraction (divided by the index)
    
    ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1            0.00         0.0     0.000000
    2            0.00         0.0     0.000000
    3            0.00         0.0     0.333333
    4            0.25         0.0     0.250000
    5            0.20         0.2     0.400000
    """

    # the observed minus random is now calculated for the fraction correct, rather than the original numbers
    df_o_minus_r = dfobs_frac - dfrand_frac
    #df_o_over_r = dfobs_frac / dfrand_frac

    """df_o_minus_r is negative where the result is lower than random

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1       -0.035714           -0.05           -0.05            0.95           -0.05       -0.045455       -0.043478       -0.034483       -0.043478       -0.034483     ...              0.95     0.958333        -0.04    -0.047619    -0.047619        -0.05    -0.047619    -0.041667     0.952381        -0.05
    2       -0.142857           -0.20           -0.20            1.80           -0.20        0.818182       -0.173913       -0.137931        0.826087       -0.137931     ...              0.80     0.833333        -0.16     0.809524    -0.190476         1.80     1.809524     0.833333     0.809524         1.80
    3        0.678571           -0.45           -0.45            2.55            0.55        0.590909       -0.391304       -0.310345        0.608696       -0.310345     ...              2.55     2.625000        -0.36     1.571429    -0.428571         1.55     1.571429     0.625000     0.571429         1.55
    4        0.428571           -0.80           -0.80            2.20            0.20        0.272727       -0.695652       -0.551724        2.304348       -0.551724     ...              2.20     3.333333         0.36     3.238095    -0.761905         1.20     2.238095     1.333333     1.238095         1.20
    5        0.107143           -1.25           -1.25            2.75            0.75        0.863636        0.913043       -0.862069        2.913043       -0.862069     ...              1.75     2.958333         0.00     2.809524    -0.190476         0.75     1.809524     1.958333     1.809524         0.75


    df_o_over_r is where the result is lower than random.
    This is probably not quite a fair comparison, as the zeros are caused by 0 overlap / signigicant random overlap

       1orqC4-crystal  1xioA4-crystal  2axtM1-crystal  2h8aA2-crystal  2j58A1-crystal  2wpdJ1-crystal  3dwwA2-crystal  3h9vA2-crystal  3rifA2-crystal  3spcA2-crystal     ...       Q12913-ETRA  Q12983-ETRA  Q16827-ETRA  Q16832-ETRA  Q6ZRP7-ETRA  Q7L4S7-ETRA  Q8NI60-ETRA  Q92729-ETRA  Q99IB8-ETRA  Q9Y286-ETRA
    1        0.000000             0.0             0.0       20.000000        0.000000        0.000000            0.00             0.0        0.000000             0.0     ...         20.000000        24.00       0.0000     0.000000         0.00     0.000000     0.000000     0.000000    21.000000     0.000000
    2        0.000000             0.0             0.0       10.000000        0.000000        5.500000            0.00             0.0        5.750000             0.0     ...          5.000000         6.00       0.0000     5.250000         0.00    10.000000    10.500000     6.000000     5.250000    10.000000
    3        3.111111             0.0             0.0        6.666667        2.222222        2.444444            0.00             0.0        2.555556             0.0     ...          6.666667         8.00       0.0000     4.666667         0.00     4.444444     4.666667     2.666667     2.333333     4.444444
    4        1.750000             0.0             0.0        3.750000        1.250000        1.375000            0.00             0.0        4.312500             0.0     ...          3.750000         6.00       1.5625     5.250000         0.00     2.500000     3.937500     3.000000     2.625000     2.500000
    5        1.120000             0.0             0.0        3.200000        1.600000        1.760000            1.84             0.0        3.680000             0.0     ...          2.400000         3.84       1.0000     3.360000         0.84     1.600000     2.520000     2.880000     2.520000     1.600000
    
    NOW df_o_minus_r EXPRESSED AS A FRACTION
    It's negative where the random value got a higher score.
    
        ...   3h9vA2-crystal  P05067-NMR  Q12983-ETRA
    1       -0.034483   -0.043478     0.933333
    2        0.431034   -0.086957     0.366667
    3        0.229885   -0.130435     0.466667
    4        0.362069   -0.173913     0.483333
    5        0.227586    0.182609     0.466667
    
    """
    #################################################################
    #          CALCULATE MEAN VALUE FOR ALL TMDS FOR BO-CURVE       #
    #################################################################
    mean_o_minus_r_df = df_o_minus_r.mean(axis=1).to_frame(name="mean_o_minus_r")

    #######################################################################################################
    #                                                                                                     #
    #                CALCULATE AREA UNDER THE BO CURVE FOR SAMPLE SIZES 1-10 (AUBOC10)                    #
    #                                                                                                     #
    #######################################################################################################

    AUBOC10_df = pd.DataFrame()
    for acc_db in df_o_minus_r.columns:
        o_minus_r_ser = df_o_minus_r[acc_db]
        AUBOC10_df.loc[acc_db, "AUBOC10"] = np.trapz(o_minus_r_ser, o_minus_r_ser.index)

        AUBOC10_df.loc[acc_db, "mean_all_sample_sizes"] = o_minus_r_ser.mean()

    mean_AUBOC10 = AUBOC10_df["AUBOC10"].mean()
    logging.info("---{: >24} mean_AUBOC10({:.2f}) n={} ---".format(predictor_name, mean_AUBOC10, AUBOC10_df.shape[0]))
    #################################################################
    #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    #################################################################

    with pd.ExcelWriter(bocurve_data_xlsx) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        dfobs_frac.to_excel(writer, sheet_name="dfobs_frac")
        dfrand_frac.to_excel(writer, sheet_name="dfrand_frac")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        #df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")
        mean_o_minus_r_df.to_excel(writer, sheet_name="mean_o_minus_r")
        AUBOC10_df.to_excel(writer, sheet_name="AUBOC10")