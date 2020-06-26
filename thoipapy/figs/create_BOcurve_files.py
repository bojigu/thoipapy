import warnings
from pathlib import Path

warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.stats import linregress
#from eccpy.tools import normalise_0_1
from thoipapy.utils import normalise_0_1
warnings.filterwarnings("ignore")


def parse_BO_data_csv_to_excel(bo_data_csv, BO_data_excel, logging, predictor_name=""):
    """

    Run using s["create_AUC_AUBOC_separate_database"]

    Parameters
    ----------
    bo_data_csv
    BO_data_excel
    logging
    predictor_name

    Returns
    -------

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

    with pd.ExcelWriter(BO_data_excel) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        dfobs_frac.to_excel(writer, sheet_name="dfobs_frac")
        dfrand_frac.to_excel(writer, sheet_name="dfrand_frac")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        #df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")
        mean_o_minus_r_df.to_excel(writer, sheet_name="mean_o_minus_r")
        AUBOC10_df.to_excel(writer, sheet_name="AUBOC10")


def parse_BO_data_csv_to_excel_DEPRECATED_NONFRACTION_VERSION(bo_data_csv, BO_data_excel, logging, predictor_name=""):
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

    df_o_minus_r = dfobs - dfrand
    df_o_over_r = dfobs / dfrand

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
    """
    #######################################################################################################
    #                                                                                                     #
    #                CALCULATE AREA UNDER THE BO CURVE FOR SAMPLE SIZES 1-10 (AUBOC10)                    #
    #                                                                                                     #
    #######################################################################################################

    AUBOC10_df = pd.DataFrame()
    for acc_db in df_o_minus_r.columns:
        o_minus_r_ser = df_o_minus_r[acc_db]
        AUBOC10_df.loc[acc_db, "AUBOC10"] = np.trapz(o_minus_r_ser, o_minus_r_ser.index)

    mean_AUBOC10 = AUBOC10_df["AUBOC10"].mean()
    logging.info("---{: >24} mean_AUBOC10({:.2f}) n={} ---".format(predictor_name, mean_AUBOC10, AUBOC10_df.shape[0]))
    #################################################################
    #           SAVE PARSED DATAFRAMES TO AN EXCEL FILE             #
    #################################################################

    with pd.ExcelWriter(BO_data_excel) as writer:
        dfobs.to_excel(writer, sheet_name="dfobs")
        dfrand.to_excel(writer, sheet_name="dfrand")
        dfp.to_excel(writer, sheet_name="dfp")
        df_o_minus_r.to_excel(writer, sheet_name="df_o_minus_r")
        df_o_over_r.to_excel(writer, sheet_name="df_o_over_r")
        AUBOC10_df.to_excel(writer, sheet_name="AUBOC10")



def save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser, plot_o_over_r=False):

    df_o_minus_r = pd.read_excel(BO_data_excel, sheet_name="df_o_minus_r", index_col=0)
    BO_scatter_png = str(BO_barchart_png)[:-12] + "scatter.png"

    #######################################################################################################
    #                                                                                                     #
    #               Create a dataframe with AUBOC10 and AUC for individual protein (df_valid_indiv)       #
    #                                                                                                     #
    #######################################################################################################
    # load AUBOC10 values as a series
    AUBOC10_ser = pd.read_excel(BO_data_excel, sheet_name="AUBOC10", index_col=0)["AUBOC10"].copy()
    # select sample sizes 5 and 10
    df_valid_indiv = df_o_minus_r.loc[[5, 10], :].T.copy()
    df_valid_indiv["AUBOC10"] = AUBOC10_ser
    df_valid_indiv["ROC AUC"] = AUC_ser
    df_valid_indiv.sort_values("AUBOC10", axis=0, ascending=False, inplace=True)

    """ df_valid_indiv should now have the results from BO curve and ROC for each protein
    
                      AUBOC10  sample size 5  sample size 10   ROC AUC
    3ij4_A-crystal  17.456522       1.913043        1.652174  0.714286
    4wit_A-crystal  16.620000       2.000000        2.000000  0.622807
    Q08345-ETRA     16.571429       2.809524        2.238095  0.842593
    P04626-ETRA     16.456522       1.913043        1.652174  0.916667
    P25189-ETRA     14.634615       2.038462        2.153846  0.812500
    """

    #######################################################################################################
    #                                                                                                     #
    #                                plot correlation between AUBOC10 and ROC                             #
    #                                                                                                     #
    #######################################################################################################
    # BO_barchart_png
    plt.close("all")
    #plt.rcParams.update({'font.size': 8})
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    #df_valid_indiv_scatter = df_valid_indiv[["AUBOC10", "ROC AUC"]]
    df_valid_indiv.plot(kind="scatter", ax=ax, x="AUBOC10", y="ROC AUC", alpha=0.7)

    # calculate linear regression for fitted line
    slope, intercept, r_value, p_value, std_err = linregress(df_valid_indiv["AUBOC10"], df_valid_indiv["ROC AUC"])
    #fit_fn = np.poly1d(linear_regression)

    #slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    x_first_last_dp = np.array([df_valid_indiv["AUBOC10"].min(), df_valid_indiv["AUBOC10"].max()])
    y_fitted = x_first_last_dp * slope + intercept
    ax.plot(x_first_last_dp, y_fitted, label="$R^2$ : {:.2f}".format(r_value**2))

    ax.set_xlabel("AUBOC10")
    ax.set_ylabel("ROC AUC")
    ax.legend()
    fig.tight_layout()
    ax.grid(False)
    #BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")

    fig.savefig(BO_scatter_png, dpi=240)

    # simply normalise all between 0 and 1
    for col in df_valid_indiv.columns:
        df_valid_indiv[col] = normalise_0_1(df_valid_indiv[col])[0] + 0.01

    BO_curve_folder = Path(s["thoipapy_data_folder"]) / "Results" / s["setname"] / "crossvalidation"
    crossvalidation_data_dir = BO_curve_folder / "data"
    if not crossvalidation_data_dir.is_dir():
        crossvalidation_data_dir.mkdir(parents=True)
    BO_data_excel = os.path.join(BO_curve_folder, "data", "{}_BO_curve_data.xlsx".format(s["setname"]))

    df_valid_indiv = df_valid_indiv.reindex(columns=["AUBOC10", 5, 10, "ROC AUC"])
    df_valid_indiv.columns = ["AUBOC10", "sample size 5", "sample size 10", "ROC AUC"]

    df_valid_indiv.to_csv(BO_data_excel[:-5] + "_valid_indiv.csv")

    """ df_valid_indiv is now normalised within each column, and sorted by AUBOC10
                          AUBOC10  sample size 5  sample size 10   ROC AUC
    3ij4_A-crystal       1.010000       0.789166        0.727758  0.724139
    4wit_A-crystal       0.980317       0.810587        0.793133  0.594927
    DDR1 [Q08345-ETRA]   0.978593       1.010000        0.837883  0.905371
    ErbB2 [P04626-ETRA]  0.974516       0.789166        0.727758  1.010000
    MPZ [P25189-ETRA]    0.909867       0.820061        0.822048  0.862866
    """

    #######################################################################################################
    #                                                                                                     #
    #                                       plot barchart                                                 #
    #                                                                                                     #
    #######################################################################################################
    # BO_barchart_png
    plt.close("all")
    #plt.rcParams.update({'font.size': 8})
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)
    # replace the protein names
    df_valid_indiv.index = pd.Series(df_valid_indiv.index).replace(namedict)
    df_valid_indiv.plot(kind="bar", ax=ax, alpha=0.7)

    ax.set_ylabel("performance value\n(observed overlap - random overlap)")
    ax.legend()#(["sample size = 5", "sample size = 10"])

    fig.tight_layout()
    ax.grid(False)
    fig.savefig(BO_barchart_png, dpi=240)


    #######################################################################################################
    #                                                                                                     #
    #                                plot linechart (combined data all proteins                           #
    #                                                                                                     #
    #######################################################################################################
    if plot_o_over_r:
        df_o_over_r = pd.read_excel(BO_data_excel, sheet_name="df_o_over_r", index_col=0)
        df_o_over_r_mean = df_o_over_r.T.mean()
    df_o_minus_r.columns = pd.Series(df_o_minus_r.columns).replace(namedict)
    df_o_minus_r_mean = df_o_minus_r.T.mean()
    # get the area under the curve
    AUBOC10 = np.trapz(y=df_o_minus_r_mean, x=df_o_minus_r_mean.index)

    # BO_linechart_png
    plt.close("all")
    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    df_o_minus_r_mean.plot(ax=ax, color="#0f7d9b", linestyle="-", label="prediction (AUBOC10 : {:0.2f}".format(AUBOC10))
    ax.plot([1, 10], [0, 0], color="#0f7d9b", linestyle="--", label="random", alpha=0.5)

    if plot_o_over_r:
        ax2 = ax.twinx()
        df_o_over_r_mean.plot(ax=ax2, color="#9b2d0f", linestyle="-", label="old method (o/r)")
        ax2.plot([1, 10], [1, 1], color="#9b2d0f", linestyle="--", label="old method random", alpha=0.5)

    # ax.set_ylim(0)
    ax.grid(False)
    ax.set_ylabel("fraction of correctly predicted residues\n(observed - random)", color="#0f7d9b")
    ax.tick_params('y', colors="#0f7d9b")

    ax.spines['left'].set_color("#0f7d9b")
    ax.legend()
    if plot_o_over_r:
        ax2.tick_params('y', colors="#9b2d0f")
        ax2.spines['right'].set_color("#9b2d0f")
        #ax.set_ylabel("performance value\n (observed / random)", color="#9b2d0f")
        ax.set_ylabel("fraction of correctly predicted residues\n(observed / random)", color="#9b2d0f")
        ax2.legend()

    ax.set_xlabel("number of TMD residues\n(sample size)")
    fig.tight_layout()
    fig.savefig(BO_linechart_png, dpi=140)

    return AUBOC10


def save_extra_BO_figs(BO_data_excel, other_figs_path):
    linechart_mean_obs_and_rand = os.path.join(other_figs_path, "1_linechart_mean_obs_and_rand.png")
    linechart_obs_indiv = os.path.join(other_figs_path, "2_linechart_obs_indiv.png")
    linechart_p_indiv = os.path.join(other_figs_path, "3_linechart_p_indiv.png")
    linechart_o_minus_r = os.path.join(other_figs_path, "4_linechart_o_minus_r.png")
    linechart_o_over_r = os.path.join(other_figs_path, "5_linechart_o_over_r.png")

    dfrand = pd.read_excel(BO_data_excel, sheet_name="dfrand", index_col=0)
    dfobs = pd.read_excel(BO_data_excel, sheet_name="dfobs", index_col=0)
    df_o_minus_r = pd.read_excel(BO_data_excel, sheet_name="df_o_minus_r", index_col=0)
    # linechart_mean_obs_and_rand

    fig, ax = plt.subplots()
    dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
    dfobs.mean(axis=1).plot(ax=ax, color="k", label="mean observed")
    ax.grid(False)
    ax.set_ylabel("mean overlap")
    ax.legend()
    fig.savefig(linechart_mean_obs_and_rand, dpi=140)

    # linechart_obs_indiv

    plt.close("all")
    fig, ax = plt.subplots()
    dfrand.mean(axis=1).plot(ax=ax, color="k", linestyle="--", label="mean random")
    dfobs.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper left", ncol=2)
    ax.set_ylabel("overlap")
    fig.savefig(linechart_obs_indiv, dpi=140)

    dfp = pd.read_excel(BO_data_excel, sheet_name="dfp", index_col=0)
    # linechart_p_indiv
    plt.close("all")
    fig, ax = plt.subplots()
    dfp.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper right", ncol=2)
    ax.set_ylabel("p-value of result")
    fig.savefig(linechart_p_indiv, dpi=140)

    # linechart_o_minus_r
    plt.close("all")
    fig, ax = plt.subplots()
    df_o_minus_r.plot(ax=ax, alpha=0.7)
    ax.legend(loc="upper left", ncol=2)
    ax.set_ylabel("observed - random")
    fig.savefig(linechart_o_minus_r, dpi=140)
