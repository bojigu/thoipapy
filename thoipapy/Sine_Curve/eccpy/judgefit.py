# -*- coding: utf-8 -*-
"""
Package: ECCpy, for EC50 calculation in python
Author: Mark Teese
License: ECCpy is free software, distributed under the GNU Lesser General Public License 3 (LGPLv3)
"""
#import thoipapy.Sine_Curve.eccpy.tools as tools
import thoipapy.Sine_Curve.eccpy.tools
import numpy as np

def judge_fit(dfe, sLet, df_settings):
    '''
    Attempt to determine if the data needs checking.
    Look for the following:
        1) Is the hillslope_orig below 1, indicating an exponential rather than a sigmoid curve?[calculated, used as a filter]
        2) Does the fitted algorithm end below -1?
            - indicates a weak slope, exponential curve rather than a sharp s-curve
        3) How well does the sigmoid curve fit the data (ful, orig)? [calculated, but not currently used as a filter]
        6) How many datapoints before the EC50?[calculated, used as a filter]
        7) What is the y-axis (response) standard deviation of datapoints before the EC50? [calculated, but not currently used as a filter]
        6) How many datapoints after the EC50?[calculated, used as a filter]
        7) What is the y-axis (response) standard deviation of datapoints after the EC50? [calculated, used as a filter]

    The judgements are conducted on two datasets, the fixed upper limit (ful) and the original data (orig).

    The output is the dfe array, with added columns "_okay", referring to whether the EC50 calculation seems to be good,
    data, or missing data. Red ("r") is used for data which is judged to be not good enough..
    '''
    # setup cutoffs for judging data quality
    # number datapoints neighbouring the EC50 that are excluded from the highdose and lowdose data selection
    # set higher if you use a large number of ampicillin concentrations
    n_neighb = df_settings.loc["n_neighb","B"]
    # maximum standard deviation of the OD600 datapoints at high ampicillin (cells should be dead, variation very small)
    max_std_resp_highdose_dp = df_settings.loc["max_std_resp_highdose_dp","B"]
    # maximum standard deviation of the OD600 datapoints at low ampicillin (final cell density is very variable)
    max_std_resp_lowdose_dp = df_settings.loc["max_std_resp_lowdose_dp","B"]
    min_flat_lowdose_dp = df_settings.loc["min_flat_lowdose_dp","B"]
    min_flat_highdose_dp = df_settings.loc["min_flat_highdose_dp","B"]

    # minimum rsquared of the fit from sigmoidal curve to the data
    min_rsquared = df_settings.loc["min_rsquared","B"]
    # minimum acceptable ampicillin concentration stepsizes. Smaller stepsizes give more accurale EC50 values!
    min_acceptable_doseconc_stepsize_at_EC50 = df_settings.loc["min_acceptable_doseconc_stepsize_at_EC50","B"]
    min_recommended_doseconc_stepsize_at_EC50 = df_settings.loc["min_recommended_doseconc_stepsize_at_EC50","B"]
    # minimum hillslope of the fit from sigmoidal curve to the data (below 1, tends not to be sigmoidal)
    weak_hillslope_range = eval(df_settings.loc["weak_hillslope_range","B"])
    # minimum value for the end of the curve, on the y-axis (below -1, tends not to be sigmoidal)
    min_curve_lowresp = df_settings.loc["min_curve_lowresp","B"]

    # create a list that contains the database suffixes (_orig for original, _ful for fixed upper limit)
    # datasets = ["_orig", "_ful"]
    datasets = eval(df_settings.loc["adjust.datasets", "B"])
    for d in datasets:
        x = np.array(dfe.loc["x{}".format(d), sLet])
        y = np.array(dfe.loc["y{}".format(d), sLet])
        # identify the datapoints at high ampicillin concentrations
        dfe.loc["indices_highdose_datapoints{}".format(d),sLet] = np.where(x > dfe.loc["EC50{}".format(d), sLet])[0]
        # remove the datapoint closest to EC50
        dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_highdose_datapoints{}".format(d),sLet][n_neighb:]
        # slice using the indices to yield the OD600 values for the highdose datapoints
        dfe.loc["response_highdose_datapoints{}".format(d),sLet] = y[dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
        # count the number of highdose datapoint
        dfe.loc["n_highdose_datapoints{}".format(d),sLet] = len(dfe.loc["response_highdose_datapoints{}".format(d),sLet])

        # identify the lowdose datapoints, count and measure standard deviation
        # identify the lowdose datapoints (x < EC50)
        dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.where(x < dfe.loc["EC50{}".format(d), sLet])[0]
        # exclude datapoint closest to the EC50
        dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_lowdose_datapoints{}".format(d),sLet][:-n_neighb]
        # use index to select the y-axis (response) data
        dfe.loc["response_lowdose_datapoints{}".format(d),sLet] = y[dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
        # count the datapoints
        dfe.loc["n_lowdose_datapoints{}".format(d),sLet] = len(dfe.loc["response_lowdose_datapoints{}".format(d),sLet])

        # indices_highdose_datapoints_excl_nearest_EC50_orig = indices_highdose_datapoints_orig[1:]
        # response_highdose_datapoints_orig = y_orig[indices_highdose_datapoints_excl_nearest_EC50_orig]
        # # count the number of highdose datapoints
        # dfe.loc["n_highdose_datapoints_orig",sLet] = len(response_highdose_datapoints_orig)

        # # identify the ful datapoints at high ampicillin concentrations, ignoring datapoint closest to EC50
        # indices_highdose_datapoints_ful = np.where(x_orig > EC50_ful)[0]
        # indices_highdose_datapoints_excl_nearest_EC50_ful = indices_highdose_datapoints_ful[1:]
        # response_highdose_datapoints_ful = y_orig[indices_highdose_datapoints_excl_nearest_EC50_ful]
        # # count the number of highdose datapoints
        # dfe.loc["n_highdose_datapoints_ful",sLet] = len(response_highdose_datapoints_ful)

        # judge whether the data contains enough high and lowdose datapoints
        if dfe.loc["n_highdose_datapoints{}".format(d),sLet] >= min_flat_highdose_dp:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'


        # evaluate as "okay" if number of highdose or lowdose datapoints is more than two
        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] >= min_flat_lowdose_dp:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'

        # judge whether the standard deviation of the high and lowdose datapoints is acceptable
        if dfe.loc["n_highdose_datapoints{}".format(d),sLet] > 1:
            # calculate std of highdose datapoints
            dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] = np.std(dfe.loc["response_highdose_datapoints{}".format(d),sLet])
            # evaluate as "okay" if std of highdose datapoints is less than a cutoff value (max_std_resp_highdose_dp)
            if dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] < max_std_resp_highdose_dp:
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_okay" % sLet] = True
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
            else:
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Replace std with 0, and colour black on the figure.
            dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'

        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] > 1:
            # calculate std of lowdose datapoints
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] = np.std(dfe.loc["response_lowdose_datapoints{}".format(d),sLet])
            # evaluate as "okay" if std of lowdose datapoints is less than a cutoff value
            if dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] < max_std_resp_lowdose_dp:
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = True
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
            else:
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Replace std with 0, and colour black.
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'

        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] > 0 and dfe.loc["n_highdose_datapoints{}".format(d),sLet] > 0:
            # identify the tested ampicillin concentration below the EC50
            indices_lowdose_datapoints = np.where(x < dfe.loc["EC50{}".format(d),sLet])[0]
            doseconc_before_EC50 = x[indices_lowdose_datapoints[-1]]
            # identify the tested ampicillin concentration after the EC50
            doseconc_after_EC50 = x[dfe.loc["indices_highdose_datapoints{}".format(d),sLet][0]]
            # add values to output dataframe, so that the plot can be annotated
            dfe.loc["doseconc_steps_at_EC50{}".format(d),sLet] = (doseconc_before_EC50, doseconc_after_EC50)
            # calculate the stepsize at the EC50. Smaller is better!
            doseconc_stepsize_at_EC50 = doseconc_after_EC50 - doseconc_before_EC50
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = doseconc_stepsize_at_EC50
            # evaluate as "okay" if the stepsize at the EC50 is smaller than the min acceptable value
            if doseconc_stepsize_at_EC50 <= min_acceptable_doseconc_stepsize_at_EC50:
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = True
                # if the stepsize is small, colour to dark red as a warning that the doseconc should be optimised
                if doseconc_stepsize_at_EC50 <= min_recommended_doseconc_stepsize_at_EC50:
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'k'
                else:
                    # colour dark red
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '#990033'
            else:
                # the stepsize is extremely high, and the data therefore has little value. doseconc needs to be optimised.
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Stepsize can't be calculated. Replace with 0, and colour grey.
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["doseconc_steps_at_EC50{}".format(d),sLet] = (0,0)
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = 0
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '0.5'

        """
        rsquared filter
        """
        # evaluate rsquared_orig as okay if above a fixed limit (min_rsquared)
        if dfe.loc["rsquared{}".format(d),sLet] > min_rsquared:
            dfe.loc["rsquared{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["rsquared{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["rsquared{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["rsquared{}".format(d),"%s_colour" % sLet] = 'r'

        """
        slope filter used to identify cases where curve does not look sigmoidal, as slope is too weak
        """
        # evaluate slope as okay if it is outside the weak_hillslope_range, which is usually nonsigmoidal
        hillslope = dfe.loc["hillslope{}".format(d), sLet]
        if weak_hillslope_range[0] < hillslope < weak_hillslope_range[1]:
            # if it's outside the range, label as "data_needs_checking"
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # if it's outside the range, label as okay
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'k'

        """
        curve_min_norm determines if the curve ends at a negative value on the y-axis, suggesting a non-sigmoidal curve
        """
        if dfe.loc["curve_min_norm{}".format(d),sLet] > min_curve_lowresp:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'r'

        #######################################################################################################
        #                                                                                                     #
        #                          Calculate the slope at the EC50                                            #
        #                                                                                                     #
        #######################################################################################################

        # calculate the slope surrounding the EC50
        x1 = dfe.loc["EC50_norm_bq{}".format(d),sLet]
        x2 = x1 + 0.01
        hill_constants = dfe.loc["hill_constants{}".format(d),sLet]
        y1 = tools.hill_eq(hill_constants, x1)
        y2 = tools.hill_eq(hill_constants, x2)
        middle_slope = (y2 - y1)/(x2 - x1)
        # print("middle_slope:", middle_slope)

        #######################################################################################################
        #                                                                                                     #
        #                          Calculate the Slope At X-axis Extremes (SAXE)
        #                  1) calculate average stepsize
        #                  2) define points left and right of the xaxis extremes
        #                  3) calculate slope between points
        #                  4) determine if slope is shallow enough to indicate an S- or Z-shaped curve
        #                                                                                                     #
        #######################################################################################################

        xnorm = dfe.loc["xnorm{}".format(d), sLet]
        # find the average stepsize by creating two truncated arrays from xnorm, and subtracting
        xnorm_left = np.array(list(xnorm)[:-1])
        xnorm_right = np.array(list(xnorm)[1:])
        xnorm_stepsizes = xnorm_right - xnorm_left
        xnorm_stepsize_mean = xnorm_stepsizes.mean()
        # define the width surrounding the datapoint for the slope measurement
        # calculated as the mean stepsize multiplied by a user value (0.001 to 1.0)
        width_lowdose_slope = xnorm_stepsize_mean/2 * df_settings.loc["width_lowdose_slope","B"]
        width_highdose_slope = xnorm_stepsize_mean/2 * df_settings.loc["width_highdose_slope","B"]
        # define the min and max of normalised datapoints (will simply be 0 and 1 for normalised data)
        xnorm_min, xnorm_max = xnorm.min(), xnorm.max()
        # define SAXE lowdose/highdose x-axis datapoints (to the left and right of the original datapoints)
        saxe_lowdose_x_dp_left = xnorm_min - width_lowdose_slope
        # if it is negative (results in nan in the sigmoidal function), replace with 0
        saxe_lowdose_x_dp_left = saxe_lowdose_x_dp_left if saxe_lowdose_x_dp_left > 0 else 0
        saxe_lowdose_x_dp_right = xnorm_min + width_lowdose_slope
        saxe_highdose_x_dp_left = xnorm_max - width_highdose_slope
        saxe_highdose_x_dp_right = xnorm_max + width_highdose_slope
        # calculate the y-values on the curve, for the x-values surrounding the min and max datapoints
        saxe_lowdose_y_dp_left = tools.hill_eq(hill_constants, saxe_lowdose_x_dp_left)
        saxe_lowdose_y_dp_right = tools.hill_eq(hill_constants, saxe_lowdose_x_dp_right)
        saxe_highdose_y_dp_left = tools.hill_eq(hill_constants, saxe_highdose_x_dp_left)
        saxe_highdose_y_dp_right = tools.hill_eq(hill_constants, saxe_highdose_x_dp_right)
        # calculate the linear slope (y2 - y1)/(x2 - x1) between the chosen datapoints from the curve
        saxe_lowdose_1 = (saxe_lowdose_y_dp_right - saxe_lowdose_y_dp_left)/(saxe_lowdose_x_dp_right - saxe_lowdose_x_dp_left)
        saxe_highdose_1 = (saxe_highdose_y_dp_right - saxe_highdose_y_dp_left)/(saxe_highdose_x_dp_right - saxe_highdose_x_dp_left)
        # convert slopes to positive numbers
        saxe_lowdose = abs(saxe_lowdose_1)
        saxe_highdose = abs(saxe_highdose_1)
        # add to output dataframe, dfe
        dfe.loc["saxe_lowdose{}".format(d), sLet] = saxe_lowdose
        dfe.loc["saxe_highdose{}".format(d), sLet] = saxe_highdose
        dfe.loc["saxe_lowdose_values{}".format(d), sLet] = [[saxe_lowdose_x_dp_left, saxe_lowdose_x_dp_right],
                                                             [saxe_lowdose_y_dp_left, saxe_lowdose_y_dp_right]]
        dfe.loc["saxe_highdose_values{}".format(d), sLet] = [[saxe_highdose_x_dp_left, saxe_highdose_x_dp_right],
                                                            [saxe_highdose_y_dp_left, saxe_highdose_y_dp_right]]
        # print("saxe_lowdose_x_dp_left", saxe_lowdose_x_dp_left)
        # print("saxe_lowdose_x_dp_right", saxe_lowdose_x_dp_right)
        # print("saxe_highdose_x_dp_left", saxe_highdose_x_dp_left)
        # print("saxe_highdose_x_dp_right", saxe_highdose_x_dp_right)
        # print("saxe_lowdose_y_dp_left", saxe_lowdose_y_dp_left)
        # print("saxe_lowdose_y_dp_right", saxe_lowdose_y_dp_right)
        # print("saxe_highdose_y_dp_left", saxe_highdose_y_dp_left)
        # print("saxe_highdose_y_dp_right", saxe_highdose_y_dp_right)
        # print("\nsaxe_lowdose", saxe_lowdose)
        # print("saxe_highdose", saxe_highdose)

        # obtain the max allowed values for the slopes
        max_lowdose_slope = df_settings.loc["max_lowdose_slope","B"]
        max_highdose_slope = df_settings.loc["max_highdose_slope","B"]

        # check that the calculated slopes of the curve do not exceed the saxe_max_slope
        if saxe_lowdose < max_lowdose_slope:
            dfe.loc["saxe_lowdose{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["saxe_lowdose{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["saxe_lowdose{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["saxe_lowdose{}".format(d),"%s_colour" % sLet] = 'r'

        if saxe_highdose < max_highdose_slope:
            dfe.loc["saxe_highdose{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["saxe_highdose{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["saxe_highdose{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["saxe_highdose{}".format(d),"%s_colour" % sLet] = 'r'

    return dfe