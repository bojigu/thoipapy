

def judge_EC50_data(dfe,sLet, df_settings):
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
    and "_colour", which indicates the colour of the text on the summary output figure. Black ("k") is used for good
    data, or missing data. Red ("r") is used for data which is judged to be not good enough..
    '''
    import numpy as np
    # setup cutoffs for judging data quality
    # number datapoints neighbouring the EC50 that are excluded from the highdose and lowdose data selection
    # set higher if you use a large number of ampicillin concentrations
    n_neighb = df_settings.loc["n_neighb","B"]
    # maximum standard deviation of the OD600 datapoints at high ampicillin (cells should be dead, variation very small)
    max_std_resp_highdose_dp = df_settings.loc["max_std_resp_highdose_dp","B"]
    # maximum standard deviation of the OD600 datapoints at low ampicillin (final cell density is very variable)
    max_std_resp_lowdose_dp = df_settings.loc["max_std_resp_lowdose_dp","B"]
    # minimum rsquared of the fit from sigmoidal curve to the data
    min_rsquared = df_settings.loc["min_rsquared","B"]
    # minimum acceptable ampicillin concentration stepsizes. Smaller stepsizes give more accurale EC50 values!
    min_acceptable_doseconc_stepsize_at_EC50 = df_settings.loc["min_acceptable_doseconc_stepsize_at_EC50","B"]
    min_recommended_doseconc_stepsize_at_EC50 = df_settings.loc["min_recommended_doseconc_stepsize_at_EC50","B"]
    # minimum hillslope of the fit from sigmoidal curve to the data (below 1, tends not to be sigmoidal)
    min_hillslope = df_settings.loc["min_hillslope","B"]
    # minimum value for the end of the curve, on the y-axis (below -1, tends not to be sigmoidal)
    min_curve_min = df_settings.loc["min_curve_min","B"]

    x_orig = np.array(dfe.loc["x_orig", sLet])
    y_orig = np.array(dfe.loc["y_orig", sLet])
    EC50_orig = dfe.loc["EC50_orig", sLet]
    EC50_ful = dfe.loc["EC50_ful", sLet]

    # create a list that contains the database suffixes (_ful for fixed upper limit, _orig for original)
    datasets = ["_orig", "_ful"]
    for d in datasets:
        # identify the datapoints at high ampicillin concentrations
        dfe.loc["indices_highdose_datapoints{}".format(d),sLet] = np.where(x_orig > dfe.loc["EC50{}".format(d), sLet])[0]
        # remove the datapoint closest to EC50
        dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_highdose_datapoints{}".format(d),sLet][1:]
        # slice using the indices to yield the OD600 values for the highdose datapoints
        dfe.loc["response_highdose_datapoints{}".format(d),sLet] = y_orig[dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
        # count the number of highdose datapoint
        dfe.loc["n_highdose_datapoints{}".format(d),sLet] = len(dfe.loc["response_highdose_datapoints{}".format(d),sLet])

    # identify the lowdose datapoints, count and measure standard deviation
    for d in datasets:
        # identify the lowdose datapoints (x < EC50)
        dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.where(x_orig < dfe.loc["EC50{}".format(d), sLet])[0]
        # exclude datapoint closest to the EC50
        dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_lowdose_datapoints{}".format(d),sLet][:-1]
        # use index to select the y-axis (response) data
        dfe.loc["response_lowdose_datapoints{}".format(d),sLet] = y_orig[dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
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
    for d in datasets:
        if dfe.loc["n_highdose_datapoints{}".format(d),sLet] >= n_neighb:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'

    for d in datasets:
        # evaluate as "okay" if number of highdose or lowdose datapoints is more than two
        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] >= n_neighb:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'

    # judge whether the standard deviation of the high and lowdose datapoints is acceptable
    for d in datasets:
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
            # the cells are overgrown, or have no growth. Stepsize can't be calculated. Replace with 0, and colour black.
            dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'


    for d in datasets:
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
            # the cells are overgrown, or have no growth. Stepsize can't be calculated. Replace with 0, and colour black.
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'

    for d in datasets:
        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] > 0 and dfe.loc["n_highdose_datapoints{}".format(d),sLet] > 0:
            # identify the tested ampicillin concentration below the EC50
            indices_lowdose_datapoints = np.where(x_orig < dfe.loc["EC50{}".format(d),sLet])[0]
            doseconc_before_EC50 = x_orig[indices_lowdose_datapoints[-1]]
            # identify the tested ampicillin concentration after the EC50
            doseconc_after_EC50 = x_orig[dfe.loc["indices_highdose_datapoints{}".format(d),sLet][0]]
            # calculate the stepsize at the EC50. Smaller is better!
            doseconc_stepsize_at_EC50_orig = doseconc_after_EC50 - doseconc_before_EC50
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = doseconc_stepsize_at_EC50_orig
            # evaluate as "okay" if the stepsize at the EC50 is smaller than the min acceptable value
            if doseconc_stepsize_at_EC50_orig <= min_acceptable_doseconc_stepsize_at_EC50:
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = True
                # if the stepsize is small, colour to dark red as a warning that the doseconc should be optimised
                if doseconc_stepsize_at_EC50_orig <= min_recommended_doseconc_stepsize_at_EC50:
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'k'
                else:
                    # colour dark red
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '#990033'
            else:
                # the stepsize is extremely high, and the data therefore has little value. doseconc needs to be optimised.
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # the cells are overgrown, or have no growth. Stepsize can't be calculated. Replace with 0, and colour grey.
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = 0
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '0.5'

    """
    rsquared filter
    """
    for d in datasets:
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
    for d in datasets:
        # evaluate slope as okay if it is higher than a fixed value
        if dfe.loc["hillslope{}".format(d), sLet] > min_hillslope:
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'r'

    # curve_min_norm_orig = dfe.loc["curve_min_norm_orig", sLet]
    """
    curve_min_norm determines if the curve ends at a negative value on the y-axis, suggesting a non-sigmoidal curve
    """
    for d in datasets:
        if dfe.loc["curve_min_norm{}".format(d),sLet] > min_curve_min:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'r'

    return dfe