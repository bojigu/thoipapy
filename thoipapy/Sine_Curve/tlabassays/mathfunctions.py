
def sigmoid(sigmoid_constants, x):
    '''
    The hill equation is used for the sigmoidal curve. This is standard for EC50 calculations.
    http://en.wikipedia.org/wiki/EC50
    The original exponential function tested was : y = c / (1 + np.exp(-k*(x-x0))) + y0
    Both equations can be used to fit the data quite well.

    NOTE:
    the equation doesn't match the typical hill equation, which is shown in python here:
    http://pydoc.net/Python/cgptoolbox/0.1.2/cgp.sigmoidmodels.doseresponse/

    '''
    #example tuple_of_constants = array([ 0.96782919,  20.67138661, -13.40275013])
    c, k, g = sigmoid_constants
    y = c * 1.0 / (1.0 + ((k/x)**g))
    #y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def sigmoid_fxn_brentq(xvalues_for_curve, sigmoid_constants, y_value_curve_center):
    '''
    This function is for the brentq root finding.
    The only difference to the original sigmoid function is the final line, "- y_value_curve_center"
    Note that if the *args are not supplied, the function will find the sigmoid_constants and y_value_curve_center in the globals.
    '''
    c, k, g = sigmoid_constants
    y = c * 1.0 / (1.0 + ((k/xvalues_for_curve)**g))
    #0.5 is because of the EC50, needs to be changed to calculate an LD25
    #check if this works for non-normalised data??
    return y - y_value_curve_center

def hill_eq_graphPad_prism (hill_constants, x):
    '''
    The hill equation is used for the EC50 calculation, but there are several variants of the formula.
    This equation uses variant 1, as used by GraphPad Prism for EC50 calculations.
    In comparison to the equation on graphpad prism webside, "top" is renamed to "lower" and is generally around zero
    "bottom" is renamed to "upper" and is generally around 1.0
    variant 1: http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    variant 2: http://en.wikipedia.org/wiki/EC50
    variant 3: http://pydoc.net/Python/cgptoolbox/0.1.2/cgp.sigmoidmodels.doseresponse/
    Sigmoid equations:
    variant 4: y = c * 1.0 / (1.0 + ((k/x)**g))
    variant 5: y = c / (1 + np.exp(-k*(x-x0))) + y0

    In theory, any formula that accurately fits a sigmoidal curve is appropriate for our purposes. The EC50 can either
    be calculated directly from the formula (e.g. in Hill equation) or by using brentq and obtaining the x-axis point
    at which the y-axis is at 50%.
    '''
    import numpy as np
    #example tuple_of_constants = array([ 0.96782919,  20.67138661, -13.40275013])
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+10**((np.log10(EC50)-x)-hillslope))
    return y

def hill_eq_graphPad_prism_brentq(xvalues_for_curve, hill_constants, y_value_curve_center):
    '''
    This is the version of the hill_eq used for the brentq root finding.
    The only difference to the original sigmoid function is the final line, "- y_value_curve_center"
    '''
    import numpy as np
    #example tuple_of_constants = array([ 0.96782919,  20.67138661, -13.40275013])
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+10**((np.log10(EC50)-xvalues_for_curve)-hillslope))
    return y - y_value_curve_center

def hill_eq(hill_constants, x):
    '''
    The hill equation is used for the EC50 calculation, but there are several variants of the formula.
    This equation uses variant 2, from wikipedia. Note that because an LD50 is the inverse of a typical EC50 graph.
    In comparison to the equation on wikipedia, "top" is renamed to "lower" and is generally around zero
    "bottom" is renamed to "upper" and is generally around 1.0
    variant 1: http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    variant 2: http://en.wikipedia.org/wiki/EC50
    variant 3: http://pydoc.net/Python/cgptoolbox/0.1.2/cgp.sigmoidmodels.doseresponse/
    Sigmoid equations:
    variant 4: y = c * 1.0 / (1.0 + ((k/x)**g))
    variant 5: y = c / (1 + np.exp(-k*(x-x0))) + y0 (http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator)

    In theory, any formula that accurately fits a sigmoidal curve is appropriate for our purposes. The EC50 can either
    be calculated directly from the formula (e.g. in Hill equation) or by using brentq and obtaining the x-axis point
    at which the y-axis is at 50%.
    '''
    import numpy as np
    #example tuple_of_constants = array([ 0.96782919,  20.67138661, -13.40275013])
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+(x/EC50)**-hillslope)
    return y

def hill_eq_brentq(xvalues_for_curve, hill_constants, y_value_curve_center):
    '''
    This is the version of the hill_eq used for the brentq root finding.
    The only difference to the original sigmoid function is the final line, "- y_value_curve_center"
    In comparison to the equation on wikipedia, "top" is renamed to "lower" and is generally around zero
    "bottom" is renamed to "upper" and is generally around 1.0
    '''
    import numpy as np
    #example tuple_of_constants = array([ 0.96782919,  20.67138661, -13.40275013])
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+(xvalues_for_curve/EC50)**-hillslope)
    return y - y_value_curve_center


def residuals(constants, function, x, y):
    """
    This is a value function (Kostenfunktion) that is used to optimise the fit of the curve to the data.
    It simply calculates the distance between y-value from real data and y-value from the function (sigmoid/sine/etc).
    """
    return y - function(constants, x)

def sine(sine_constants, x):
    ''' Sine equation used to fit sine curve to data.
    f (x) = a*sin (bx + c) + d
    see http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    '''
    import numpy as np
    a, b, c, d = sine_constants
    y = a * np.sin (b * x + c) + d
    return y

def sine_perfect_helix(sine_constants_cd, x):
    ''' Sine equation which is forced to retain alpha helical periodicity, with 3.6 residues per turn.
    f (x) = a*sin (bx + c) + d
    a = 0.2
    b = 1.745
    Why is a fixed to 0.2?
        This is arbitrary, resulting in a curve that is 0.4 in height
    Why is b fixed to 1.745?
        Since periodicity = 2*np.pi/a, for a perfect alpha helix b = 2*np.pi/3.6 = 1.745
    '''
    import numpy as np
    a = 0.2
    b = 1.745
    c, d = sine_constants_cd
    y = a * np.sin (b * x + c) + d
    return y


def resize(arr,lower=0.0,upper=1.0):
    '''
    The resize function normalises the data between the given max and min values.
    The normalisation of the Y-axis in particular seems to help in the curve fitting,
    perhaps because it increases the accuracy of the initial guess of the constants.
    :rtype : object
    ALEX AND I STILL DON'T UNDERSTAND THIS SCRIPT. RECOMMEND USE "NORMALISE_0_1" INSTEAD.
    '''
    #copy the array "arr"
    arr = arr.copy()
    #if someone has mixed up upper and lower, switch the values
    if lower > upper: lower, upper = upper,lower
    #remove the smallest value from the data points
    arr -= arr.min()
    #normalise by multiplying each value in the array by (upper-lower)/arr.max()
    arr *= (upper-lower)/arr.max()
    #add the bottom (lowest value. e.g. 0.3) to all values in the array
    arr += lower
    #return the array
    return arr

def normalise_0_1(arraylike):
    """ Normalises an array to values between 0 and 1.
    The formula for normalisation was norm = (x - array_min)/(array_max - array_min)
    Returns a tuple, consisting of a normalised numpy array of floats, the min of orig data, and the max of the orig data
    """
    import numpy as np
    array_min = np.min(arraylike)
    array_max = np.max(arraylike)
    normalised = (arraylike - array_min)/(array_max - array_min)
    # convert to float
    normalised = np.array(normalised).astype(float)
    return normalised, array_min, array_max

def denormalise_0_1(value_or_array, array_min, array_max):
    """ reverts normalised array back to original values.
    Returns an arraylike if an array is used as input. Returns a single value if a float or integer is used as input.
    The input is the normalised array, and the max and min values in the original array.
    The formula for normalisation was norm = (x - array_min)/(array_max - array_min)
    The formula for denormalisation is therefore x = norm * (array_max - array_min) + array_min
    """
    if isinstance(value_or_array, list):
        raise ValueError('this function accepts arraylike data, not a list. '
                         'Please check data or convert list to numpy array')
    else:
        denormalised = value_or_array*(array_max - array_min) + array_min
    return denormalised


def normalise(arr):
    '''
    Instead of normalising between the highest and lowest values, this function simply divides all values by the highest value (max = 1)
    '''
    arr = arr.copy()
    arr /= arr.max()
    return arr

def test():
    print("congratulations, tlabassays module seems to be installed")