#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities file containing useful functions.
Authors: Mark Teese, Rimma Jenske
test text
"""

# Fibonacci numbers module. Use this to test that the utilities are working
def fib(n):    # write Fibonacci series up to n
    global a
    a, b = 0, 1
    while b < n:
        print(b),
        a, b = b, a+b

def savefig_if_necessary(savefig, fig, fig_nr, base_filepath, formats = 'png', dpi = 400):
    '''
    Function to save figure with multiple subplots. (i.e., a canvas containing multiple figures)
    Designed to work with the function create_dict_organising_subplots(), which creates a bool object "savefig".
    Automatically names the figure based on the figure number (fig_nr), using a previously defined file path as a base.
    for multiple formats, try formats = ['png','pdf']
    '''
    import matplotlib.pyplot as plt
    if savefig:
        if 'png' in formats:
            fig.savefig(base_filepath + '_%01d.png' % fig_nr, format='png', dpi=dpi)
        if 'pdf' in formats:
            fig.savefig(base_filepath + '_%01d.pdf' % fig_nr, format='pdf')
        #close any open figures
        plt.close('all')

def convert_listlike_cols_to_str(df, list_cols, convert_nan=False, nanreplacement = "[]"):
    ''' Converts arraylike or listlike data within pandas dataframes to stringlists, which can be saved easily in
    csv or excel, and re-converted back to lists or arrays later with the eval function.
    :param df: pandas DataFrame
    :param list_cols: list of columns that contain the arraylike data
    :param convert_nan: if you want the columns without data (np.nan) to also be converted to stringlists "[]"
    :return df: returns the same dataframe, with the converted columns
    '''
    for col in list_cols:
        if col in df:
            # convert each np.array to a list, and then a stringlist
            example_data1 = df[col].dropna()[0]
            # if the datatype is a numpy array or pandas series, convert to a list
            if "ndarray" in str(type(example_data1)) or "Series" in str(type(example_data1)):
                df[col] = df[col].dropna().apply(lambda x: list(x))
            example_data2 = df[col].dropna()[0]
            # check that the first nonnan datapoint is now a list
            if "list" in str(type(example_data2)):
                df[col] = df[col].dropna().apply(lambda x: str(x))
            else:
                raise TypeError("datatype for col {a}, ({b}) is not listlike".format(a=col, b=str(type(example_data2))))
            # if desired, convert np.nan to empty stringlists
            if convert_nan == True:
                df[col] = df[col].fillna(nanreplacement)
        else:
            raise KeyError("The column {} is not in the dataframe".format(col))
    return df


def reindex_df_so_selected_cols_are_first(df, selected_cols):
    """ Takes a dataframe, and rearranges the columns so that the selected columns are placed
    first (on the left)
    Inputs:
    :param df: input dataframe
    :param selected_cols: list of columns to place first
    :return: reindexed dataframe
    """
    # convert columns to list
    col_list_orig = list(df.columns)
    # remove the list_cols_to_place_first from the original columns
    for col in selected_cols:
        if col in col_list_orig:
            col_list_orig.remove(col)
        else:
            raise ValueError("Error, reindex_df_so_selected_cols_are_first, %s not in columns" % col)
    # join to create desired list of columns, and reindex the dataframe
    col_list_final = selected_cols + col_list_orig
    return df.reindex(columns = col_list_final)



class DatafileError(Exception):
    pass

def calc_av_over_3aa_window(df, data_col, new_col, n_aa = 4):
    ''' Function to average the properties of an amino acid (conservation, disruption, etc) across a window, consisting
    of three amino acids on the same side of an alpha helix (i, i-4, i+4)
    INPUTS:
    df: dataframe
    data_col: string showing the column name, which contains the data to be averaged
    new_col: string for the new column name, containing the average
    n_aa: number of amino acids separating the aa. Default = 4 (e.g. i, i-4, i+4)
    OUTPUTS:
    The original dataframe, with an extra column.
    '''
    import numpy as np
    # iterate over each row of the dataframe
    for aa_pos in df.index:
        # define the central datapoint
        disrupt_centre = df.loc[aa_pos,data_col]
        # define the "left" N-term datapoint (e.g. i-4)
        if aa_pos-n_aa in df.index:
            disrupt_left = df.loc[aa_pos-n_aa,data_col]
        else:
            disrupt_left = np.nan
        # define the "right" C-term datapoint (e.g. i+4)
        if aa_pos+n_aa in df.index:
            disrupt_right = df.loc[aa_pos+n_aa,data_col]
        else:
            disrupt_right = np.nan
        # calc average
        av_disrupt = np.nanmean([disrupt_centre, disrupt_left, disrupt_right])
        #print(aa_pos, disrupt_centre, disrupt_left, disrupt_right)
        # add average to dataframe
        df.loc[aa_pos,new_col] = av_disrupt
    return df

def sleep(number, units):
    '''
    sleeps for a certain amount of time
    e.g. sleep(number=15, units="seconds") will sleep for 15 seconds
    units can be either "seconds", "minutes","hours" or "days"
    a single dot "." is printed after each unit of time, whether seconds, minutes or hours
    '''
    import sys
    import time
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    #check that the number is an integer
    if isinstance(number,int):
        if units == "seconds":
            for i in range(number):
                time.sleep(1)
                sys.stdout.write(" .")
                sys.stdout.flush()
        elif units == "minutes":
            for i in range(number):
                time.sleep(60)
                sys.stdout.write(" .")
                sys.stdout.flush()
        elif units == "hours":
            for i in range(number):
                time.sleep(360)
                sys.stdout.write(" .")
                sys.stdout.flush()
        elif units == "days":
            for i in range(number):
                time.sleep(360*24)
                sys.stdout.write(" .")
                sys.stdout.flush()
        else:
            raise ValueError('units can only be seconds, minutes, hours or days')
    else:
        raise ValueError('number must be an integer')
    print(' end.\n')

def create_colour_lists():
    '''
    Converts several lists of rgb colours to the python format (normalized to between 0 and 1)
    Returns a dictionary that contains dictionaries of palettes with named colours (eg. TUM blues)
    and also lists of unnamed colours (e.g. tableau20)
    '''
    import numpy as np
    import matplotlib.colors as colors
    output_dict = {}

    matplotlib_150 = list(colors.cnames.values())
    output_dict['matplotlib_150'] = matplotlib_150

    #define colour dictionaries. TUM colours are based on the style guide.
    colour_dicts = {
                    'TUM_colours' : {
                                    'TUMBlue':(34,99,169),
                                    'TUM1':(100,160,200),
                                    'TUM2':(1,51,89),
                                    'TUM3':(42,110,177),
                                    'TUM4':(153,198,231),
                                    'TUM5':(0,82,147)
                                    },
                    'TUM_accents' : {
                                    'green':(162,183,0),
                                    'orange':(227,114,34),
                                    'ivory':(218,215,203),
                                    }
                    }

    #convert the nested dicts to python 0 to 1 format
    for c_dict in colour_dicts:
        for c in colour_dicts[c_dict]:
            #define r, g, b as ints
            r, g, b = colour_dicts[c_dict][c]
            #normalise r, g, b and add to dict
            colour_dicts[c_dict][c] = (r / 255., g / 255., b / 255.)
        #add normalised colours to output dictionary
        output_dict[c_dict] = colour_dicts[c_dict]

    #define colour lists
    colour_lists = {
                    'tableau20' : [
                                 (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
                                    ],
                    'tableau20blind' : [
                                         (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
                                         (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
                                         (255, 188, 121), (207, 207, 207)
                                          ]
                    }
    #normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list])/255.
        colour_array_tup = tuple(map(tuple,colour_array))
        colour_lists[rgb_list] = colour_array_tup
        #add normalised colours to output dictionary
        output_dict[rgb_list] = colour_lists[rgb_list]
    #create a mixed blue/grey colour list, with greys in decreasing darkness
    TUM_colours_list_with_greys = []
    grey = 0.7
    for c in colour_dicts['TUM_colours'].values():
        TUM_colours_list_with_greys.append('%0.2f' % grey)
        TUM_colours_list_with_greys.append(c)
        grey -= 0.1
    output_dict['TUM_colours_list_with_greys'] = TUM_colours_list_with_greys
    return output_dict