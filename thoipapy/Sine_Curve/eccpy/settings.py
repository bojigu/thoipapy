#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: ECCpy, for EC50 calculation in python
Author: Mark Teese
License: ECCpy is free software, distributed under the GNU Lesser General Public License 3 (LGPLv3)
"""
import os
import pandas as pd
from time import strftime
from thoipapy.Sine_Curve.eccpy.tools import convert_truelike_to_bool, convert_nonelike_to_none

def read_settings_file(settings_excel_file):
    """ Opens the settings excel file tabs as individual dataframes.

    Also creates paths for output files, and adds them to the output dff dataframe.

    Parameters
    ----------
    settings_excel_file : settings file containing the list of datafiles for analysis, and also chosen parameters

    Returns
    -------
    df_settings : DataFrame
        Dataframe containing user settings for EC50 calculation and data analysis.
        Created from the "settings" tab of the settings excel file.
    dff : DataFrame
        Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.
        Contains the "True" / "False" list of input files to analyse in that run.
    shortnames_dict : dict
        Dictionary to convert long sample names to shorter ones that are easier to fit into figures.
        Created from the "shortnames" tab of the settings excel file.
    """
    # convert settings file to pandas dataframe, set the first column "A" as the index
    df_settings = pd.read_excel(settings_excel_file, sheetname = "settings").set_index("A")
    # read the settings file tab that contains a list of short names to describe the data
    df_shortnames = pd.read_excel(settings_excel_file, sheetname="shortnames")
    # convert to a dictionary
    shortnames_dict = dict(zip(df_shortnames.long_name, df_shortnames.short_name))
    # open tab with list of files for analysis as a pandas dataframe (data frame files, dff)
    dff = pd.read_excel(settings_excel_file, sheetname = "files")

    # convert true-like objects (TRUE, true, WAHR, etc) to python bool True
    dff["run curvefit"] = dff["run curvefit"].apply(convert_truelike_to_bool)
    dff["run analysis"] = dff["run analysis"].apply(convert_truelike_to_bool)

    # replace np.nan values with the textstring "None"
    dff["textstring identifier in datafile"] = dff["textstring identifier in datafile"].fillna(value="None")
    dff["textstring identifier in datafile"] = dff["textstring identifier in datafile"].apply(convert_nonelike_to_none)


    # the "output_directory" is optional. replace blank "Not a number, NaN" values with an empty string ""
    dff["output file directory"].fillna("", inplace=True)
    # define Series as the output directory given in the settings file
    ofd = dff.loc[:, "output file directory"]
    # select only empty rows in the output file directory (ofd), obtain index
    ofd_empty_index = ofd.loc[ofd == ""].index
    # replace empty rows with the input file directory
    dff.loc[ofd_empty_index,"output file directory"] = dff.loc[ofd_empty_index,"input file directory"]
    dff.loc[:,"data_file_path"] = dff["input file directory"] + '/' + dff["response data file"]

    # define an output file directory (ofd), normalise the path so that it is os independent
    dff.loc[:,"ofd"] = dff.loc[:,"output file directory"].apply(lambda x: os.path.normpath(x))
    # create the "output_folder" as the directory plus a new folder with the orig response data filename
    dff.loc[:,"data_file_base"] = dff.loc[:,"response data file"].apply(lambda name: name[:-4])
    dff.loc[:,"output_folder"] = dff.loc[:,"ofd"] + "/" + dff.loc[:,"data_file_base"]
    dff.loc[:,"ofd_pdfs"] = dff.loc[:,"output_folder"] + "/" + "pdfs"
    dff.loc[:,"ofd_csv"] = dff.loc[:,"output_folder"] + "/" + "csv"
    dff.loc[:,"ofd_EC50_eval_excel"] = dff.loc[:,"output_folder"] + "/" + dff.loc[:,"data_file_base"] + "_EC50_evaluation.xlsx"
    dff.loc[:,"ofd_EC50_eval_csv"] = dff.loc[:,"ofd_csv"]  + "/" + dff.loc[:,"data_file_base"] + "_EC50_evaluation.csv"
    dff.loc[:,"ofd_EC50_eval_tabsep_csv"] = dff.loc[:,"ofd_csv"]  + "/" + dff.loc[:,"data_file_base"] +  "_EC50_evaluation_tabsep.csv"
    dff.loc[:,"EC50_analysis_fig_basename"] = dff.loc[:,"output_folder"] + "/" + dff.loc[:,"data_file_base"] + "EC50_analysis_fig"
    dff.loc[:,"EC50_analysis_fig_basename_pdf"] = dff.loc[:,"ofd_pdfs"] + "/" + dff.loc[:,"data_file_base"] + "EC50_analysis_fig"

    list_paths_to_normalise = ["data_file_path", "ofd", "output_folder", "ofd_pdfs", "ofd_csv", "ofd_EC50_eval_excel",
                               "ofd_EC50_eval_csv","ofd_EC50_eval_tabsep_csv"]
    # normalise the paths for selected columns, so that they are appropriate for the operating system
    for path in list_paths_to_normalise:
        dff.loc[:,path] = dff.loc[:,path].apply(lambda x: os.path.normpath(x))

    return df_settings, dff, shortnames_dict

def setup_analysed_data_folder(settings_excel_file):
    """ Creates the folder to save the files after analysis.

    Simply creates a subfolder "analysed" in the location of the settings_excel_file.

    The folder is created as follows:
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_datestring/

    The analysed_data_basename for output files includes the datestring again in the filename, and "analysed"
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_datestring/todays_datestring_analysed

    todays_datestring is represented as YEAR|MONTH|DATE, e.g. 20151215.

    Parameters
    ----------
    settings_excel_file : settings file containing the list of datafiles for analysis, and also chosen parameters

    Returns
    -------
    analysed_data_basename : String
        The path for the subfolder, where the data from multiple experiments is saved after analysis.
    """

    """ Create base folder and filename for the output files from analysis of analysed data

    """
    # create a string with the current date
    date_string = strftime("%Y%m%d")

    settings_path, settings_excel_filename = os.path.split(settings_excel_file)

    # create a folder for all the analysed data
    analysed_data_folder = os.path.join(settings_path, "analysed", date_string)

    print(analysed_data_folder)
    if not os.path.exists(analysed_data_folder):
        os.makedirs(analysed_data_folder)
    #print(analysed_data_folder)
    # create an output file for the analysed data from multiple experiments
    analysed_data_basename = os.path.join(analysed_data_folder, date_string + "_analysed")
    return analysed_data_basename
