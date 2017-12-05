#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: ECCpy, for EC50 calculation in python
Author: Mark Teese
License: ECCpy is free software, distributed under the GNU Lesser General Public License 3 (LGPLv3)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import sys
from thoipapy.Sine_Curve.eccpy.settings import setup_analysed_data_folder, read_settings_file
from thoipapy.Sine_Curve.eccpy.tools import setup_t20_colour_list, hill_eq, denormalise_0_1
#import warnings
#warnings.simplefilter("error")

def run_analysis(settings_excel_file):
    """ Analysis of EC50 data from multiple experiments.

    Processes the datafiles in settings excel file marked as "TRUE" for "run analysis."
    Collects and analyses the output from the "run_curvefit" program.

    The output of the analysis is saved in the following folder:
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_date/

    Running this script will overwrite any previous files with the same name (i.e., analysed on the same day)

    Parameters
    ----------
    settings_excel_file : settings file containing the list of datafiles for analysis, and also chosen parameters

    Saved Files and Figures
    -------
        EC50_barchart : Barchart of EC50 values from all experiments.
            Only data is included that is judged as good quality (i.e. not labelled as "data_needs_checking").
            Four barcharts are created.
            1) Original data, long sample names
            2) Original data, short sample names
            3) Adjusted data (e.g. fixed upper limit), long sample names
            4) Adjusted data (e.g. fixed upper limit), short sample names
        EC50_datapoints: Scattergram, sample names on x-axis, EC50 on y-axis.
            Effectively the same data as the EC50_barchart, except that the datapoints from each
            experiment are displayed individually, as in a scattergram. A legend with colours for each experiment
            is included. Very useful to determine if the EC50 values from one particular day were uniformly higher or
            lower than the EC50 values calculated on the other days.
            As with the barchart, four variations are created.
            1) Original data, long sample names
            2) Original data, short sample names
            3) Adjusted data (e.g. fixed upper limit), long sample names
            4) Adjusted data (e.g. fixed upper limit), short sample names

    Returns
    -------
    dfc_uniq_summ: DataFrame
        For developers only. The dataframe showing the evaluated LD50 values by sample, for all experiments combined.
    dff : DataFrame
        Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.
        Contains the "True" / "False" list of input files to analyse in that run.

    Note
    -------
    The "run_curvefit" program needs to be run to create the relevent output files, before the analysis program
    can be started. Shifting the location of the output files will result in an error.

    The output files begin as follows:
    todays_date_analysed_
    todays_date is represented as YEAR|MONTH|DAY, e.g. 20151215.

    Note that figures are usually better displayed if a dictionary of long to short sample names is created, and saved
    in the settings_excel_file.

    Currently, the analysis automatically runs for original and adjusted datasets (e.g. fixed upper limit dataset).
    However summary graphs are created separately for each dataset.
    """
    print("\nStarting run_analysis\n")
    # create output folder for analysed data, define basename
    analysed_data_basename = setup_analysed_data_folder(settings_excel_file)
    # add the relevant paths to the data files to the dataframe for files (dff)
    df_settings, dff, shortnames_dict = read_settings_file(settings_excel_file)
    # create t20 colour list
    t20 = setup_t20_colour_list()
    # extract list of adjusted datasets for analysis
    datasets = eval(df_settings.loc["adjust.datasets", "B"])
    """
    COLLECT THE EC50 VALUES FROM ALL THE OUTPUT FILES
    """
    # fix the suffixes denoting the datasets (_orig for original data, _ful for fixed upper limit)
    # if any of the files are labelled True for "run analysis"
    if True in list(dff.loc[:, "run analysis"]):
        # create an empty dataframe to hold the average EC50 values
        dfc = pd.DataFrame()
        # create another dataframe to hold all the boolean data
        dfcb = pd.DataFrame()
        # create a new dataframe for all individual EC50 datapoints, including replicates
        df_allp = pd.DataFrame()
        print("Analysing data from multiple experiments. List of experiments analysed :")
        # iterate through only the files labelled "True" for "run analysis", and join all output dataframes together
        for fn in dff.loc[dff["run analysis"] == True].index:
            # define the response data file
            data_file = dff.loc[fn, "response data file"]
            print(data_file)
            # if it is a real file, open
            if os.path.isfile(dff.loc[fn,"ofd_EC50_eval_excel"]):
                # open  as a new pandas dataframe
                df_eval_values = pd.read_excel(dff.loc[fn,"ofd_EC50_eval_excel"], sheetname="v_" + data_file[:20])
                # add the sample name to the dataframe, so it can be identified later
                df_eval_values["file"] = dff.loc[fn,"ofd_EC50_eval_excel"]
                # convert the sample_name column to a string datatype (forcing np.nan to be "nan")
                df_eval_values["sample_name"] = df_eval_values["sample_name"].astype(str)
                # drop any rows that contain "nan" as the sample name
                df_eval_values = df_eval_values.loc[df_eval_values["sample_name"] != "nan"]
                # join the dataframe with all previously joined dataframes
                dfc = pd.concat([dfc,df_eval_values], axis=0)
                # open the tab of the summary excel file that contains all the boolean values
                df_eval_bool = pd.read_excel(dff.loc[fn,"ofd_EC50_eval_excel"], sheetname="b_" + data_file[:20])
                # join the dataframe with all previously joined dataframes
                dfcb = pd.concat([dfcb,df_eval_bool], axis=0)
                # set the sample_name as the index
                df_eval_values = df_eval_values.set_index("sample_name")
                # iterate through datasets (save in the same df, by adding suffix to the column name)
                for d in datasets:
                    # define the column name in the dataframe
                    col_EC50 = "EC50" + d
                    # select data which "seems okay" according the the automatic data analysis
                    df_eval_values_OK = df_eval_values.loc[df_eval_values["data_seems_okay" + d] == True]
                    # create a list of unique sample names
                    unique_names = list(df_eval_values_OK.index.unique())
                    # create a new dataframe, called df_eval_uniq, which has a single row for each unique sample
                    df_eval_uniq = pd.DataFrame().astype(object)
                    for sn in unique_names:
                        # select only the data for that sample
                        df_sel = df_eval_values_OK.loc[sn,:]
                        # if there is only one sample, the selected data will form a series
                        if isinstance(df_sel, pd.Series):
                            # add the n, the EC50, and the std
                            df_eval_uniq.loc[sn,data_file + d] = df_sel["EC50{}".format(d)]
                        # if the name is not unique, the selected data will form a dataframe.
                        elif isinstance(df_sel, pd.DataFrame):
                            # transfer the EC50 values as a stringlist
                            df_eval_uniq.loc[sn,data_file + d] = str(["%0.2f"%l for l in df_sel[col_EC50]])
                        else:
                            raise TypeError("expected a series or dataframe.")
                    # add the dataframe containing _orig and other dataset columns for that exp to the final df with all data
                    df_allp = pd.concat([df_allp,df_eval_uniq], axis=1)
            else:
                print("File not found! {}".format(dff.loc[fn,"ofd_EC50_eval_excel"]))

        print("\nPercentage data okay:")
        if dfc.empty:
            raise ValueError("No data collected for analysis! Double-check that run_curvefit program has been"
                             " carried out for all the relevant files")

        for d in datasets:
            vc = dfc["data_seems_okay{}".format(d)].value_counts()
            if True in vc:
                n_data_okay = vc[True]
            else:
                n_data_okay = 0
            if False in vc:
                n_data_not_okay = vc[False]
            else:
                n_data_not_okay = 0
            perc_data_okay = n_data_okay / (n_data_okay + n_data_not_okay)*100
            print("{b:0.0f}% ({a} dataset)".format(a=d[1:], b=perc_data_okay))
        # select only the data labeled as "data_seems_okay"
        # save the current index as the sample letter
        dfc["sLet"] = dfc.index
        # convert the index to the sample name
        dfc.index = dfc.sample_name
        # create a new dataframe, called dfc_uniq_summ, which has a single row for each unique sample
        dfc_uniq_summ = pd.DataFrame()
        for d in datasets:
             # select only the data labeled as "data_seems_okay"
            series_data_seems_okay = dfc["data_seems_okay{}".format(d)] == True
            dfc_ok = dfc.loc[series_data_seems_okay]
            # create list of unique sample names, where data is available
            list_unique_sample_names = list(dfc_ok.sample_name.dropna().unique())

            for sn in list_unique_sample_names:
                # select data for that one sample, resulting in either a series, or dataframe (indicating mult. exper.)
                data_1_sample_name = dfc_ok.loc[sn,:]
                # if there is only one sample, df_sel will be a series.
                if isinstance(data_1_sample_name, pd.Series):
                    # add the n, the EC50, and the std
                    dfc_uniq_summ.loc[sn,"n{}".format(d)] = 1
                    dfc_uniq_summ.loc[sn,"mean{}".format(d)] = data_1_sample_name["EC50{}".format(d)]
                    dfc_uniq_summ.loc[sn,"std{}".format(d)] = 0
                # if there are multiple datapoints with the same name, df_sel will be a DataFrame
                elif isinstance(data_1_sample_name, pd.DataFrame):
                    # add the n, the mean EC50 of the samples, and the std
                    dfc_uniq_summ.loc[sn,"n{}".format(d)] = data_1_sample_name["EC50{}".format(d)].shape[0]
                    dfc_uniq_summ.loc[sn,"mean{}".format(d)] = data_1_sample_name["EC50{}".format(d)].mean()
                    dfc_uniq_summ.loc[sn,"std{}".format(d)] = data_1_sample_name["EC50{}".format(d)].std()

        dfc_uniq_summ["longname"] = dfc_uniq_summ.index
        dfc_uniq_summ["shortname"] = dfc_uniq_summ["longname"].apply(lambda x : shortnames_dict[x] if x in list(shortnames_dict.keys()) else x)

        # save the dataframe with all mean data from all experiments to a csv
        dfc_uniq_summ.to_csv(analysed_data_basename + ".csv", sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save both dataframes (mean data and indiv datapoints) from all experiments to excel
        writer = pd.ExcelWriter(analysed_data_basename + ".xlsx")#engine='xlsxwriter'
        dfc_uniq_summ.to_excel(writer, sheet_name = "EC50_mean_all_exp")
        df_allp.to_excel(writer, sheet_name="EC50_indiv_exp")
        writer.save()
        writer.close()

        # # split the path of the settings file to get the filename
        # settings_excel_filename = os.path.split(settings_excel_file)[1]
        # set the fontsize
        fontsize = 6
        plt.rcParams['font.size'] = fontsize
        '''
        AllExperimentsBarcharts_01-04: longnames & shortnames, original and ful
        '''
        # iterate through the original and fixed-upper-limit datasets
        for d in datasets:
            col_mean = "mean" + d
            col_std = "std" + d
            # create a subset to plot that contain data
            df_to_plot = dfc_uniq_summ.dropna(subset = [col_mean]).copy()
            # use the dictionary of long to short names to obtain the short name for the datapoints
            df_to_plot["longname"] = df_to_plot.index
            df_to_plot["shortname"] = df_to_plot.longname.apply(lambda x : shortnames_dict[x] if x in list(shortnames_dict.keys()) else x)
            # iterate through the long and the short names, and produce figures for each
            for namecol in ["longname", "shortname"]:
                sys.stdout.write(".")
                sys.stdout.flush()
                #df_to_plot = df_to_plot.sort_values(by=namecol)
                # define name on the x-axis
                x_names = df_to_plot[namecol]
                # define height of bar or y-axis
                y_data = df_to_plot[col_mean]
                # define error-bars
                yerr = df_to_plot[col_std]
                # define the indices of the boxes on the x-axis
                x_n_boxes = df_to_plot.shape[0]
                box_indices = range(x_n_boxes)
                # define error bar parameters
                error_kw = dict(ecolor='k', lw=1, capsize=2, capthick=1)
                # close all plots, create a new figure, add a barchart
                plt.close("all")
                fig, ax = plt.subplots()
                barcontainer = ax.bar(box_indices, y_data, yerr=yerr, align="center", error_kw=error_kw, color = '#1b9e77')
                # set the xticks
                ax.set_xticks(box_indices)
                # set the labels of the x-axis
                ax.set_xticklabels(x_names, rotation=90)
                # set the limits of the x-axis
                ax.set_xlim([-1, x_n_boxes])
                # set the limit of the y-axis
                ax.set_ylim(0)
                # set the y-axis title
                # ax.set_ylabel("EC50 (ug/ml)")
                ax.set_ylabel("{a}{b}, {c}".format(a=df_settings.loc["calculation_type","B"],
                                                  b=str(df_settings.loc["percentage_response","B"]),
                                                  c=df_settings.loc["x-axis (dose) label","B"]))
                # ax.annotate(s="%s%s" % (namecol,d), xy=(0.015,0.93), fontsize=af, xycoords=xyc)
                ax.set_title("analysed data ({e} experiments),  {a}{b},  {c}".format(a=namecol,b=d,c=os.path.split(settings_excel_file)[1],
                                                                  e=dff.loc[dff["run analysis"] == True].shape[0]))
                # automatically tighten the layout and save figure
                fig.tight_layout()
                # save the figure
                fig.savefig(analysed_data_basename + "_bar_" +  namecol + d + '.png', format='png', dpi=150)
                plt.close('all')

        '''
        AllExperimentsPlots_01-04: Comparison EC50 different experiments. (longnames & shortnames, original and ful)
        '''
        plt.rcParams['legend.numpoints'] = 3
        # sort index
        df_allp.sort_index(inplace=True)
        # copy the index with the full sample names to to a new column
        df_allp["longname"] = df_allp.index
        # create a column with the shortened sample names
        df_allp["shortname"] = df_allp.longname.apply(lambda x : shortnames_dict[x] if x in list(shortnames_dict.keys()) else x)
        # replace index with range
        df_allp.index = range(df_allp.shape[0])
        # give the index a name
        df_allp.index.name = "sSnum"
        # sort the columns
        df_allp.sort_index(axis=1, inplace=True)
        for d in datasets:
            sys.stdout.write(".")
            sys.stdout.flush()

            # identify the colums associated with that dataset
            col_contains_d = list(pd.Series(df_allp.columns).apply(lambda x: d in x))
            # select only data for that dataset (e.g. only orig data)
            sel_df_allp = df_allp.loc[:,col_contains_d]
            # create separate figures for the long names, and the short names
            for name in ["longname", "shortname"]:
                # create a new figure with a single plot
                plt.close("all")
                fig, ax = plt.subplots()
                # set the transparency
                alpha = 0.5
                # define the x-axis labels (long or short)
                xticklabels = list(df_allp[name])
                # iterate through the columns in the dataframe, each representing a single experiment
                for n, c in enumerate(sel_df_allp.columns):
                    # select data within column
                    series = df_allp.loc[:,c].dropna()
                    # create a series of bools, describing if the data is a float
                    isfloat_series = series.apply(lambda x: isinstance(x,float))
                    # use boolean series to select the float data
                    values_single_sample = series.loc[isfloat_series]
                    # create a series of bools, describing if the data is a string (in this case, a stringlist)
                    isstring_series = series.apply(lambda x: isinstance(x,str))
                    if True in list(isstring_series):
                        # use boolean series to select the stringlist data (mult datapoints for same sample in 1 day)
                        values_mult_sample = series.loc[isstring_series]
                        # the multiple EC50 values belong to a single sample number. Identify list of sample numbers.
                        index_ssNum_with_lists = values_mult_sample.index.tolist()
                        # iterate through each selected sample number, associated with multiple EC50 values
                        for sSnum in index_ssNum_with_lists:
                            # convert stringlist to list of floats
                            list_EC50_values = [float(sf) for sf in eval(values_mult_sample[sSnum])]
                            # create corresponding x-axis sample numbers (array of sSnum of length len(values_mult_sample))
                            list_index_values = np.ones(len(list_EC50_values))*sSnum
                            # plot datapoints with same colour as the single-replicate floats
                            ax.scatter(list_index_values, list_EC50_values, color = t20[n], s=40, alpha=alpha, label="_nolegend_")
                    # plot the float data as a scattergram (x-axis is the range, resembling a bar or line chart)
                    ax.scatter(values_single_sample.index, values_single_sample, color=t20[n], s=40, alpha=alpha, label=c[:-5])
                # set the xticks and labels to match the index of df_allp
                ax.set_xticks(np.arange(df_allp.shape[0]))
                ax.set_xticklabels(xticklabels, rotation=90)
                # set the grid to go in between the sample names, as minor xticks
                ax.set_xticks(np.arange(df_allp.shape[0])+0.5, minor=True)
                ax.grid(which='minor', alpha=0.9)
                # set the x axis limits
                ax.set_xlim(-0.5, df_allp.shape[0])
                # set the y-axis title
                ax.set_ylabel("{a}{b}, {c}".format(a=df_settings.loc["calculation_type","B"],
                                                  b=str(df_settings.loc["percentage_response","B"]),
                                                  c=df_settings.loc["x-axis (dose) label","B"]))
                ax.set_title("analysed data ({e} experiments),  {a}{b},  {c}".format(a=name,b=d,c=os.path.split(settings_excel_file)[1],
                                                                  e=sel_df_allp.shape[1]))
                #### FAILED ATTEMPTS TO PREVENT LEGEND CONTAINING EVERYTHING DOUBLED ####
                # set the legend on the top left, with two columns, and a single coloured scatterpoint
                #lgd = ax.legend(loc=2, ncol=2,borderaxespad=0., scatterpoints=1, fontsize=5)
                #lgd = ax.legend(loc='upper left' , ncol=2,scatterpoints=1, fontsize=5)
                # handles, labels = ax.get_legend_handles_labels()
                # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=2, scatterpoints=1)
                #plt.rcParams['legend.handlelength'] = 0
                #ax.subplots_adjust(right=0.7, top=0.5)
                # # change the ylim so that the legend does not cover datapoints
                # ymax = ax.get_ylim()[1]
                # ax.set_ylim(0, ymax*1.25)
                # set legend (currently contains a bug, sometimes all samples there doubled)
                lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=2, scatterpoints=1, numpoints =1)
                # automatically tighten the layout and save figure
                fig.tight_layout()
                fig.savefig(analysed_data_basename + "_datapoints_" + name + d + '.png', format='png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close("all")

        print('\nData analysis is finished.\nOutput files start with the following:\n{}'.format(analysed_data_basename))

        return dfc_uniq_summ, dff
    else:
        print("\nNo files are selected for analysis! Double-check TRUE/FALSE columns in settings file.")
        return "no files selected", "no files selected"

def compare_rawdata(settings_excel_file, list_sample_names, list_tuple_names=None):
    """ Compare raw dose-response curves between selected samples, for selected experiments.

    Processes the datafiles in settings excel file marked as "TRUE" for "run analysis."
    Collects output from the "run_curvefit" program, but only for the selected samples.
    Recreates the fitted curves from the four-parameter Hill equation, with the previously calculated hill_constants.

    The output of the compare_rawdata is saved in the same folder as the "run_analysis":
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_date/

    Running this script will overwrite any previous files with the same name (i.e., analysed on the same day)

    Parameters
    ----------
    settings_excel_file : String
        Path to the settings file containing the list of datafiles for analysis, and also chosen parameters
    list_sample_names : list of strings, or list of tuples
        A list of sample names for analysis. Use the long, original sample names, exactly as
        written in the original data files.
        For multiple pairwise analyses, input a list of tuples containing the relevant sample names.

    Saved Files and Figures
    -------
        Dose_Response_Curve : Scattergram, Line_Chart
            Original dose and response datapoints, plotted as scattergram, with dose on x-axis and response on y-axis.
            Fitted curve, recreated from previously established hill_constants
            1) Original data, long sample names
            2) Adjusted data (e.g. fixed upper limit), long sample names


    Note
    -------
    The compare_rawdata function is best used for the comparison of 2-3 samples, with 5-10 experiments.
    It currently can accept 8 different samples.
    Increasing the number of samples or experiments is likely to result in a very cluttered graph.
    """
    if list_tuple_names != None:
        # check that the list of sample tuples and list of names has the same length
        if len(list_tuple_names) != len(list_sample_names):
            raise IndexError("The length of list_sample_names does not match the "
                             "list_tuple_names. Please check your list of samples.")

    # create output folder for analysed data, define basename
    analysed_data_basename = setup_analysed_data_folder(settings_excel_file)
    compare_raw_data_path = os.path.join(os.path.split(analysed_data_basename)[0], "compare_raw")
    if not os.path.exists(compare_raw_data_path):
        os.mkdir(compare_raw_data_path)
    # fig_raw_data_sel_samples_basename = analysed_data_basename + "_raw_data_sel_samples"
    # # open tab with list of files for analysis as a pandas dataframe (data frame for files, dff)
    # dff = pd.read_excel(settings_excel_file, sheetname = "files")
    t20 = setup_t20_colour_list()
    # add black (k) to the front of the list
    t20.insert(0,"k")
    # add the relevant paths to the data files to the dataframe for files (dff)
    df_settings, dff, shortnames_dict = read_settings_file(settings_excel_file)
    markers = [".",",","o","v","^","<",">","1","2","3","4","8","s","p","*","h","H","+","x","D","d","|","_"]
    # extend the list, in the unlikely case that someone has many replicates!
    markers = markers + markers + markers
    # set transparency of datapoints
    alpha = 0.5
    # set default fontsize
    plt.rcParams['font.size'] = 6
    # define xycoordinates for later annotations
    xyc = "axes fraction"
    # extract list of adjusted datasets for analysis
    datasets = eval(df_settings.loc["adjust.datasets", "B"])

    print("starting compare_curves. Data files:")
    # create boolean
    at_least_one_sample_found_in_selected_datafiles = False

    if not isinstance(list_sample_names[0], tuple):
        # the list of sample names contains only strings, and therefore only a single raw analysis is performed
        # convert to a list containing only one tuple, with the sample names for comparison
        list_sample_names = [tuple(list_sample_names)]

    for stn, sample_tuple in enumerate(list_sample_names):
        print("(Files analysed first 10 characters):")
        if True in list(dff.loc[:, "run analysis"]):
            n_files_to_analyse = dff.loc[dff["run analysis"] == True].shape[0]
            for d in datasets:
                # close any open plots
                plt.close("all")
                # create new canvas (figure) containing a single plot (ax)
                fig, ax = plt.subplots()
                # create a counter for the number of files
                file_counter = 0
                # iterate through all of the data files labeled for analysis
                for fn in dff.loc[dff["run analysis"] == True].index:
                    file_counter += 1
                    data_file = dff.loc[fn, "response data file"]
                    # print the data file name
                    sys.stdout.write("{}, ".format(data_file[:8]))
                    sys.stdout.flush()
                    # open output summary file with LD50 values as pandas dataframe
                    ofd_EC50_eval_csv = dff.loc[fn,"ofd_EC50_eval_csv"]
                    if os.path.isfile(ofd_EC50_eval_csv):
                        filename = os.path.split(ofd_EC50_eval_csv)[1]
                        df = pd.read_csv(ofd_EC50_eval_csv)
                        # set the index as the sample_name (long name)
                        df.set_index("sample_name", inplace=True)
                        # redefine to only include data that is labelled as "data_seems_okay"
                        df = df.loc[df["data_seems_okay{}".format(d)] == True]
                        sample_counter = 0
                        for sample_name in sample_tuple:
                            # counter = 0
                            if sample_name in df.index:
                                at_least_one_sample_found_in_selected_datafiles = True
                                # obtain the bool, or series of bools that say if the data is okay
                                data_seems_okay_X = df.loc[sample_name,"data_seems_okay{}".format(d)]
                                # if it's not a series, the sample name was only found once in that experiment
                                if not isinstance(data_seems_okay_X, pd.Series):
                                    # counter += 1
                                    # convert the x_orig data from a stringlist to a numpy array
                                    x = np.array(eval(df.loc[sample_name,"x{}".format(d)]))
                                    # convert the y_orig data from sample_name stringlist to a numpy array
                                    y = np.array(eval(df.loc[sample_name,"y{}".format(d)]))
                                    # plot the datapoints for that set of data
                                    if sample_counter == 0:
                                        # if it's the first datapoint from that file, set a label for the legend
                                        ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label=filename[:8])
                                    else:
                                        # otherwise, do not write another legend label
                                        ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label="_nolabel_")
                                    # retrieve the hill constants for the curve
                                    hill_constants = eval(df.loc[sample_name,"hill_constants{}".format(d)])
                                    # create 500 datapoints on the x-axis to plot the curve
                                    x_fitted_norm = np.linspace(0, 1, 500)
                                    # create the y datapoints using the sigmoid equation
                                    y_fitted_norm = hill_eq(hill_constants, x_fitted_norm)
                                    # denormalise the x datapoints to the original concentrations
                                    x_fitted = denormalise_0_1(x_fitted_norm, x.min(), x.max())
                                    # denormalise the y datapoints to the original concentrations
                                    y_fitted = denormalise_0_1(y_fitted_norm, y.min(), y.max())
                                    # plot the curve of the fitted data, using the same colours as the datapoints
                                    ax.plot(x_fitted, y_fitted, color = t20[sample_counter], alpha=alpha)
                                    # sample_counter += 1

                                # if it is a series, the sample name was found more than once in that experiment
                                elif isinstance(data_seems_okay_X, pd.Series):
                                    # retrieve the list of x datapoints, y datapoints, and hill constants from curve
                                    x_list_replicates = list(df.loc[sample_name,"x{}".format(d)])
                                    y_list_replicates = list(df.loc[sample_name,"y{}".format(d)])
                                    hill_constants_reps = list(df.loc[sample_name,"hill_constants{}".format(d)])
                                    for i in range(len(x_list_replicates)):
                                        # counter += 1
                                        # convert the x, y and hill constants from a stringlists to numpy arrays
                                        x = np.array(eval(x_list_replicates[i]))
                                        y = np.array(eval(y_list_replicates[i]))
                                        hill_constants = np.array(eval(hill_constants_reps[i]))
                                        # plot the datapoints for that set of data
                                        if sample_counter == 0:
                                            # if it's the first datapoint from that file, set a label for the legend
                                            ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                       marker=markers[file_counter], label=filename[:8])
                                        else:
                                            # otherwise, do not write another legend label
                                            ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                       marker=markers[file_counter], label="_nolabel_")
                                        # retrieve the hill constants for the curve
                                        # print(sample_name, type(df.loc[sample_name,"hill_constants{}".format(d)]), df.loc[sample_name,"hill_constants{}".format(d)])
                                        # hill_constants = eval(df.loc[sample_name,"hill_constants{}".format(d)])
                                        # create 500 datapoints on the x-axis to plot the curve
                                        x_fitted_norm = np.linspace(0, 1, 500)
                                        # create the y datapoints using the sigmoid equation
                                        y_fitted_norm = hill_eq(hill_constants, x_fitted_norm)
                                        # denormalise the x datapoints to the original concentrations
                                        x_fitted = denormalise_0_1(x_fitted_norm, x.min(), x.max())
                                        # denormalise the y datapoints to the original concentrations
                                        y_fitted = denormalise_0_1(y_fitted_norm, y.min(), y.max())
                                        # plot the curve of the fitted data, using the same colours as the datapoints
                                        ax.plot(x_fitted, y_fitted, color = t20[sample_counter], alpha=alpha)
                                else:
                                    raise TypeError("data_seems_okay_X is neither bool nor series")
                            sample_counter += 1
                            # print("{a} : {b}".format(a=sample_name, b=counter))

                if not at_least_one_sample_found_in_selected_datafiles:
                    raise ValueError("No samples found in the selected datasets!\nSamples: {}".format(list_sample_names))
                xaxis_pos = 0.02
                yaxis_pos = np.linspace(0.95,0.7,8)
                for n, sample_name in enumerate(sample_tuple):
                    ax.annotate(s=sample_name,  xy=(xaxis_pos,yaxis_pos[n]),
                                xycoords=xyc,
                                color = t20[n])
                ymax = ax.get_ylim()[1]
                ax.set_ylim(0,ymax*1.3)
                xmax = ax.get_xlim()[1]
                ax.set_xlim(-10,xmax*1.1)
                ax.legend(ncol=2, scatterpoints=1)
                ax.set_title("comparison of raw data for selected samples ({e} experiments),  "
                             "{b},  {c}".format(b=d,c=os.path.split(settings_excel_file)[1],e=n_files_to_analyse))
                # set xlabel, ylabel
                ax.set_xlabel(df_settings.loc["x-axis (dose) label","B"])
                ax.set_ylabel(df_settings.loc["y-axis (response) label","B"],rotation='vertical')
                # set tuple name "t" for saving. If a list of tuple names is not given, use the sample_tuple number, "n"
                t = list_tuple_names[stn] if list_tuple_names != None else stn
                # save the figure in png format
                figpath = os.path.join(compare_raw_data_path, "{n}{d}.png".format(n=t,d=d))
                fig.savefig(figpath, format = "png", dpi = 150)