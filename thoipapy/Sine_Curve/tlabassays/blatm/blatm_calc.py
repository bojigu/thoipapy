# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import brentq
import os
import csv
import tlabassays.mathfunctions as mf
from tlabassays.blatm.judge import judge_EC50_data
from tlabtools.figures import create_dict_organising_subplots
import tlabtools.tools as tools
from time import strftime
import sys
import warnings
# warnings.simplefilter("error")

def calc_all_EC50(settings_excel_file):
    '''
    Calculates EC50 values for the BlaTM assay. Users no longer need to alter anything in the script - it should
    recognise from the data file whether it is a 96-well or 12-well format.
    :param settings_excel_file: excel file with list of data files for analysis
    Output files:
     - a figure showing the datapoints, fitted curve and predicted EC50 from the data.
     - a summary figure with all data from that experiment
     - a summary excel file with the EC50 values
    Example of usage:
    from tlabassays.blatm.blatm_calc import calc_BlaTM_EC50
    settings_excel_file = r"D:\Schweris\Projects\Programming\Python\files\20150920_sBla_96_well_multi_figures\20150920_MT_new_templates_96_well_sBla\sBla_96well_script_settings.xlsx"
    calc_BlaTM_EC50(settings_excel_file)
    '''

    # create output folder for collected data, define basename
    collected_data_basename = setup_collected_data_folder(settings_excel_file)
    # add the relevant paths to the data files to the dataframe for files (dff)
    df_settings, dff, shortnames_dict = read_settings_file(settings_excel_file)
    # create t20 colour list
    t20 = setup_t20_colour_list()

    """
    ANALYSE THE RAW DATA
    """
    # if any of the files are labelled True for "calculate EC50 values"
    if True in list(dff.loc[:, "calculate EC50 values"]):
        for fn in dff.loc[dff["calculate EC50 values"] == True].index:
            # define the location of the microplate reader data file
            data_file = dff.loc[fn, "data file"]
            df_eval_values = calc_EC50(fn, data_file, dff, df_settings, t20)
        return df_eval_values
    else:
        print("None of the datafiles are marked TRUE for 'calculate EC50 values'. Suggest checking the excel settings file.")

def collect_results(settings_excel_file):
    # create output folder for collected data, define basename
    collected_data_basename = setup_collected_data_folder(settings_excel_file)
    # add the relevant paths to the data files to the dataframe for files (dff)
    df_settings, dff, shortnames_dict = read_settings_file(settings_excel_file)
    # create t20 colour list
    t20 = setup_t20_colour_list()

    """
    COLLECT THE EC50 VALUES FROM ALL THE OUTPUT FILES
    """
    # fix the suffixes denoting the datasets (_orig for original data, _ful for fixed upper limit)
    datasets = ["_orig", "_ful"]
    # if any of the files are labelled True for "collect results"
    if True in list(dff.loc[:, "collect results"]):
        # create an empty dataframe to hold the average EC50 values
        dfc = pd.DataFrame()
        # create another dataframe to hold all the boolean data
        dfcb = pd.DataFrame()
        # create a new dataframe for all individual EC50 datapoints, including replicates
        df_allp = pd.DataFrame()
        print("collecting data from multiple experiments:")
        # iterate through only the files labelled "True" for "collect results", and join all output dataframes together
        for fn in dff.loc[dff["collect results"] == True].index:
            # define the data file
            data_file = dff.loc[fn, "data file"]
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
                # iterate through _orig & _ful datasets (save in the same df, by adding suffix to the column name)
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
                    # add the dataframe containing _orig and _ful columns for that exp to the final df with all data
                    df_allp = pd.concat([df_allp,df_eval_uniq], axis=1)

        print("Percentage data okay:")
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
            print("{b:0.0f}% ({a} dataset)".format(a=d, b=perc_data_okay))
        # select only the data labeled as "data_seems_okay"
        dfc = dfc.loc[dfc.data_seems_okay_ful == True]
        # save the current index as the sample letter
        dfc["sLet"] = dfc.index
        # convert the index to the sample name
        dfc.index = dfc.sample_name

        # create a list of unique sample names from all the data
        list_unique_sample_names = list(dfc.sample_name.dropna().unique())

        # create a new dataframe, called dfc_uniq_summ, which has a single row for each unique sample
        dfc_uniq_summ = pd.DataFrame()
        for sn in list_unique_sample_names:
            # create a selected dataframe containing only the data for that sample
            df_sel = dfc.loc[sn,:]
            for d in datasets:
                # if there is only one sample, df_sel will be a series.
                if isinstance(df_sel, pd.Series):
                    # add the n, the EC50, and the std
                    dfc_uniq_summ.loc[sn,"n{}".format(d)] = 1
                    dfc_uniq_summ.loc[sn,"mean{}".format(d)] = df_sel["EC50{}".format(d)]
                    dfc_uniq_summ.loc[sn,"std{}".format(d)] = 0
                # if there are multiple samples with the same name, df_sel will be a DataFrame
                elif isinstance(df_sel, pd.DataFrame):
                    # there should be multiple samples with the same name, returning a DataFrame
                    # add the n, the mean EC50, and the std
                    dfc_uniq_summ.loc[sn,"n{}".format(d)] = df_sel.shape[0]
                    dfc_uniq_summ.loc[sn,"mean{}".format(d)] = df_sel["EC50{}".format(d)].mean()
                    dfc_uniq_summ.loc[sn,"std{}".format(d)] = df_sel["EC50{}".format(d)].std()

        dfc_uniq_summ["longname"] = dfc_uniq_summ.index
        dfc_uniq_summ["shortname"] = dfc_uniq_summ["longname"].apply(lambda x : shortnames_dict[x] if x in list(shortnames_dict.keys()) else x)

        # save the dataframe with all mean data from all experiments to a csv
        dfc_uniq_summ.to_csv(collected_data_basename + ".csv", sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save both dataframes (mean data and indiv datapoints) from all experiments to excel
        writer = pd.ExcelWriter(collected_data_basename + ".xlsx")#engine='xlsxwriter'
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
            df_to_plot = dfc_uniq_summ.dropna(subset = [col_mean])
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
                                                  c=df_settings.loc["x-axis (dose)","B"]))
                # ax.annotate(s="%s%s" % (namecol,d), xy = (0.015,0.93), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction')
                ax.set_title("collected data ({e} experiments),  {a}{b},  {c}".format(a=namecol,b=d,c=os.path.split(settings_excel_file)[1],
                                                                  e=dff.loc[dff["collect results"] == True].shape[0]))
                # automatically tighten the layout and save figure
                fig.tight_layout()
                # save the figure
                fig.savefig(collected_data_basename + "_bar_" +  namecol + d + '.png', format='png', dpi=150)
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
            # create a new figure with a single plot
            fig, ax = plt.subplots()
            # set the transparency
            alpha = 0.5
            # identify the colums associated with that dataset
            col_contains_d = list(pd.Series(df_allp.columns).apply(lambda x: d in x))
            # select only data for that dataset (e.g. only orig data)
            sel_df_allp = df_allp.loc[:,col_contains_d]
            # create separate figures for the long names, and the short names
            for name in ["longname", "shortname"]:
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
                    ax.scatter(values_single_sample.index, values_single_sample, color=t20[n], s=40, alpha=alpha, label=c)
                # set the xticks and labels to match the index of df_allp
                ax.set_xticks(np.arange(df_allp.shape[0]))
                ax.set_xticklabels(xticklabels, rotation=90)
                # set the grid to go in between the sample names, as minor xticks
                ax.set_xticks(np.arange(df_allp.shape[0])+0.5, minor=True)
                ax.grid(which='minor', alpha=0.9)
                # set the x axis limits
                ax.set_xlim(-0.5, df_allp.shape[0])
                # set the y-axis title
                # ax.set_ylabel("EC50 (ug/ml)")
                ax.set_ylabel("{a}{b}, {c}".format(a=df_settings.loc["calculation_type","B"],
                                                  b=str(df_settings.loc["percentage_response","B"]),
                                                  c=df_settings.loc["x-axis (dose)","B"]))
                # ax.annotate(s="%s%s" % (namecol,d), xy = (0.015,0.93), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction')
                ax.set_title("collected data ({e} experiments),  {a}{b},  {c}".format(a=name,b=d,c=os.path.split(settings_excel_file)[1],
                                                                  e=sel_df_allp.shape[1]))
                # set the legend on the top left, with two columns, and a single coloured scatterpoint
                #lgd = ax.legend(loc=2, ncol=2,borderaxespad=0., scatterpoints=1, fontsize=5)
                #lgd = ax.legend(loc='upper left' , ncol=2,scatterpoints=1, fontsize=5)
                # handles, labels = ax.get_legend_handles_labels()
                # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=2, scatterpoints=1)
                plt.rcParams['legend.handlelength'] = 0
                lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=2, scatterpoints=1, numpoints =1)
                #ax.subplots_adjust(right=0.7, top=0.5)
                # # change the ylim so that the legend does not cover datapoints
                # ymax = ax.get_ylim()[1]
                # ax.set_ylim(0, ymax*1.25)
                # automatically tighten the layout and save figure
                fig.tight_layout()
                fig.savefig(collected_data_basename + "_datapoints_" + name + d + '.png', format='png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close("all")


        # # for Tanja's dataset: calculate the disruption index at particular positions
        # # choose first amino acid in graph
        # start_aa = 233
        # # choose last amino acid to appear in graph
        # end_aa = 251
        #
        # for i in range(start_aa,end_aa+1):
        #     # for each sample name in the dataframe with unique samples
        #     for sn in dfc_uniq_summ.index:
        #         # if the amino acid number (233, 234, etc) is in the sample name
        #         if str(i) in sn:
        #             # add the amino acid number to a new column
        #             dfc_uniq_summ.loc[sn,"pos"] = i
        # # sort the samples by the amino acid position in the dataframe
        # dfc_uniq_summ = dfc_uniq_summ.sort_values(by="pos")
        # # copy the current index to a new column, "sample_name"
        # dfc_uniq_summ["sample_name"] = dfc_uniq_summ.index
        # # create a regex string for searching all the sample names for the particular amino acid numbers
        # dfc_uniq_summ["pos_re"] = dfc_uniq_summ.pos.apply(lambda x: "\w%s\w" % x)

        # # create a list of unique amino acid positions (for which there could be multiple mutations)
        # list_unique_pos = dfc_uniq_summ.pos.dropna().unique()
        # dfc_uniq_summ_posindex = dfc_uniq_summ.copy()
        # dfc_uniq_summ_posindex.index = dfc_uniq_summ_posindex.pos
        # dfc_mult = pd.DataFrame()
        # for pos in list_unique_pos:
        #     #convert to integer
        #     pos = int(pos)
        #     # select the rows with that position
        #     df_sum_sel = dfc_uniq_summ_posindex.loc[pos,:]
        #     # calculate mean and std
        #     if "mean" in df_sum_sel.index:
        #         dfc_mult.loc[pos,"mean_mult"] = df_sum_sel["mean"]
        #         dfc_mult.loc[pos,"std_mult"] = 0
        #     else:
        #         dfc_mult.loc[pos,"mean_mult"] = df_sum_sel["mean"].mean()
        #         dfc_mult.loc[pos,"std_mult"] = df_sum_sel["mean"].std()
        #
        # list_long_names = ["N(P01920_B2_wt)_C(P01909_A1_wt)_sBla1.1","N(P02724_0_wt)_C(P02724_0_wt)_sBla1.1"]
        # list_short_names = ["wt","GpA"]
        # df_select = dfc_uniq_summ.loc[list_long_names,:]
        # df_select.index = list_short_names
        # df_select.columns = ['n', 'mean_mult', 'std_mult', 'pos', 'sample_name', 'pos_re']
        # df_select = df_select.reindex(columns = ['mean_mult', 'std_mult'])
        # dfc_mult = pd.concat([dfc_mult,df_select])
        #
        # #outfig = r"D:\Schweris\Projects\Programming\Python\files\20150920_sBla_96_well_multi_figures\20150920_MT_new_templates_96_well_sBla\quick_summary_mutations.png"
        # yerr = dfc_mult["std_mult"]
        # dfc_mult.mean_mult.plot(kind="bar", color = "#0489B1", yerr = [yerr,yerr])
        # plt.savefig(collected_data_basename + "scanmut.png", dpi=150)
        #
        # from tlabassays.sin_disrupt_fit import fit_sin_to_data
        # for name in list_short_names:
        #     dfc_mult.drop(name, inplace=True)
        # x_sin = dfc_mult.index
        # y_sin = mf.normalise_0_1(dfc_mult.mean_mult)[0]
        # title = "HLA II beta scanmut"
        # #basename = r"D:\Schweris\Projects\Programming\Python\files\20150920_sBla_96_well_multi_figures\20150920_MT_new_templates_96_well_sBla\quick_summary"
        # # HLA_sine_constants, HLA_periodicity = fit_sin_to_data(x_sin,y_sin, title, collected_data_basename + "scanmut_sine")
        print('\n Data analysis is finished')

        return dfc_uniq_summ, dff
    else:
        raise ValueError ("No files are selected for analysis! Double-check TRUE/FALSE columns in settings file.")

def compare_raw_data_selected_samples(settings_excel_file, list_sample_names):
    # # slice the path of the settings file to get the filename
    # settings_excel_filename = os.path.split(settings_excel_file)[1]
    # # create output folder for collected data, define basename
    collected_data_basename = setup_collected_data_folder(settings_excel_file)
    # fig_raw_data_sel_samples_basename = collected_data_basename + "_raw_data_sel_samples"
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
    datasets = ["_orig", "_ful"]
    print("starting compare_raw_data_selected_samples. Data files:")
    if True in list(dff.loc[:, "collect results"]):
        n_files_to_analyse = dff.loc[dff["collect results"] == True].shape[0]
        for d in datasets:
            plt.close("all")
            fig, ax = plt.subplots()
            for fn in dff.loc[dff["collect results"] == True].index:
                data_file = dff.loc[fn, "data file"]
                sys.stdout.write("{}, ".format(data_file[:8]))
                sys.stdout.flush()
                ofd_EC50_eval_csv = dff.loc[fn,"ofd_EC50_eval_csv"]
                if os.path.isfile(ofd_EC50_eval_csv):
                    filename = os.path.split(ofd_EC50_eval_csv)[1]
                    df = pd.read_csv(ofd_EC50_eval_csv)
                    df.set_index("sample_name", inplace=True)
                    sample_counter = 0
                    for sample_name in list_sample_names:
                        if sample_name in df.index:
                            if not isinstance(df.loc[sample_name,"data_seems_okay{}".format(d)], pd.Series):
                                if df.loc[sample_name,"data_seems_okay{}".format(d)] == True:
                                    # convert the x_orig data from a stringlist to a numpy array
                                    x_orig = np.array(eval(df.loc[sample_name,"x_orig"]))
                                    # convert the y_orig data from sample_name stringlist to a numpy array
                                    y_orig = np.array(eval(df.loc[sample_name,"y{}".format(d)]))
                                    # plot the datapoints for that set of data
                                    if sample_counter == 0:
                                        # if it's the first datapoint from that file, set a label for the legend
                                        ax.scatter(x_orig, y_orig, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label=filename[:8])
                                    else:
                                        # otherwise, do not write another legend label
                                        ax.scatter(x_orig, y_orig, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label="_nolabel_")
                                    # retrieve the hill constants for the curve
                                    hill_constants = eval(df.loc[sample_name,"hill_constants{}".format(d)])
                                    # create 500 datapoints on the x-axis to plot the curve
                                    x_fitted_norm = np.linspace(0, 1, 500)
                                    # create the y datapoints using the sigmoid equation
                                    y_fitted_norm = mf.hill_eq(hill_constants, x_fitted_norm)
                                    # denormalise the x datapoints to the original concentrations
                                    x_fitted_non_normalised = mf.denormalise_0_1(x_fitted_norm, x_orig.min(), x_orig.max())
                                    # denormalise the y datapoints to the original concentrations
                                    yvalues_for_curve_non_normalised = mf.denormalise_0_1(y_fitted_norm, y_orig.min(), y_orig.max())
                                    # plot the curve of the fitted data, using the same colours as the datapoints
                                    ax.plot(x_fitted_non_normalised, yvalues_for_curve_non_normalised, color = t20[sample_counter], alpha=alpha)
                            elif isinstance(df.loc[sample_name,"data_seems_okay{}".format(d)], pd.Series) and True in df.loc[sample_name,"data_seems_okay{}".format(d)]:
                                x_orig_list_replicates = list(df.loc[sample_name,"x_orig"])
                                y_orig_list_replicates = list(df.loc[sample_name,"y{}".format(d)])
                                for i in range(len(x_orig_list_replicates)):
                                    # convert the x_orig data from a stringlist to a numpy array
                                    x_orig = np.array(eval(x_orig_list_replicates[i]))
                                    # convert the y_orig data from sample_name stringlist to a numpy array
                                    y_orig = np.array(eval(y_orig_list_replicates[i]))
                                    # plot the datapoints for that set of data
                                    if sample_counter == 0:
                                        # if it's the first datapoint from that file, set a label for the legend
                                        ax.scatter(x_orig, y_orig, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label=filename[:8])
                                    else:
                                        # otherwise, do not write another legend label
                                        ax.scatter(x_orig, y_orig, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[fn], label="_nolabel_")
                                    # retrieve the hill constants for the curve
                                    hill_constants = eval(df.loc[sample_name,"hill_constants{}".format(d)])
                                    # create 500 datapoints on the x-axis to plot the curve
                                    x_fitted_norm = np.linspace(0, 1, 500)
                                    # create the y datapoints using the sigmoid equation
                                    y_fitted_norm = mf.hill_eq(hill_constants, x_fitted_norm)
                                    # denormalise the x datapoints to the original concentrations
                                    x_fitted_non_normalised = mf.denormalise_0_1(x_fitted_norm, x_orig.min(), x_orig.max())
                                    # denormalise the y datapoints to the original concentrations
                                    yvalues_for_curve_non_normalised = mf.denormalise_0_1(y_fitted_norm, y_orig.min(), y_orig.max())
                                    # plot the curve of the fitted data, using the same colours as the datapoints
                                    ax.plot(x_fitted_non_normalised, yvalues_for_curve_non_normalised, color = t20[sample_counter], alpha=alpha)
                            sample_counter = sample_counter + 1
            xaxis_pos = 0.02
            yaxis_pos = np.linspace(0.95,0.7,8)
            for n, sample_name in enumerate(list_sample_names):
                ax.annotate(s=sample_name,  xy = (xaxis_pos,yaxis_pos[n]),
                            xytext=None, xycoords='axes fraction',
                            color = t20[n])
            ymax = ax.get_ylim()[1]
            ax.set_ylim(0,ymax*1.3)
            xmax = ax.get_xlim()[1]
            ax.set_xlim(-10,xmax*1.1)
            ax.legend(ncol=2, scatterpoints=1)
            ax.set_title("comparison of raw data for selected samples ({e} experiments),  {b},  {c}".format(b=d,c=os.path.split(settings_excel_file)[1],
                                                                               e=n_files_to_analyse))
            # set xlabel, ylabel
            # axarr[row_nr].set_xlabel('ampicillin concentration (µg/ml)', fontsize = fig_fontsize)
            ax.set_xlabel(df_settings.loc["x-axis (dose)","B"])
            # axarr[row_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
            ax.set_ylabel(df_settings.loc["y-axis (response)","B"],rotation='vertical')
            # # set the y-axis title
            # ax.set_xlabel("ampicillin concentration (ug/ml)")
            # ax.set_ylabel("cell density (A600)")
            fig.savefig(collected_data_basename + "_raw_data_sel_samples" + d + ".png", format = "png", dpi = 300)

def setup_collected_data_folder(settings_excel_file):
    """ Create base folder and filename for the output files from analysis of collected data

    """
    # create a string with the current date
    date_string = strftime("%Y%m%d")

    settings_path, settings_excel_filename = os.path.split(settings_excel_file)

    # create a folder for all the collected data
    collected_data_folder = os.path.join(settings_path, "collected", date_string)

    print(collected_data_folder)
    if not os.path.exists(collected_data_folder):
        os.makedirs(collected_data_folder)
    #print(collected_data_folder)
    # create an output file for the collected data from multiple experiments
    collected_data_basename = os.path.join(collected_data_folder, date_string + "_collected")
    return collected_data_basename

def setup_t20_colour_list():
    # create a colour list with 30 colours for the graphs, based on tableau20
    colour_lists = tools.create_colour_lists()
    # extract the tableau20 lists, join together
    t20 = list(colour_lists["tableau20"] + colour_lists["tableau20blind"])
    # extend the list(colours will be redundant for large numbers of samples/experiments!)
    t20 = t20 + t20
    return t20

def read_settings_file(settings_excel_file):
    """ Opens the settings excel file tabs as individual dataframes. Creates paths for output files.

    :param settings_excel_file: user input excel file, containing all the files for analysis and chosen parameters

    :return df_settings, dff, shortnames_dict
    df_settings: dataframe for settings, contains parameters used to analyse all files
    dff: dataframe for files, contains all the paths for input and output files
    shortnames_dict: dictionary to convert long sample names to shorter ones that are easier to fit into figures
    """
    # convert settings file to pandas dataframe, set the first column "A" as the index
    df_settings = pd.read_excel(settings_excel_file, sheetname = "settings").set_index("A")
    # read the settings file tab that contains a list of short names to describe the data
    pd_shortnames = pd.read_excel(settings_excel_file, sheetname="shortnames")
    # convert to a dictionary
    shortnames_dict = dict(zip(pd_shortnames.long_name, pd_shortnames.short_name))
    # open tab with list of files for analysis as a pandas dataframe (data frame files, dff)
    dff = pd.read_excel(settings_excel_file, sheetname = "files")
    # the "output_directory" is optional. replace blank "Not a number, NaN" values with an empty string ""
    dff["output file directory (leave blank to use input directory)"].fillna("", inplace=True)
    # define Series as the output directory given in the settings file
    ofd = dff.loc[:, "output file directory (leave blank to use input directory)"]

    # select only empty rows in the output file directory (ofd), obtain index
    ofd_empty_index = ofd.loc[ofd == ""].index
    # replace empty rows with the input file directory
    dff.loc[ofd_empty_index,"output file directory (leave blank to use input directory)"] = dff.loc[ofd_empty_index,"input file directory"]
    dff.loc[:,"data_file_path"] = dff["input file directory"] + '/' + dff["data file"]

    # define an output file directory (ofd), normalise the path so that it is os independent
    dff.loc[:,"ofd"] = dff.loc[:,"output file directory (leave blank to use input directory)"].apply(lambda x: os.path.normpath(x))
    # create the "output_folder" as the directory plus a new folder with the orig data filename
    dff.loc[:,"data_file_base"] = dff.loc[:,"data file"].apply(lambda name: name[:-4])
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

def calc_EC50(fn, data_file, dff, df_settings, t20):
    # define the datasets (i.e. data adjusted before fitting, such as the "fixed upper limit" dataset)
    datasets = ["_orig", "_ful"]

    # create new output file directories, if they don't exist already
    list_paths = [dff.loc[fn, "output_folder"], dff.loc[fn, "ofd_pdfs"], dff.loc[fn, "ofd_csv"]]
    for path in list_paths:
        if not os.path.exists(path):
            os.makedirs(path)

    assay_type, input_data_file_seems_okay, sample_96_well, \
    sample_12_well_plate = examine_input_datafile(fn, data_file, dff.loc[fn,"data_file_path"], df_settings.loc["identifier_VERSAmax", "B"],
                                                   eval(df_settings.loc["list_96_well_doseconc_types", "B"]))

    # for the 96-well samples, obtain path for file with the ampicillin concentrations and sample names
    amp_conc_excelfile = dff.loc[fn, "file with concentrations and sample names"]
    # replace np.nan with an empty string, if the sample is from the 12-well platereader and no excel file is given
    if isinstance(amp_conc_excelfile, float):
        if np.isnan(amp_conc_excelfile):
            amp_conc_excelfile = ""
    amp_conc_excel_path = os.path.join(dff.loc[fn, "input file directory"], amp_conc_excelfile)


    if sample_96_well == True:
        df_resp_orig, df_resp_all, df_dose_orig = read_96well_inputfiles(assay_type, input_data_file_seems_okay,
                                                                dff.loc[fn,"data_file_path"], amp_conc_excel_path)
    else:
        df_dose_orig = "not yet created"
        df_resp_all = "not yet created"
        # do nothing. the data file is not from VersaMax, is probably from the 12-well platereader
        pass

    dfAC_all, df_resp_all = standardise_doseconc_data(assay_type, df_dose_orig, df_resp_all, dff.loc[fn,"data_file_path"])

    if sample_96_well == True:
        # the y-value (response) for A01 will be placed in a new dataframe at df_resp_orig.loc["A","01"]
        # identify index for datapoints
        df_resp_orig['index'] = df_resp_orig.Sample.apply(lambda xa : xa[0])
        # identify columns for datapoints
        df_resp_orig['column'] = df_resp_orig.Sample.apply(lambda xb : xb[1:3])
        # create empyty dataframe to hold the data
        df_resp_all = pd.DataFrame(index = dfAC_all.index, columns = dfAC_all.columns)
        # add the "Contains_Data" column
        df_resp_all["Contains_Data"] = dfAC_all["Contains_Data"]

        # iterate through the rows, dispersing the data throughout the new dataframe df_resp_all
        for row in df_resp_orig.index:
            index_in_df_resp_all = df_resp_orig.loc[row, "index"]
            column_in_df_resp_all = df_resp_orig.loc[row, "column"]
            y_value_response = df_resp_orig.loc[row, "MeanOD600"]
            #transfer data point
            df_resp_all.loc[index_in_df_resp_all,column_in_df_resp_all] = y_value_response

        # create new dataframe from excel file containing Sample names
        dfS = pd.read_excel(amp_conc_excel_path, sheetname = "samples", index_col=0)
        # reindex so the selected columns appear first. change all content to strings.
        selected_cols = ["samples", "Contains_Data"]
        dfS = tools.reindex_df_so_selected_cols_are_first(dfS, selected_cols).astype(str)
        # transfer sample names to dataframes with dose concentrations and response values
        dfAC_all['samples'] = dfS.samples
        df_resp_all['samples'] = dfS.samples

    # create output path for summary figure
    if sample_96_well == True:
        data_file_base = data_file[:-4]
    elif assay_type == "12_ampconc_12_well_platereader":
        data_file_base = data_file[:-5]
    # EC50_analysis_fig_basename = os.path.join(dff.loc[fn,"ofd"], data_file_base) + "EC50_analysis_fig_basename"
    # EC50_analysis_fig_basename_pdf = os.path.join(dff.loc[fn,"ofd_pdfs"], data_file_base) + "EC50_analysis_fig_basename"

    # create a view on the dataframes, so that it only shows the microplate data manually marked as "Contains_Data" = True
    dfdose = dfAC_all[dfAC_all.Contains_Data == True]
    dfresp = df_resp_all[df_resp_all.Contains_Data == True]
    # delete the two columns containing text, rather than data for plotting ("Contains_Data" and "samples")
    dfdose = dfdose.drop(["Contains_Data", "samples"], axis = 1)
    dfresp = dfresp.drop(["Contains_Data", "samples"], axis = 1)
    dict_dfe = {}

    #set up a single figure to contain 4 subplots
    n_plots_per_fig = 4
    nrows_in_each_fig = 2
    ncols_in_each_fig = 2
    dict_organising_subplots = create_dict_organising_subplots(n_plots_per_fig,n_rows=nrows_in_each_fig,n_cols=ncols_in_each_fig)
    #set the fontsize for the figure
    fig_fontsize = 6
    #set the default font for the figures
    plt.rc('font', family='sans-serif')
    plt.rc('font', serif='Helvetica Neue')
    plt.rc('text', usetex='false')
    plt.rcParams.update({'font.size': fig_fontsize})

    # iterate through all of the samples marked for analysis within the relevant AmpConc excel file
    for sNum in range(dfdose.shape[0]):
        #set up an empty nested dic
        # create a new dataframe for the evaluation of the EC50 calculations
        dfe = pd.DataFrame()
        # obtain the index letter for that sample
        sLet = dfdose.index[sNum]
        dfe.loc["sLet", sLet] = sLet
        dfe.loc["sNum", sLet] = sNum
        #add sample letter to final dictionary
        # obtain the name for that sample
        sample_name = str(dfAC_all.samples[sNum])
        dfe.loc["sample_name", sLet] = sample_name

        # set up the path for the image files to be saved in
        fig0_single_sample_png = os.path.join(dff.loc[fn,"output_folder"], "%s " % sLet + sample_name) + ".png"
        fig0_single_sample_pdf  = os.path.join(dff.loc[fn,"ofd_pdfs"], "%s " % sLet + sample_name) + ".pdf"

        #reindex so that only rows that contain data in both dataframes (i.e. dose and response data) are kept for analysis
        #take the index of both dataframes after NaN is removed.
        index_dfs = dfresp.loc[sLet,:].dropna().index
        #Find common elements using the set intersection function.
        cols_with_data_in_both_x_and_y = index_dfs.intersection(dfdose.loc[sLet,:].dropna().index)
        # reindex to drop the columns with text or boolean values
        x_orig = dfdose.loc[sLet,:].reindex(index = cols_with_data_in_both_x_and_y)
        dfe.loc["x_orig",sLet] = list(x_orig)

        y_orig = dfresp.loc[sLet,:].reindex(index = cols_with_data_in_both_x_and_y)
        # fixed upper limit (ful) data
        # when using the fixed upper limit in the response data (y-values), convert all higher values to the limit
        #y_ful = y_orig.apply(lambda x: dff.loc[fn, "yaxis fixed upper limit min"] if x > dff.loc[fn, "yaxis fixed upper limit min"] else x)

        # find where the first datapoint drops below the dff.loc[fn, "yaxis fixed upper limit min"]
        if y_orig.min() < dff.loc[fn, "yaxis fixed upper limit min"]:
            index_y_ful = np.min(np.where(y_orig < dff.loc[fn, "yaxis fixed upper limit min"]))
        else:
            # cells are overgrown, all datapoints are above ful, np.where returns a tuple, replace index with len(array)
            index_y_ful = len(y_orig)
        # select values above index_y_ful, which are above the fixed upper limit
        y_orig_data_above_ful = y_orig[:index_y_ful]
        y_orig_data_below_ful = y_orig[index_y_ful:]
        # normalise these datapoints between 0 and 1
        y_orig_data_above_ful_norm_0_1 = mf.normalise_0_1(y_orig_data_above_ful)[0]
        # calculate the width for the normalisation of the adjusted datapoints
        normalisation_width = dff.loc[fn, "yaxis fixed upper limit max"] - dff.loc[fn, "yaxis fixed upper limit min"]
        # convert the data which is normalised from 0-1, to normalised between the fixed_upper_limit_min and max
        y_orig_data_above_ful_norm = y_orig_data_above_ful_norm_0_1*normalisation_width + dff.loc[fn, "yaxis fixed upper limit min"]
        y_ful = np.append(y_orig_data_above_ful_norm, y_orig_data_below_ful)

        dfe.loc["y_orig",sLet] = list(y_orig)
        dfe.loc["y_ful",sLet] = list(y_ful)

        # in case there are problems with the curve fitting, plot the raw data alone
        # close any previous figures
        plt.close('all')
        # set a fontsize for figure annotation
        anno_fontsize = 10
        '''
        Fig01: Raw data, without any attempt at fitting.
        '''
        Plot_Nr = 1
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #Because this is a new sample, create a new figure (i.e. a new 2x2 canvas)
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)
        #add plot, set labels, set titles and grid
        axarr[row_nr, col_nr].scatter(x_orig, y_ful, color = t20[2], s=15)
        axarr[row_nr, col_nr].scatter(x_orig, y_orig, color = t20[0], s=5, label='_nolegend_')
        axarr[row_nr, col_nr].set_xlabel('ampicillin concentration (µg/ml)', fontsize = fig_fontsize)
        axarr[row_nr, col_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
        axarr[row_nr, col_nr].set_title("%s   %s" %(sLet, sample_name), fontsize = fig_fontsize)
        axarr[row_nr, col_nr].grid(True, color = '0.75')
        #set y-axis intercept
        ymin = -0.05
        # set the limit of the y-axis to 1 + 0.1, if all datapoints are very low
        ylim_max = y_orig.max() + 0.1 if y_orig.max() > 1.0 else 1.0
        axarr[row_nr, col_nr].set_ylim(ymin, ylim_max)
        # set the x-axis limit so that the legend does not hide too many data points
        # find the maximum Amp Conc in the whole experiment for that day
        maxAC = x_orig.max()
        # obtain the variable altering the extension of the x-axis
        x_axis_extension_after_ACmax_in_plot1 = dff.loc[fn, "x-axis extension in summary fig_0"]
        #define the limit of the x-axis as the maximum amp conc
        xlim_max_plot1 = maxAC + x_axis_extension_after_ACmax_in_plot1
        # set the x-axis limits
        xmin = -10
        axarr[row_nr, col_nr].set_xlim(xmin,xlim_max_plot1)
        #axarr[row_nr, col_nr].annotate(s=sLet, xy = (0.05,0.9), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction')
        axarr[row_nr, col_nr].annotate(s="fixed upper limit (ful) data", xy = (0.43,0.9), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', color = t20[2])
        axarr[row_nr, col_nr].annotate(s="original data", xy = (0.71,0.8), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', color = t20[0])

        #normalise the data to the highest value to improve results from the curve fitting algorithm
        #convert to float from object
        xnorm, x_orig_min,x_orig_max = mf.normalise_0_1(x_orig)
        ynorm, y_orig_min,y_orig_max = mf.normalise_0_1(y_orig)

        # I want the normalization to give similar values for the orig and the _ful datasets
        # therefore, if 0.2 is in the array, replace one of the values of 0.2 with the y_orig max value.
        # the normalization should then be quite similar
        #index_of_max_y_ful = y_ful.argmax()
        #ynorm_ful[index_of_max_y_ful] = y_orig_max
        # normalise the data with a fixed upper limit
        ynorm_ful, y_orig_min_ful,y_orig_max_ful = mf.normalise_0_1(y_ful)

        # #now replace that altered datapoint, which will be the fixed_upper_limit_min, de
        # ynorm_ful[index_of_max_y_ful] = mf.denormalise_0_1(dff.loc[fn, "yaxis fixed upper limit min"],y_orig_min_ful,y_orig_max_ful)

        #make an array of >250 datapoints representing the x-axis of the curve
        min = 0
        max = df_settings.loc["fitted_curve_xaxis_max","B"]
        n_datapoints = df_settings.loc["fitted_curve_n_datapoints","B"]
        x_fitted = np.linspace(min, max, n_datapoints)

        dfe.loc["x_fitted", sLet] = list(x_fitted)
        dfe.loc["amp_max", sLet] = x_orig_max
        dfe.loc["n_doseconc_tested", sLet] = len(x_orig)
        dfe.loc["x_orig", sLet] = list(x_orig)

        dfe.loc["y_orig_min", sLet] = y_orig_min
        dfe.loc["y_orig_max", sLet] = y_orig_max

        # count the number of datapoints in y_orig that are above and below the "live-dead cutoff" value
        dfe.loc["n_y_orig_dp_below_ldc", sLet] = len(np.where(y_orig < dff.loc[fn, "yaxis upper-lower cutoff"])[0])
        dfe.loc["n_y_orig_dp_above_ldc", sLet] = len(np.where(y_orig > dff.loc[fn, "yaxis upper-lower cutoff"])[0])

        if dfe.loc["n_y_orig_dp_below_ldc", sLet] < 2:
            """The cells are "overgrown" if any of the following are true
                - there are less than two datapoints above the live-dead cutoff value
            """
            EC50_ful = "overgrown"
            datasets = ["_orig", "_ful"]
            for d in datasets:
                dfe.loc["EC50{}".format(d), sLet], dfe.loc["rsquared{}".format(d), sLet] = "overgrown", 0
                dfe.loc["EC50_hill_eq{}".format(d), sLet], dfe.loc["n_highdose_datapoints{}".format(d), sLet] = 0, 0
                dfe.loc["n_lowdose_datapoints{}".format(d), sLet], dfe.loc["EC50_calculable{}".format(d), sLet] = 0, False
                dfe.loc["residuals_mean{}".format(d),sLet], dfe.loc["yvalues_for_curve_non_normalised{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["y_fitted_norm{}".format(d),sLet], dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = np.nan
            dfe.loc["x_fitted_non_normalised", sLet], dfe.loc["response_lowdose_datapoints_orig", sLet] = np.nan, np.nan

            # label the sample as not okay
            dfe.loc["y_orig_min", "%s_okay" % sLet] = False
            data_seems_okay_orig, data_seems_okay_ful = False, False
            # create empty values for the output of the fitting
            x_fitted_non_normalised, EC50_orig, EC50_hill_eq_orig, yvalues_for_curve_non_normalised_orig, \
            EC50_hill_eq_ful, yvalues_for_curve_non_normalised_ful = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            EC50_calculable_ful, EC50_calculable_orig = False, False

        elif any([y_orig[1] < df_settings.loc["min_responsevalue_second_doseconc","B"],
                 dfe.loc["n_y_orig_dp_above_ldc", sLet] < 2]):
            """For high-throughput LD50 calculations, the cells have "no growth" if any of the following are true:
                - the y-value of the second datapoint (second dose) is smaller than a fixed minimum value (min_responsevalue_second_doseconc)
                - there are less than two datapoints above a fixed value (yaxis upper-lower cutoff)
            """
            datasets = ["_orig", "_ful"]
            for d in datasets:
                dfe.loc["EC50{}".format(d), sLet], dfe.loc["rsquared{}".format(d), sLet] = "no_growth", 0
                dfe.loc["EC50_hill_eq{}".format(d), sLet], dfe.loc["n_highdose_datapoints{}".format(d), sLet] = 0, 0
                dfe.loc["n_lowdose_datapoints{}".format(d), sLet], dfe.loc["EC50_calculable{}".format(d), sLet] = 0, False
                dfe.loc["residuals_mean{}".format(d),sLet], dfe.loc["yvalues_for_curve_non_normalised{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["y_fitted_norm{}".format(d),sLet], dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = np.nan
            dfe.loc["x_fitted_non_normalised", sLet], dfe.loc["response_lowdose_datapoints_orig", sLet] = np.nan, np.nan

            # label the sample as not okay
            dfe.loc["y_orig_max", "%s_okay" % sLet] = False
            data_seems_okay_orig, data_seems_okay_ful = False, False
            # create empty values for the output of the fitting
            x_fitted_non_normalised, EC50_orig, EC50_hill_eq_orig, yvalues_for_curve_non_normalised_orig, \
            EC50_hill_eq_ful, yvalues_for_curve_non_normalised_ful = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            EC50_calculable_ful, EC50_calculable_orig = False, False
        else:
            #as a starting point, guess the sigmoidal constants
            hill_constants_guess = (1.0,0.0,0.5,10.0)
            #use the scipy optimise function to fit a curve to the data points
            #the out-of-date "warning" parameter probably relates to the -c:12: RuntimeWarning: overflow encountered in exp
            #see http://stackoverflow.com/questions/4359959/overflow-in-exp-in-scipy-numpy-in-python
            hill_constants_orig, cov_orig, infodict_orig, mesg_orig, ier_orig = scipy.optimize.leastsq(
                mf.residuals, hill_constants_guess, args=(mf.hill_eq, xnorm, ynorm), full_output=1) #warning=True

            # repeat for the "fixed upper limit" data
            hill_constants_ful, cov_ful, infodict_ful, mesg_ful, ier_ful = scipy.optimize.leastsq(
                mf.residuals, hill_constants_guess, args=(mf.hill_eq, xnorm,ynorm_ful), full_output=1) #warning=True

            # save the hill constants for later use
            dfe.loc["hill_constants_orig",sLet] = list(hill_constants_orig)
            dfe.loc["hill_constants_ful",sLet] = list(hill_constants_ful)

            # obtain the rsquared value for the fit of the curve to the data
            # code is from http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator
            ss_err = np.sum(np.array(infodict_orig['fvec'])**2)
            ss_err_ful = np.sum(np.array(infodict_ful['fvec'])**2)
            ss_tot = np.sum((np.array(ynorm)-ynorm.mean())**2)
            ss_tot_ful = np.sum((np.array(ynorm_ful)-ynorm.mean())**2)
            rsquared_orig = 1 - (ss_err/ss_tot)
            rsquared_ful = 1 - (ss_err_ful/ss_tot_ful)
            dfe.loc["rsquared_orig",sLet] = rsquared_orig
            dfe.loc["rsquared_ful",sLet] = rsquared_ful

            # also calculate the average residual as a rough "goodness of fit"
            # first apply the optimised function to the original datapoints
            y_fitted_xnorm = mf.hill_eq(hill_constants_orig, xnorm)
            y_fitted_xnorm_ful = mf.hill_eq(hill_constants_ful, xnorm)
            # calculate the residuals, as the (observed y-values)-(fitted y-values). Convert to positive floats.
            residuals_norm = abs(ynorm - y_fitted_xnorm)
            residuals_norm_ful = abs(ynorm_ful - y_fitted_xnorm_ful)
            # calculate mean residual
            residuals_norm_mean = residuals_norm.mean()
            residuals_norm_mean_ful = residuals_norm_ful.mean()
            # denormalise to the original y-value scale
            residuals_mean_orig = mf.denormalise_0_1(residuals_norm_mean, y_orig_min, y_orig_max)
            residuals_mean_ful = mf.denormalise_0_1(residuals_norm_mean_ful, y_orig_min_ful, y_orig_max_ful)
            # add to output series
            dfe.loc["residuals_mean_orig",sLet] = residuals_mean_orig
            dfe.loc["residuals_mean_ful",sLet] = residuals_mean_ful

            '''
            # Attempt to fit the hill equation to the data without normalization
            # DOES NOT FIT! Currently only the normalized values can be fit.
            x_orig_arr = np.array(x_orig).astype(float)
            y_orig_arr = np.array(y_orig).astype(float)
            hill_constants_guess_nonnorm = (2.0,0.0,100,5)
            hill_constants_nonnorm, cov_nonnorm, infodict_nonnorm, mesg_nonnorm, ier_nonnorm = scipy.optimize.leastsq(
                mf.residuals, hill_constants_guess_nonnorm, args=(mf.hill_eq, x_orig_arr,y_orig_arr), full_output=1) #warning=True
            yvalues_from_non_norm_Hill_fit = mf.hill_eq(hill_constants_orig_nonnorm, x_fitted)
            yvalues_from_Hill_guess = mf.hill_eq(hill_constants_guess, x_fitted)
            EC50_HillEq_nonnorm = hill_constants_orig_nonnorm[2]
            #printx("EC50_HillEq_nonnorm : %s" % EC50_HillEq_nonnorm)
            '''
            # obtain the constants of the optimised sigmoidal Hill function
            upper_orig, lower_orig, EC50_hill_eq_norm_orig, hillslope_orig = hill_constants_orig
            upper_ful, lower_ful, EC50_norm_ful, hillslope_ful = hill_constants_ful

            dfe.loc["EC50_hill_eq_orig",sLet] = mf.denormalise_0_1(hill_constants_orig[2], x_orig_min, x_orig_max)
            dfe.loc["EC50_hill_eq_ful",sLet] = mf.denormalise_0_1(hill_constants_ful[2], x_orig_min, x_orig_max)

            dfe.loc['upper_orig', sLet] = upper_orig
            dfe.loc['lower_orig', sLet] = lower_orig
            dfe.loc['EC50_hill_eq_norm_orig', sLet] = EC50_hill_eq_norm_orig
            dfe.loc['hillslope_orig', sLet] = hillslope_orig

            dfe.loc['upper_ful', sLet] = upper_ful
            dfe.loc['lower_ful', sLet] = lower_ful
            dfe.loc['EC50_norm_ful', sLet] = EC50_norm_ful
            dfe.loc['hillslope_ful', sLet] = hillslope_ful

            #calculate the value for y for the 1500 points
            y_fitted_norm_orig = mf.hill_eq(hill_constants_orig, x_fitted)
            dfe.loc['y_fitted_norm_orig', sLet] = list(y_fitted_norm_orig)
            y_fitted_norm_ful = mf.hill_eq(hill_constants_ful, x_fitted)
            dfe.loc['y_fitted_norm_ful', sLet] = list(y_fitted_norm_ful)

            #determine the EC50
            #the y-value of 50% cell density is calculated as the middle position in the curve
            #if the curve is perfectly symmetrical, the EC50 should equal the constant 'k' from the hill_constants
            curve_max_norm = y_fitted_norm_orig.max()
            curve_max_norm_ful = y_fitted_norm_ful.max()
            curve_min_norm_orig = y_fitted_norm_orig.min()
            curve_min_norm_ful = y_fitted_norm_ful.min()
            dfe.loc["curve_min_norm_orig",sLet] = curve_min_norm_orig
            dfe.loc["curve_min_norm_ful",sLet] = curve_min_norm_ful
            #y_value_curve_center_norm = (curve_max_norm - curve_min_norm_orig)/2 + curve_min_norm_orig
            # fix the minimum of the curve to be zero in all cases,  and calculate the center based on the max of the curve
            y_value_curve_center_norm = (curve_max_norm - 0)/2 + 0
            y_value_curve_center_norm_ful = (curve_max_norm_ful - 0)/2 + 0

            EC50_norm_bq_orig, EC50_calculable_orig, EC50_norm_bq_ful, \
            EC50_calculable_ful = calc_EC50_brent_eq(sLet, sample_name, hill_constants_orig, y_value_curve_center_norm,
                                                    hill_constants_ful, y_value_curve_center_norm_ful)

            # add if the EC50 was calculable to the summary Series
            if EC50_calculable_ful == True:
                dfe.loc["EC50_calculable_ful", sLet] = True
            else:
                dfe.loc["EC50_calculable_ful", sLet] = False
                dfe.loc["EC50_norm_bq_orig","%s_okay" % sLet] = False

            if EC50_calculable_orig == True:
                dfe.loc["EC50_calculable_orig", sLet] = True
            else:
                dfe.loc["EC50_calculable_orig", sLet] = False
                dfe.loc["EC50_norm_bq_orig","%s_okay" % sLet] = False

            #Plot the results
            '''
            Fig02: Normalised data, with fitted curve.
            '''
            Plot_Nr = Plot_Nr + 1
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
            #add plot, set labels, set titles and grid
            axarr[row_nr, col_nr].scatter(xnorm, ynorm_ful, color = t20[2], s=15)
            axarr[row_nr, col_nr].scatter(xnorm, ynorm, color = t20[0], s=5)
            axarr[row_nr, col_nr].plot(x_fitted, y_fitted_norm_ful, '-', color = t20[2], alpha = 0.8)
            axarr[row_nr, col_nr].plot(x_fitted,y_fitted_norm_orig,'-',color=t20[0], alpha = 0.8)
            # set xlabel, ylabel, title, grid, etc
            axarr[row_nr, col_nr].set_xlabel(df_settings.loc["x-axis (dose)","B"] + " normalised", fontsize = fig_fontsize)
            axarr[row_nr, col_nr].set_ylabel(df_settings.loc["y-axis (response)","B"],rotation='vertical', fontsize = fig_fontsize)
            # axarr[row_nr, col_nr].set_xlabel('ampicillin concentration normalised', fontsize = fig_fontsize)
            # axarr[row_nr, col_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
            axarr[row_nr, col_nr].set_title("%s   %s" %(sLet, sample_name), fontsize = fig_fontsize)
            axarr[row_nr, col_nr].grid(True, color = '0.75')
            #set y-axis intercept
            ymin_norm = -0.2
            xmin_norm = -0.1
            axarr[row_nr, col_nr].set_ylim(ymin_norm, 1.2)
            axarr[row_nr, col_nr].set_xlim(xmin_norm, 1.2)

            #horizontal line at global_half_absorption from 0 to EC50
            axarr[row_nr, col_nr].hlines(y=y_value_curve_center_norm_ful, xmin=xmin_norm, xmax=EC50_norm_bq_ful, colors = t20[2], linestyles='dashed', label='')
            axarr[row_nr, col_nr].hlines(y=y_value_curve_center_norm, xmin=xmin_norm, xmax=EC50_norm_bq_orig, colors = t20[0], linestyles='dashed', label='')
            #vertical line at EC50 from ymin_norm to global_half_absorption
            axarr[row_nr, col_nr].vlines(x=EC50_norm_bq_ful, ymin=ymin_norm, ymax=y_value_curve_center_norm_ful, colors = t20[2], linestyles='dashed', label='EC50 ful')
            axarr[row_nr, col_nr].vlines(x=EC50_norm_bq_orig, ymin=ymin_norm, ymax=y_value_curve_center_norm, colors = t20[0], linestyles='dashed', label='EC50_orig')
            # add the predicted EC50 from the hill equation to the graph (not implemented, brent eq works much better!)
            # if hill_constants_orig[2] < 1:
            #     axarr[row_nr, col_nr].hlines(y=y_value_curve_center_norm, xmin=0, xmax=hill_constants_orig[2], colors = t20[0], linestyles='dotted', label='', alpha = 0.9)
            #     axarr[row_nr, col_nr].hlines(y=y_value_curve_center_norm_ful, xmin=0, xmax=hill_constants_ful[2], colors = t20[2], linestyles='dotted', label='', alpha = 0.9)
            #     axarr[row_nr, col_nr].vlines(x=hill_constants_orig[2], ymin=ymin_norm, ymax=y_value_curve_center_norm, colors = t20[0], linestyles='dotted', alpha = 0.9, label='EC50 hill eq')
            #     axarr[row_nr, col_nr].vlines(x=hill_constants_ful[2], ymin=ymin_norm, ymax=y_value_curve_center_norm_ful, colors = t20[2], linestyles='dotted', alpha = 0.9, label='EC50 hill eq ful')
            lg = axarr[row_nr, col_nr].legend(loc=7)
            lg.draw_frame(False)

            # set annotation in top right
            axarr[row_nr, col_nr].annotate(s="normalised (ful) data", xy = (0.53,0.9), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', color=t20[2])
            axarr[row_nr, col_nr].annotate(s="normalised data", xy = (0.63,0.8), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', color=t20[0])
            #calculate the curve for the original x-values by multiplying by the maximum values
            #x_fitted_non_normalised = x_fitted * x_orig.max()
            x_fitted_non_normalised = mf.denormalise_0_1(x_fitted, x_orig_min, x_orig_max)
            dfe.loc["x_fitted_non_normalised",sLet] = x_fitted_non_normalised
            # EC50_orig = float(mf.denormalise_0_1(EC50_norm_bq_orig, x_orig_min, x_orig_max))
            # EC50_ful = float(mf.denormalise_0_1(EC50_norm_bq_ful, x_orig_min, x_orig_max))

            dfe.loc["EC50_norm_bq_orig",sLet] = EC50_norm_bq_orig
            dfe.loc["EC50_norm_bq_ful",sLet] = EC50_norm_bq_ful
            dfe.loc["EC50_orig",sLet] = float(mf.denormalise_0_1(EC50_norm_bq_orig, x_orig_min, x_orig_max))
            dfe.loc["EC50_ful",sLet] = float(mf.denormalise_0_1(EC50_norm_bq_ful, x_orig_min, x_orig_max))
            # dfe.loc["EC50_hill_eq_orig",sLet] = float('%02f' % EC50_hill_eq_orig)
            # dfe.loc["EC50_hill_eq_ful",sLet] = float('%02f' % EC50_hill_eq_ful)

            #calculate the curve for the original y-values
            yvalues_for_curve_non_normalised_orig = mf.denormalise_0_1(y_fitted_norm_orig, y_orig_min, y_orig_max)
            dfe.loc["yvalues_for_curve_non_normalised_orig",sLet] = yvalues_for_curve_non_normalised_orig
            yvalues_for_curve_non_normalised_ful = mf.denormalise_0_1(y_fitted_norm_ful, y_orig_min_ful, y_orig_max_ful)
            dfe.loc["yvalues_for_curve_non_normalised_ful",sLet] = yvalues_for_curve_non_normalised_ful

            #calculate the centre of the y-curve, used for the EC50
            y_value_curve_center_nonnorm = mf.denormalise_0_1(y_value_curve_center_norm, y_orig_min, y_orig_max)
            y_value_curve_center_nonnorm_ful = mf.denormalise_0_1(y_value_curve_center_norm_ful, y_orig_min_ful, y_orig_max_ful)

            # plot on the subplot with the original y values
            axarr[0, 0].plot(x_fitted_non_normalised, yvalues_for_curve_non_normalised_orig, '-', color = t20[0])
            axarr[0, 0].plot(x_fitted_non_normalised, yvalues_for_curve_non_normalised_ful, '-', color = t20[2], alpha = 0.8)

            #plot lines for the EC50 calculation
            #horizontal line at global_half_absorption from 0 to EC50
            axarr[0, 0].vlines(x=dfe.loc["EC50_ful",sLet], ymin=ymin, ymax=y_value_curve_center_nonnorm_ful, colors = t20[2], linestyles='dashed', label='EC50 ful')
            axarr[0, 0].vlines(x=dfe.loc["EC50_orig",sLet], ymin=ymin, ymax=y_value_curve_center_nonnorm, colors = t20[0], linestyles='dashed', label='EC50_orig')
            #vertical line at EC50 from ymin to global_half_absorption
            axarr[0, 0].hlines(y=y_value_curve_center_nonnorm_ful, xmin=xmin, xmax=dfe.loc["EC50_ful",sLet], colors = t20[2], linestyles='dashed', label='')
            axarr[0, 0].hlines(y=y_value_curve_center_nonnorm, xmin=xmin, xmax=dfe.loc["EC50_orig",sLet], colors = t20[0], linestyles='dashed', label='')
            # put in the location of the EC50 according to the hill equation (not currently used, Hill Eq constant is less reliable than the brent eq root-finding method)
            # if EC50_hill_eq_orig < x_orig_max:
            #     axarr[0, 0].hlines(y=y_value_curve_center_nonnorm, xmin=xmin, xmax=EC50_hill_eq_orig, colors = t20[0], linestyles='dotted', label='', alpha = 0.9)
            #     axarr[0, 0].vlines(x=EC50_hill_eq_orig, ymin=ymin, ymax=y_value_curve_center_nonnorm, colors = t20[0], linestyles='dotted', alpha = 0.9, label='EC50 hill eq')
            # if EC50_hill_eq_ful < x_orig_max:
            #     axarr[0, 0].hlines(y=y_value_curve_center_nonnorm_ful, xmin=xmin, xmax=EC50_hill_eq_ful, colors = t20[2], linestyles='dotted', label='', alpha = 0.9)
            #     axarr[0, 0].vlines(x=EC50_hill_eq_ful, ymin=ymin, ymax=y_value_curve_center_nonnorm_ful, colors = t20[2], linestyles='dotted', alpha = 0.9, label='EC50 hill eq ful')
            lg = axarr[0, 0].legend(loc=7)
            lg.draw_frame(False)

            # analyse the curve fit and data to judge whether the EC50 value is accurate
            dfe = judge_EC50_data(dfe, sLet, df_settings)

            dfe_index = pd.Series(dfe.index)
            dfe_index_ful = dfe_index[dfe_index.apply(lambda x : "_ful" in x)]
            dfe_index_orig = dfe_index[dfe_index.apply(lambda x : x[-5:] == "_orig")]
            dfe_ful = dfe.loc[dfe_index_ful]
            dfe_orig = dfe.loc[dfe_index_orig]
            # label data if any of the analysed variables are not okay
            if False in list(dfe_ful.loc[:,"%s_okay" % sLet]):
                data_seems_okay_ful = False
            else:
                data_seems_okay_ful = True
            if False in list(dfe_orig.loc[:,"%s_okay" % sLet]):
                data_seems_okay_orig = False
            else:
                data_seems_okay_orig = True

        # add the final judgement to the list, for use in the dataframe later
        dfe.loc["data_seems_okay_orig",sLet] = data_seems_okay_orig
        dfe.loc["data_seems_okay_ful",sLet] = data_seems_okay_ful

        '''
        Fig03: Notes on data quality 01
        '''
        Plot_Nr = Plot_Nr + 1
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        # write text on the empty figure
        yaxis_pos = np.linspace(0.9,0.1,10)
        xaxis_left = 0.05
        # set the xaxis position of the annotation for the ful samples
        xful = 0.6
        # set the xaxis position of the annotation for the orig samples
        xori = 0.8
        data_evaluation = "data seems good" if data_seems_okay_ful else "ful data needs checking"
        title_colour = "k" if data_seems_okay_ful else "r"
        title_summ = "Sample %s, %s" % (sLet, data_evaluation)
        axarr[row_nr, col_nr].set_title(title_summ, fontsize = anno_fontsize, color = title_colour, alpha=0.75)

        if sample_96_well == True:
            data_file_base = data_file[:-4]
            # combine columns to create the NBLA sample name
            NBLA = "".join(dfS.loc[sLet,"Uniprot#1":"Mutant#1"])
            # combine columns to create the CBLA sample name
            CBLA = "".join(dfS.loc[sLet,"Uniprot#2":"Mutant#2"])
            # extract version (sBla1.1, sBla1.2)
            Version = dfS.loc[sLet,"notes"].strip('_')
            axarr[row_nr, col_nr].annotate(s="NBLA: %s" % NBLA, xy = (xaxis_left,yaxis_pos[0]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="CBLA: %s" % CBLA, xy = (xaxis_left,yaxis_pos[1]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)

            # axarr[row_nr, col_nr].annotate(s=data_evaluation, xy = (0.6,yaxis_pos[0]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s=Version, ha='right', xy = (0.98,yaxis_pos[0]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)

        #add a table showing the rsquared and other aspects of the fit and dataset
        axarr[row_nr, col_nr].annotate(s="ful", xy = (xful,yaxis_pos[2]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
        axarr[row_nr, col_nr].annotate(s="orig", xy = (0.8,yaxis_pos[2]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
        axarr[row_nr, col_nr].annotate(s="EC50(µg/ml)", xy = (xaxis_left,yaxis_pos[3]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)

        if EC50_calculable_ful == True:
            if data_seems_okay_ful == True:
                axarr[row_nr, col_nr].annotate(s="%0.0f" % dfe.loc["EC50_ful",sLet], xy = (xful,yaxis_pos[3]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            else:
                if isinstance(dfe.loc["EC50_ful",sLet], str):
                    EC50_to_insert_ful = dfe.loc["EC50_ful",sLet]
                elif isinstance(dfe.loc["EC50_ful",sLet], float):
                    EC50_to_insert_ful = "%0.0f" % dfe.loc["EC50_ful",sLet]
                # insert error string or EC50, coloured red to indicate likely poor data
                axarr[row_nr, col_nr].annotate(s="%s" % EC50_to_insert_ful, xy = (xful,yaxis_pos[3]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = "r")
        if EC50_calculable_orig == True:
            if data_seems_okay_orig == True:
                axarr[row_nr, col_nr].annotate(s="%0.0f" % dfe.loc["EC50_orig",sLet], xy = (0.8,yaxis_pos[3]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            else:
                if isinstance(dfe.loc["EC50_orig",sLet], str):
                    EC50_to_insert_orig = dfe.loc["EC50_orig",sLet]
                elif isinstance(dfe.loc["EC50_orig",sLet], float):
                    EC50_to_insert_orig = "%0.0f" % dfe.loc["EC50_orig",sLet]
                # insert error string or EC50, coloured red to indicate likely poor data
                axarr[row_nr, col_nr].annotate(s="%s" % EC50_to_insert_orig, xy = (0.8,yaxis_pos[3]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = "r")
            # rsquared of the fit to the data
            axarr[row_nr, col_nr].annotate(s="rsquared", xy = (xaxis_left,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%0.2f"% rsquared_ful, xy = (xful,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75,color = dfe.loc["rsquared_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%0.2f"% rsquared_orig, xy = (xori,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75,color = dfe.loc["rsquared_orig","%s_colour" % sLet])

            # hillslope_orig of the fit to the data
            axarr[row_nr, col_nr].annotate(s="hillslope_orig", xy = (xaxis_left,yaxis_pos[5]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%0.1f"% hillslope_ful, xy = (xful,yaxis_pos[5]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["hillslope_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%0.1f"% hillslope_orig, xy = (xori,yaxis_pos[5]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["hillslope_orig","%s_colour" % sLet])
             # number of highdose datapoints
            axarr[row_nr, col_nr].annotate(s="n_highdose_datapoints", xy = (xaxis_left,yaxis_pos[6]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%i"% dfe.loc["n_highdose_datapoints_ful",sLet], xy = (xful,yaxis_pos[6]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["n_highdose_datapoints_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%i"% dfe.loc["n_highdose_datapoints_orig",sLet], xy = (xori,yaxis_pos[6]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["n_highdose_datapoints_orig","%s_colour" % sLet])
            # std of highdose datapoints
            axarr[row_nr, col_nr].annotate(s="std_highdose_datapoints", xy = (xaxis_left,yaxis_pos[7]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%0.2f"% dfe.loc["std_resp_highdose_datapoints_ful",sLet], xy = (xful,yaxis_pos[7]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["std_resp_highdose_datapoints_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%0.2f"% dfe.loc["std_resp_highdose_datapoints_orig",sLet], xy = (xori,yaxis_pos[7]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["std_resp_highdose_datapoints_orig","%s_colour" % sLet])

            # number of lowdose datapoints
            axarr[row_nr, col_nr].annotate(s="n_lowdose_datapoints", xy = (xaxis_left,yaxis_pos[8]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%i"% dfe.loc["n_lowdose_datapoints_ful",sLet], xy = (xful,yaxis_pos[8]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["n_lowdose_datapoints_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%i"% dfe.loc["n_lowdose_datapoints_orig",sLet], xy = (xori,yaxis_pos[8]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["n_lowdose_datapoints_orig","%s_colour" % sLet])
            # std of lowdose datapoints
            axarr[row_nr, col_nr].annotate(s="std_lowdose_datapoints", xy = (xaxis_left,yaxis_pos[9]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            axarr[row_nr, col_nr].annotate(s="%0.2f"% dfe.loc["std_resp_lowdose_datapoints_ful",sLet], xy = (xful,yaxis_pos[9]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["std_resp_lowdose_datapoints_ful","%s_colour" % sLet])
            axarr[row_nr, col_nr].annotate(s="%0.2f"% dfe.loc["std_resp_lowdose_datapoints_orig",sLet], xy = (xori,yaxis_pos[9]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color = dfe.loc["std_resp_lowdose_datapoints_orig","%s_colour" % sLet])

        '''
        Fig04: Notes on data quality 02
        '''
        Plot_Nr = Plot_Nr + 1
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]

        #add a table showing the rsquared and other aspects of the fit and dataset
        axarr[row_nr, col_nr].annotate(s="ful", xy = (xful,yaxis_pos[2]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
        axarr[row_nr, col_nr].annotate(s="orig", xy = (0.8,yaxis_pos[2]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)

        # add the stepsize near the EC50, which determines whether more dose concentrations are necessary
        axarr[row_nr, col_nr].annotate(s="amp. conc. stepsize", xy = (xaxis_left,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
        if EC50_calculable_ful == True:
            axarr[row_nr, col_nr].annotate(s=dfe.loc["doseconc_stepsize_at_EC50_ful",sLet],xy=(xful,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75, color=dfe.loc["doseconc_stepsize_at_EC50_ful","%s_colour" % sLet])
        if EC50_calculable_orig == True:
            axarr[row_nr, col_nr].annotate(s=dfe.loc["doseconc_stepsize_at_EC50_orig",sLet],xy=(xori,yaxis_pos[4]), fontsize=anno_fontsize, xytext=None, xycoords='axes fraction', alpha=0.75,color = dfe.loc["doseconc_stepsize_at_EC50_orig","%s_colour" % sLet])

        # Print data name and data evaluation
        data_evaluation_ful = "ful data seems good" if data_seems_okay_ful else "data needs checking"
        EC50_str_dict = {}
        for d in datasets:
            EC50_str_dict[d] =  "%.01f" % dfe.loc["EC50{}".format(d),sLet] if isinstance(dfe.loc["EC50{}".format(d),sLet], float) else dfe.loc["EC50{}".format(d),sLet]
        # EC50_for_summary = "%.01f" % dfe.loc["EC50_ful",sLet] if isinstance(dfe.loc["EC50_orig",sLet], float) else dfe.loc["EC50_orig",sLet]
        #print('%s %s : %s, EC50 (ful) = %s' % (sLet, sample_name, data_evaluation_ful, EC50_for_summary))
        print("{sLet} {s} : {EC}{p}({d1}) = {EC1},  {EC}{p}({d2}) = {EC2}".format(sLet=sLet, s=sample_name,
              EC=df_settings.loc["calculation_type","B"], p=str(df_settings.loc["percentage_response","B"]),
              d1=datasets[0][1:], EC1 = EC50_str_dict[datasets[0]], d2=datasets[1][1:],
              EC2=EC50_str_dict[datasets[1]]))
        # A sample_name : EC50(orig) = EC50_orig, EC50 (ful) = EC50

        #save figure with the fitted curve and calculated EC50 value
        fig.tight_layout()
        fig.savefig(fig0_single_sample_png, format='png', dpi=140)
        fig.savefig(fig0_single_sample_pdf, format='pdf')
        plt.close('all')
        dict_dfe[sLet] = dfe

    df_eval = pd.DataFrame()
    for sLet in dict_dfe:
        indiv_df = dict_dfe[sLet]
        df_eval = pd.concat([df_eval,indiv_df], axis=1)

    # transpose, (exchange column and rows)
    df_eval = df_eval.T
    # reindex so that selected columns are visible first
    df_eval_selected_cols = ['sample_name', 'EC50_ful', 'data_seems_okay_ful', 'data_seems_okay_orig','rsquared_orig', 'EC50_hill_eq_ful','n_highdose_datapoints_orig','n_lowdose_datapoints_orig','sNum']
    df_eval = tools.reindex_df_so_selected_cols_are_first(df_eval, df_eval_selected_cols)
    # create a column containing the sample letter and sample name combined
    df_eval['sLet_plus_sample_name'] = df_eval['sLet'] + " " + df_eval['sample_name']

    # divide the dataframe into two new dataframes, one for the values (_val) and another for the booleans
    # related to whether the data is okay or not. simply copy orig dataframe and drop unwanted columns.
    df_eval_values = df_eval.copy()
    for row in df_eval_values.index:
        if "_okay" in row:
            df_eval_values.drop(row, inplace=True)
        if "_colour" in row:
            df_eval_values.drop(row, inplace=True)
    # copy orig dataframe and drop unwanted columns.
    df_eval_bool = df_eval.copy()
    for row in df_eval_bool.index:
        if "_okay" not in row:
            df_eval_bool.drop(row, inplace=True)
    # drop empty columns
    df_eval_bool.dropna(how="all", axis=1, inplace=True)

    # sort by index (sLet)
    df_eval_values.sort_index(inplace=True)
    # # reindex so that selected columns are visible first
    # df_eval_values_selected_cols = ['sample_name', 'EC50_ful', 'data_seems_okay_ful', 'data_seems_okay_orig', 'rsquared_orig', 'EC50_hill_eq_ful','n_highdose_datapoints_orig','n_lowdose_datapoints_orig','sNum']
    # df_eval_values = tools.reindex_df_so_selected_cols_are_first(df_eval_values, df_eval_values_selected_cols)
    # sort by index (sLet)
    df_eval_bool.sort_index(inplace=True)
    # # reindex so that selected columns are visible first
    # df_eval_bool = tools.reindex_df_so_selected_cols_are_first(df_eval_bool, df_eval_values_selected_cols)

    #set up the summary figure to contain 2 subplots
    n_plots_per_fig = 2
    nrows_in_each_fig = 2
    ncols_in_each_fig = 1
    dict_organising_subplots = create_dict_organising_subplots(n_plots_per_fig,n_rows=nrows_in_each_fig,n_cols=ncols_in_each_fig)
    #set the fontsize for the figure
    fig_fontsize = 6

    for d in datasets:
        # go through all of the data, and set the label colour to red if the data_seems_okay is false
        df_eval_values["xlabel_colour{}".format(d)] = df_eval_values["data_seems_okay{}".format(d)].apply(lambda c: "k" if c == True else "r")

    '''
    Summ_Plot01, summ_Fig00: scattergram, original data with fitted curve
    '''
    Plot_Nr = 1
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    # create new figure (i.e., canvas)
    fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)
    if True in list(df_eval_values.loc[:,"EC50_calculable_orig"]):
        # plot the curves first, as this is used in the legend
        for sLet in df_eval_values.index:
            sNum = df_eval_values.loc[sLet,'sNum']
            axarr[row_nr].plot(df_eval_values.loc[sLet,"x_fitted_non_normalised"],
                                       df_eval_values.loc[sLet,"yvalues_for_curve_non_normalised_orig"],
                                       '-',
                                       color = t20[sNum],
                                       alpha = 0.8,
                                       label = sLet)
        # set the legend. Note that this is based on the last samples plotted
        lg = axarr[row_nr].legend(df_eval_values['sLet_plus_sample_name'], loc = 'upper right', fontsize = fig_fontsize)
        lg.draw_frame(False)
    # set xlabel, ylabel, title, grid, etc
    # axarr[row_nr].set_xlabel('ampicillin concentration (µg/ml)', fontsize = fig_fontsize)
    axarr[row_nr].set_xlabel(df_settings.loc["x-axis (dose)","B"], fontsize = fig_fontsize)
    # axarr[row_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
    axarr[row_nr].set_ylabel(df_settings.loc["y-axis (response)","B"],rotation='vertical', fontsize = fig_fontsize)
    # axarr[row_nr].set_title('Original Data', fontsize = fig_fontsize)
    axarr[row_nr].set_title('Original Data                              %s' % data_file,
                            fontsize = fig_fontsize, x = 0.22)
    axarr[row_nr].grid(True, color = '0.75')
    # plot the raw datapoints last, as these interfere with the legend
    for sLet in df_eval_values.index:
        sNum = df_eval_values.loc[sLet,'sNum']
        axarr[row_nr].scatter(df_eval_values.loc[sLet,"x_orig"],
                            df_eval_values.loc[sLet,"y_orig"],
                            color = t20[sNum],
                            alpha = 0.8,
                            s = 15,
                            label = sLet)
    if not True in list(df_eval_values.loc[:,"EC50_calculable_orig"]):
        # if none of the curves could be fitted, base the legend on the datapoints rather than the curves
        lg = axarr[row_nr].legend(df_eval_values['sLet_plus_sample_name'], loc = 'upper right', fontsize = fig_fontsize)
        lg.draw_frame(False)
    # set the x-axis limit so that the legend does not hide too many data points
    # find the maximum Amp Conc in the whole experiment for that day
    maxAC = dfdose.max().max()
    # obtain the variable altering the extension of the x-axis
    x_axis_extension_after_ACmax_in_summ_plot = dff.loc[fn, "x-axis extension in summary fig_0"]
    #define the limit of the x-axis as the maximum amp conc
    xlim_max = maxAC + x_axis_extension_after_ACmax_in_summ_plot
    # set the x-axis limits
    axarr[row_nr].set_xlim(-10,xlim_max)

    '''
    Summ_Plot02, summ_Fig01: barchart original data, EC50_orig (sLet on x-axis)
    '''
    Plot_Nr = Plot_Nr + 1
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    # use the sample letter alone as the name on the x-axis
    x_names = df_eval_values.sLet
    # the number of boxes in the bar-chart is the length of the initial dataset
    x_n_boxes = x_names.shape[0]
    # the position of the boxes in the bar-chart is the range(n_boxes)
    box_indices = range(x_n_boxes)
    # define the y-axis data
    y_EC50_orig = pd.to_numeric(df_eval_values.EC50_orig, errors="coerce")
    # use the 100 * mean of the residuals as the yerr
    yerr = df_eval_values.residuals_mean_orig.fillna(0)*100
    # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
    if True in list(df_eval_values.loc[:,"EC50_calculable_orig"]):
        bar_container = axarr[row_nr].bar(box_indices, y_EC50_orig, color = t20, yerr = yerr, align = "center",
                                          error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))
    # set the xticks
    axarr[row_nr].set_xticks(box_indices)
    # set the labels of the x-axis
    axarr[row_nr].set_xticklabels(x_names)
    for xtick, colour in zip(axarr[row_nr].get_xticklabels(), df_eval_values["xlabel_colour_orig"]):
        xtick.set_color(colour)
    # set the limits of the x-axis
    axarr[row_nr].set_xlim([-1, x_n_boxes])
    # set the x-axis title
    axarr[row_nr].set_xlabel("sample letter")
    # set the y-axis title
    axarr[row_nr].set_ylabel("EC50 (ug/ml)")
    #save figure
    fig.tight_layout()
    tools.savefig_if_necessary(savefig, fig, fig_nr, dff.loc[fn,"EC50_analysis_fig_basename"], formats = 'png', dpi=150)
    tools.savefig_if_necessary(savefig, fig, fig_nr, dff.loc[fn,"EC50_analysis_fig_basename_pdf"], formats = 'pdf')

    if True in list(df_eval_values.loc[:,"EC50_calculable_orig"]):
        '''
        Barchart01: barchart original data, EC50_orig (full name on x-axis)
        '''
        bar_fig_nr = 1
        # create new figure (i.e., canvas)
        fig, ax = plt.subplots(dpi=300)
        # use the sample letter plus the sample name as the name on the x-axis
        x_names = df_eval_values.sLet_plus_sample_name
        # use the mean residual as the yerr. replace np.nan with 0. multiply by 100 to give a visible error.
        yerr_ful = df_eval_values.residuals_mean_ful.fillna(0)*100
        ax.set_title('Original Data                              %s' % data_file,
                                fontsize = fig_fontsize, x = 0.22)
        # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
        bar_container = ax.bar(box_indices, y_EC50_orig, color = t20, yerr = yerr_ful, align = "center",
                                          error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))
        # set the xticks. apply the appropriate colour, as decided by the "judge" script
        ax.set_xticks(box_indices)
        # df_eval_values["xlabel_colour_orig"] = df_eval_values.data_seems_okay_orig.apply(lambda c: "k" if c == True else "r")
        for xtick, colour in zip(ax.get_xticklabels(), df_eval_values["xlabel_colour_orig"]):
            xtick.set_color(colour)
        # set the labels of the x-axis
        ax.set_xticklabels(x_names, rotation = 90)
        # set the limits of the x-axis
        ax.set_xlim([-1, x_n_boxes])
        # set the y-axis title
        ax.set_ylabel("EC50 (ug/ml)")

        #save figure
        fig.tight_layout()
        tools.savefig_if_necessary(savefig, fig, bar_fig_nr, dff.loc[fn,"EC50_analysis_fig_basename"] + "_bar", formats = 'png', dpi=150)
        tools.savefig_if_necessary(savefig, fig, bar_fig_nr, dff.loc[fn,"EC50_analysis_fig_basename_pdf"] + "_bar", formats = 'pdf')
    '''
    Summ_Plot04, summ_Fig03: scattergram, fixed upper limit (ful) data with fitted curve
    '''
    Plot_Nr = Plot_Nr + 1
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    # create new figure (i.e., canvas)
    fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)
    # if "yvalues_for_curve_non_normalised_ful" in df_eval_values.columns:
    # plot the curves first, as this is used in the legend)
    if True in list(df_eval_values.loc[:,"EC50_calculable_ful"]):
        for sLet in df_eval_values.index:
            sNum = df_eval_values.loc[sLet,'sNum']
            axarr[row_nr].plot(df_eval_values.loc[sLet,"x_fitted_non_normalised"],
                                       df_eval_values.loc[sLet,"yvalues_for_curve_non_normalised_ful"],
                                       '-',
                                       color = t20[sNum],
                                       alpha = 0.8,
                                       label = sLet)
        # set the legend. Note that this is based on the last samples plotted
        lg = axarr[row_nr].legend(df_eval_values['sLet_plus_sample_name'], loc = 'upper right', fontsize = fig_fontsize)
        lg.draw_frame(False)
    # set xlabel, ylabel, title, grid, etc
    # axarr[row_nr].set_xlabel('ampicillin concentration (µg/ml)', fontsize = fig_fontsize)
    axarr[row_nr].set_xlabel(df_settings.loc["x-axis (dose)","B"], fontsize = fig_fontsize)
    # axarr[row_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
    axarr[row_nr].set_ylabel(df_settings.loc["y-axis (response)","B"],rotation='vertical', fontsize = fig_fontsize)
    # axarr[row_nr].set_xlabel('ampicillin concentration (µg/ml)', fontsize = fig_fontsize)
    # axarr[row_nr].set_ylabel('cell density (A600)',rotation='vertical', fontsize = fig_fontsize)
    axarr[row_nr].set_title('Fixed Upper Limit (ful) Data                      %s' % data_file,
                            fontsize = fig_fontsize, x = 0.22)
    axarr[row_nr].grid(True, color = '0.75')
    # plot the raw datapoints last, as these interfere with the legend
    for sLet in df_eval_values.index:
        sNum = df_eval_values.loc[sLet,'sNum']
        axarr[row_nr].scatter(df_eval_values.loc[sLet,"x_orig"],
                            df_eval_values.loc[sLet,"y_ful"],
                            color = t20[sNum],
                            alpha = 0.6,
                            s=15)
    if not True in list(df_eval_values.loc[:,"EC50_calculable_orig"]):
        # if none of the curves could be fitted, base the legend on the datapoints rather than the curves
        lg = axarr[row_nr].legend(df_eval_values['sLet_plus_sample_name'], loc = 'upper right', fontsize = fig_fontsize)
        lg.draw_frame(False)
    # set the x-axis limit so that the legend does not hide too many data points
    # find the maximum Amp Conc in the whole experiment for that day
    maxAC = dfdose.max().max()
    # obtain the variable altering the extension of the x-axis
    x_axis_extension_after_ACmax_in_summ_plot = dff.loc[fn, "x-axis extension in summary fig_0"]
    #define the limit of the x-axis as the maximum amp conc
    xlim_max = maxAC + x_axis_extension_after_ACmax_in_summ_plot
    # set the x-axis limits
    axarr[row_nr].set_xlim(-10,xlim_max)
    axarr[row_nr].set_ylim(0,dff.loc[fn, "yaxis fixed upper limit max"])


    '''
    Summ_Plot05, summ_Fig03: barchart fixed upper limit (ful) data, EC50_ful (sLet on x-axis)
    '''
    Plot_Nr = Plot_Nr + 1
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    # use the sample letter alone as the name on the x-axis
    x_names = df_eval_values.sLet
    # the number of boxes in the bar-chart is the length of the initial dataset
    x_n_boxes = x_names.shape[0]
    # the position of the boxes in the bar-chart is the range(n_boxes)
    box_indices = range(x_n_boxes)
    # define the y-axis data. replace strings with np.nan
    y_EC50_ful = pd.to_numeric(df_eval_values.EC50_ful, errors = "coerce")
    # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
    error_kw={'ecolor':'k', 'linewidth':1, 'capsize':2, 'capthick':1}
    if True in list(df_eval_values.loc[:,"EC50_calculable_ful"]):
        bar_container = axarr[row_nr].bar(left=box_indices, height=y_EC50_ful, color = t20, yerr = yerr_ful, align = "center",
                                          error_kw=error_kw)
    # set the xticks
    axarr[row_nr].set_xticks(box_indices)
    # set the labels of the x-axis
    axarr[row_nr].set_xticklabels(x_names)
    for xtick, colour in zip(axarr[row_nr].get_xticklabels(), df_eval_values["xlabel_colour_ful"]):
        xtick.set_color(colour)
    # set the limits of the x-axis
    axarr[row_nr].set_xlim([-1, x_n_boxes])
    # set the x-axis title
    axarr[row_nr].set_xlabel("sample letter")
    # set the y-axis title
    axarr[row_nr].set_ylabel("EC50 (ug/ml)")
    #save figure
    fig.tight_layout()
    tools.savefig_if_necessary(savefig, fig, fig_nr, dff.loc[fn,"EC50_analysis_fig_basename"], formats = 'png', dpi=150)
    tools.savefig_if_necessary(savefig, fig, fig_nr, dff.loc[fn,"EC50_analysis_fig_basename_pdf"], formats = 'pdf')

    if True in list(df_eval_values.loc[:,"EC50_calculable_ful"]):
        '''
        Barchart02: barchartfixed upper limit (ful) data, EC50_ful (full name on x-axis)
        '''
        bar_fig_nr = 2
        # create new figure (i.e., canvas)
        fig, ax = plt.subplots(dpi=300)
        # use the sample letter plus the sample name as the name on the x-axis
        x_names = df_eval_values.sLet_plus_sample_name
        # create the new canvas (fig), figure (ax)
        ax.set_title("EC50 (fixed upper limit data)", fontsize = fig_fontsize)#
        ax.set_title('Fixed Upper Limit (ful) Data                      %s' % data_file,
                                fontsize = fig_fontsize, x = 0.22)
        # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
        bar_container = ax.bar(box_indices, y_EC50_ful, color = t20, yerr = yerr_ful, align = "center",
                                          error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))
        # set the xticks. apply the appropriate colour, as decided by the "judge" script
        ax.set_xticks(box_indices)
        for xtick, colour in zip(ax.get_xticklabels(), df_eval_values["xlabel_colour_ful"]):
            xtick.set_color(colour)
        # set the labels of the x-axis
        ax.set_xticklabels(x_names, rotation = 90)
        # set the limits of the x-axis
        ax.set_xlim([-1, x_n_boxes])
        # set the y-axis title
        ax.set_ylabel("EC50 (ug/ml)")

        #save figure
        fig.tight_layout()
        tools.savefig_if_necessary(savefig, fig, bar_fig_nr, dff.loc[fn,"EC50_analysis_fig_basename"] + "_barp", formats = 'png', dpi=150)
        tools.savefig_if_necessary(savefig, fig, bar_fig_nr, dff.loc[fn,"EC50_analysis_fig_basename_pdf"] + "_bar", formats = 'pdf')

    # convert arraylike and listlike data to strings
    list_arraylike_cols = ["x_orig", "y_orig","x_fitted_non_normalised","yvalues_for_curve_non_normalised_ful",
                           "yvalues_for_curve_non_normalised_orig","y_ful", "y_fitted_norm_orig",
                            "y_fitted_norm_ful", "x_fitted","indices_lowdose_datapoints_ful",
                            "indices_lowdose_datapoints_excl_nearest_EC50_ful","indices_lowdose_datapoints_orig",
                           "indices_lowdose_datapoints_excl_nearest_EC50_orig", "response_lowdose_datapoints_orig"]
    df_eval_values = tools.convert_listlike_cols_to_str(df_eval_values, list_arraylike_cols)

    # df_eval_bool = tools.convert_listlike_cols_to_str(df_eval_bool, list_arraylike_cols)
    # save evaluation dataframe to csv
    df_eval_values.to_csv(dff.loc[fn,"ofd_EC50_eval_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    df_eval_values.to_csv(dff.loc[fn,"ofd_EC50_eval_tabsep_csv"], sep="\t", quoting=csv.QUOTE_NONNUMERIC)
    # save evaluation dataframe to excel
    writer = pd.ExcelWriter(dff.loc[fn,"ofd_EC50_eval_excel"])#engine='xlsxwriter'
    df_eval_values.to_excel(writer, sheet_name="v_" + data_file[:20])
    df_eval_bool.to_excel(writer, sheet_name="b_" + data_file[:20])
    writer.save()
    writer.close()
    print('-------------------------------------')
    return df_eval_values

def calc_EC50_brent_eq(sLet, sample_name, hill_constants_orig, y_value_curve_center_norm, hill_constants_ful, y_value_curve_center_norm_ful):
    '''Use the brentq function to find the position on the x-axis at a particular y-value (for example #y_value_curve_center_norm),
    in some cases, this will result in ValueError: f(a) and f(b) must have different signs. This is due to a value
    outside of the initial range (between 0 and 1), and usually indicates a very poor fit
    http://stackoverflow.com/questions/22277982/how-to-find-50-point-after-curve-fitting-using-numpy/22288826#22288826
    '''
    try:
        EC50_norm_bq_orig = brentq(mf.hill_eq_brentq, 0.0, 1.0,
                                 args = (hill_constants_orig,
                                         y_value_curve_center_norm))
        EC50_calculable_orig = True
    except ValueError:
        try:
            print("ValueError encountered in brentq. Attempting wider scope for EC50. HIGHLY UNUSUAL! "
                  "Check data for %s_%s" % (sLet, sample_name))
            EC50_norm_bq_orig = brentq(mf.hill_eq_brentq, -1.0, 2.0,
                                     args = (hill_constants_orig,
                                             y_value_curve_center_norm))
            # dfe.loc["EC50_norm_bq_orig","%s_okay" % sLet] = False
            EC50_calculable_orig = True
        except ValueError:
            EC50_calculable_orig = False
            print("ValueError encountered in brentq_ful, even with wider scope! "
                  "Check data for %s_%s"
                  "EC50_calculable_orig = False" % (sLet, sample_name))

    try:
        EC50_norm_bq_ful = brentq(mf.hill_eq_brentq, 0.0, 1.0,
                                 args = (hill_constants_ful,
                                         y_value_curve_center_norm_ful))
        EC50_calculable_ful = True
    except ValueError:
        try:
            print("ValueError encountered in brentq_ful. Attempting wider scope for EC50. HIGHLY UNUSUAL! "
                  "Check data for %s_%s" % (sLet, sample_name))
            EC50_norm_bq_ful = brentq(mf.hill_eq_brentq, -1.0, 2.0,
                                     args = (hill_constants_ful,
                                             y_value_curve_center_norm_ful))
            # dfe.loc["EC50_norm_bq_ful","%s_okay" % sLet] = False
            EC50_calculable_ful = True
        except ValueError:
            EC50_calculable_ful = False
            print("ValueError encountered in brentq_ful, even with wider scope! "
                  "Check data for %s_%s"
                  "EC50_calculable = False" % (sLet, sample_name))
    return EC50_norm_bq_orig, EC50_calculable_orig, EC50_norm_bq_ful, EC50_calculable_ful

def standardise_doseconc_data(assay_type, df_dose_orig, df_resp_all, data_file_path):
    ''' Reformats the data from the different assay types so that they can all be analysed with the same followup
     scripts. Currently compatible with the 96-well data from the versamax microplate reader (8, 12 and 24 doseconc
     variations) and also the data from the Fluostar, which an read 12-well plates directly. For the 12-well format,
     the file contains both dose and response data (e.g. ampicillin concentration, and OD600), and will be loaded here.

    :param assay_type: string describing which assay type used
    :param df_dose_orig: unformatted ampicillin concentrations
    :param df_resp_all: formatted response data. ignored for 96-well samples, will be created for 12-well samples
    :param data_file_path: file path to load data for 12-well samples, ignored for 96-well samples.
    :return dfAC_all, df_resp_all: The two finished, formatted pandas dataframes that contain the dose concentration
    and response data respectively. The index and most columns should be identical.
    '''
    if assay_type == "8_ampconc":
        # transpose (flip index and column), and create a copy of dataframe
        dfAC_8 = df_dose_orig.T.copy()
        # isolate data from plate 1
        dfAC_plate1 = dfAC_8.iloc[:,:9].copy()
        # isolate data from plate 2
        dfAC_plate2 = dfAC_8.iloc[:,9:].copy()
        # join dataframes containing AmpConc for plates 1 and 2
        dfAC_all = pd.concat([dfAC_plate1,dfAC_plate2])
        # now convert to simple index
        # NOTE: after conversion, the letters on the microplate and now represented as AmpConc numbers (01-08)
        # The numbers across the top of the microplate are now represented as sample letters (A = sample 1, B= sample 2)
        # create numbers for new index (simple A01 style)
        AC_letters_8AC = list("ABCDEFGHIJKLMNOPQRSTUVWX")
        #AC_letters_8AC.append("orig_cols")
        dfAC_all.index = AC_letters_8AC
        # now convert to simple columns
        # create numbers for new index (simple A01 style)
        AC_numbers_8AC = ["%02d" % n for n in range(1,9)]
        AC_numbers_8AC.append('Contains_Data')
        #AC_numbers_8AC.append('orig_index')
        dfAC_all.columns = AC_numbers_8AC
        # because the columns in the orig excel file contain a mix of int and bool, the dtype is "object"
        # change the "Contains_Data" column to bool, rather than object
        dfAC_all['Contains_Data'] = dfAC_all.Contains_Data.apply(lambda xx: True if xx == "TRUE" else False)

    elif assay_type == "12_ampconc":
        reformatted_cols = ["%02d" % col for col in df_dose_orig.columns[:-1]]
        reformatted_cols.append(df_dose_orig.columns[-1])
        #create a copy of the orig dataframe
        dfAC_all = df_dose_orig.copy()
        dfAC_all.columns = reformatted_cols

    elif assay_type == "12_ampconc_12_well_platereader":
        #load the data into three dataframes, one with the amp concentrations, and the other with the A600 data and one samples names
        df_A600 = pd.read_excel(data_file_path, sheetname='A600', skiprows = 2)
        df_amp_conc = pd.read_excel(data_file_path, sheetname='amp_conc', skiprows = 2)
        number_of_initial_A600_columns = df_A600.shape[1]
        # df_samples = pd.read_excel(data_file_path, sheetname='samples', skiprows = 1, index_col = 0)

        #select only the useful data, assuming that the first 4 columns are ignored
        df_A600_data = df_A600.iloc[:,3:]
        # select useful data. ignore any columns longer than the A600 data
        df_amp_conc_data = df_amp_conc.iloc[:,3:number_of_initial_A600_columns]
        # drop any columns that are completely empty
        df_A600_data = df_A600_data.dropna(how="all", axis = 1)
        df_amp_conc_data = df_amp_conc_data.dropna(how="all", axis=1)
        # drop any rows that are completely empty
        df_A600_data = df_A600_data.dropna(how="all")
        df_amp_conc_data = df_amp_conc_data.dropna(how="all")
        # check if the column names are exactly the same
        if list(df_A600_data.columns) != list(df_amp_conc_data.columns):
            print("The sample names are different between the A600 and amp_conc tabs. "
                  "The 'samples' tab in the excel file is no longer used in the script, but is still "
                  "useful for generating accurate and consistent sample names. The names should be the same for both the "
                  "A600 and amp_conc tabs. Strongly recommend double-checking the names, and deleting any remaining "
                  "junk in empty columns. The OD600 names will be taken as TRUE, and amp_conc names ignored. "
                  "A600\n%sAmpconc\n%s " % (list(df_A600_data.columns), list(df_amp_conc_data.columns)))
        # df_samples_data = df_samples.iloc[:,0:].dropna(how="all", axis=1)
        temp_index = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789")
        # The excel files often contain junk, with undeleted data from previous experiments
        # Assume the df_A600_data has the "correct" sample names
        # count the number of ampicillin concentrations (according to the A600 tab)
        n_doseconc = df_A600_data.shape[0]
        # count the number of samples (according to the A600 tab)
        n_samples = df_A600_data.shape[1]
        # slice the sample data, transpose (exchange columns and rows) to standardise format weth 96-well data
        df_resp_all = df_A600_data.iloc[:n_doseconc + 1,:n_samples + 1].T
        # assume the A600 data has the correct names, and label the columns likewise
        df_resp_all['samples'] = df_A600_data.columns
        # label all columns (after the dropna step) as containing data
        df_resp_all['Contains_Data'] = True
        # label the index with letters (sLet)
        index_letters = temp_index[0:n_samples]
        df_resp_all.index = index_letters
        # slice the sample data, transpose (exchange columns and rows) to standardise format weth 96-well data
        dfAC_all = df_amp_conc_data.iloc[:n_doseconc + 1,:n_samples + 1].T
        # add the samples columnn to the amp conc dataframe, based on the "df_samples_data" dataframe.
        # Ignore any columns in the samples dataframe longer than the doseconc dataframe.
        dfAC_all['samples'] = df_A600_data.columns
        # label all columns (after the dropna step) as containing data
        dfAC_all['Contains_Data'] = True
        dfAC_all.index = index_letters

    elif assay_type == "24_ampconc":
        reformatted_cols = ["%02d" % col for col in df_dose_orig.columns[:-1]]
        reformatted_cols.append(df_dose_orig.columns[-1])
        #create a copy of the orig dataframe
        dfAC_24 = df_dose_orig.copy()
        dfAC_24.columns = reformatted_cols
        temp_index = list("ABCDEFGHIJKLMNOP")
        dfAC_24.index = temp_index
        rows_containing_01to12_AC = dfAC_24.index[0::2]
        rows_containing_13to24_AC = dfAC_24.index[1::2]
        final_index = list("ABCDEFGH")
        dfAC_01to12 = dfAC_24.loc[rows_containing_01to12_AC]
        dfAC_01to12.index = final_index
        dfAC_01to12 = dfAC_01to12.drop("Contains_Data", axis=1)
        dfAC_13to24 = dfAC_24.loc[rows_containing_13to24_AC]
        dfAC_13to24.index = final_index
        dfAC_13to24_cols = ["%02d" % r for r in range(13,25)]
        dfAC_13to24_cols.append(dfAC_24.columns[-1])
        dfAC_13to24.columns = dfAC_13to24_cols
        dfAC_all = pd.concat([dfAC_01to12, dfAC_13to24], axis=1)

    else:
        raise tools.DatafileError('assay type not identified')

    return dfAC_all, df_resp_all


def examine_input_datafile(fn, data_file, data_file_path, identifier_VERSAmax, list_96_well_doseconc_types):
    '''  Checks the file integrity of the input datafiles. If the datafile is .xlsx, it assumes that the input is
    from the 12-well platereader. For 96-well samples, where the input is a text file, it first looks for the keyword
    identifying as a VERSAmax file from the langosch lab (identifier_VERSAmax), and then looks for the keyword
    indicating the BlaTM 96-well format(either 8,12,24 amp concentrations, from list_96_well_doseconc_types).

    NOTE: This means the assay_type is DIRECTLY read from the input text datafile. The assay column on the settings file
    is for personal use only.

    :param fn: File number in list of files to analyse
    :param data_file: input file name
    :param data_file_path: input file path
    :param identifier_VERSAmax: string identifying data
    :param list_96_well_doseconc_types: list strings identifying assay_type
    :return assay_type (str), input_data_file_seems_okay (bool), sample_96_well (bool), sample_12_well_plate (bool):
    '''
    # checks the file integrity. First look for the keyword identifying as a VERSAmax file, and then looks for
    # text indicating the sBla 96-well format(either 8,12,24 amp concentrations)
    if os.path.exists(data_file_path):
        input_file_is_from_VERSAmax = False
        input_file_is_from_xlsx_template = False
        # if the data file ends in .xlsx, assume it is from the 12-well plate reader. Otherwise, open as textfile.
        if data_file_path[-5:] != ".xlsx":
            with open(data_file_path,'r') as f:
                for line in f:
                    # check if the identifying text is in the line
                    # NOTE: changes to the VERSAmax template introduction could delete this text,
                    # making the file unrecognisable to the python script
                    if identifier_VERSAmax in line:
                        input_file_is_from_VERSAmax = True
                    # check if any of the doseconc types are in the line.
                    if any (doseconctype in line for doseconctype in list_96_well_doseconc_types):
                        #define the assay type (8, 12 or 24 amp conc)
                        assay_type = line.strip('\n')
            if input_file_is_from_VERSAmax == False:
                raise tools.DatafileError("The input datafile does not contain the expected identifier"
                                               " phrase.\n'%s'\n Suggest checking the filepath and file integrity. "
                                               "For VersaMAX datafiles, double-check that when the data was"
                                               " exported, that both check-boxes were selected, so that"
                                                  " 'all sections on plate' were exported.\n"
                                                  "Relevant file:\n%s" % (identifier_VERSAmax,data_file_path))
        else:
            try:
                df_A600 = pd.read_excel(data_file_path, sheetname='A600', skiprows = 2)
                if "A" in df_A600.iloc[0,0] and "X" in df_A600.iloc[0,1]:
                    input_file_is_from_xlsx_template = True
                    assay_type = "12_ampconc_12_well_platereader"
                else:
                    raise LookupError("input file with data does not seem to match the xlsx template for the 12-well platereader :\n%s" %df_A600)
            except:
                raise LookupError("input file with data does not seem to match the xlsx template for the 12-well platereader.")
    else:
        raise IOError('File {0} does not exist'.format(data_file_path))

    # create a bool describing that the input files seem to be correct
    input_data_file_seems_okay, sample_96_well, sample_12_well_plate = False, False, False
    if isinstance(input_file_is_from_VERSAmax,bool) and input_file_is_from_VERSAmax == True and assay_type in list_96_well_doseconc_types:
        input_data_file_seems_okay = True
        sample_96_well = True
        print("Input file #%d, assay type = %s.\n%s. " % (fn, assay_type, data_file))
    elif isinstance(input_file_is_from_xlsx_template,bool) and input_file_is_from_xlsx_template == True:
        input_data_file_seems_okay = True
        print("Input file #%d, assay type = %s.\n%s. " % (fn, assay_type, data_file))
    else:
        print("Input file #%d. %s. Error. Data not readable as VERSAmax or Excel template file. "
              "Please check the file directory, filename and/or template." % (fn, data_file))

    return assay_type, input_data_file_seems_okay, sample_96_well, sample_12_well_plate

def read_96well_inputfiles(assay_type, input_data_file_seems_okay, data_file_path, amp_conc_excel_path):
    # identify where the data starts and ends in the csv file.
    if input_data_file_seems_okay:
        identifier_start_of_data = "Data_for_export"
        #identifier_end_of_data = "Instrument type: Instrument serial number:"
        with open(data_file_path,'r') as f:
            linenumber = 0
            # Create the initial local variable to keep the code inspection algorithm happy.
            linenumber_header = 0
            for line in f:
                # the data starts one line after the "data_for_export" text header
                if identifier_start_of_data in line:
                    linenumber_header = linenumber + 1
                # the data ends 3 lines before the "Instrument type:" text
                # not currently used in script, which assumes that there are always
                # 2 lines with text after the end of the data
                #if identifier_end_of_data in line:
                #    linenumber_end_data = linenumber - 3
                linenumber += 1
            # open the VERSAmax microplate reader exported csv file as a new dataframe, df_resp_orig (MICROPLATE READER DATA)
            df_resp_orig = pd.read_csv(data_file_path, sep = '\t', skiprows = linenumber_header, decimal = ",")
            #remove the last two rows, which should contain "~End" and "Instrument type: Instrument serial number:"
            df_resp_orig = df_resp_orig.dropna()
    # create new dataframe "df_dose_orig" to contain the dose concentrations (e.g. ampicillin conc, in LD50 assay)
    df_dose_orig = pd.read_excel(amp_conc_excel_path, sheetname = assay_type)
    df_resp_all = "temporary object, will be overwritten by standardisation of df_resp_orig, 12-well or 96well"
    return df_resp_orig, df_resp_all, df_dose_orig

