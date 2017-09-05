#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         Mark Teese
Created:        06.11.2015
Dependencies:   Python 3.x
                Numpy
                SciPy
                Pandas
                tlabassays tools (install from Bitbucket)
Purpose:        Fit a sine wave to a ToxR disruption index
Credits:        All sections by Mark Teese.
"""
from __future__ import unicode_literals
import matplotlib.pyplot as plt
import numpy as np
from thoipapy.Sine_Curve.tlabassays.mathfunctions import sine, sine_perfect_helix, residuals
import scipy.optimize
import os
import pandas as pd
import thoipapy.Sine_Curve.tlabtools.tools as tools

# test of VCS, 17.03.2016, *2

def fit_sin_to_data(x,y,title, output_folder, output_base_filename, n_aa=4):
    """ Takes experimental data (e.g. a ToxR disruption index) and attempts to fit the data to a sine curve.

    1) fits orig data to unbiased sine curve
    2) fits orig data to sine curve biased for alpha helical periodicity (3.6)
     - calculates disruption over a 3 aa window (e.g. AxxxAxxxA)
    3) fits 3aa window data to unbiased sine curve
    4) fits 3aa window data to sine curve biased for alpha helical periodicity (3.6)

    Parameters
    ----------
        x : array-like data of the x axis (e.g. aa position)
        y : array-like data of the y axis (e.g. disruption)
        title : title of the output figure
        output_basename : base filepath of output figures
        n_aa : number of amino acids separating 3 residues in the averaged window. Default is 4 (i, i+4).

    Saved Files and Figures
    -------
    basename.xlsx : excel file with the following tabs
        orig_disruption : original data, plus data averaged over a window
        sine_constants : sine constants from fitted curves
        yvalues_fitted : x and y datapoints used to draw the curves
    basename.png : line graph of original input data,
        fitted unbiased sine curve, and fitted sine curve with fixed periodicity
    basename + "_av_3aa_window.png" : line graph of original input data, averaged over a window of 3 aa (AxxxAxxxA)
        fitted unbiased sine curve, and fitted sine curve with fixed periodicity

    Returns
    -------
        figure plotting original data, and fitted curve
        excel summary
        sine_constants (tuple of a,b,c & d)
        periodicity_fitted_curve
        sine_constants : tuple of parameters a,b,c & d.
            From the unbiased fit to the original data (no window, no fixed 3.6 periodicity).
        periodicity_fitted_curve : the calculated periodicity of the unbiased curve fitted to the original data
    """
    #global orig_data_df, df_sine_constants
    output_folder_pdf = os.path.join(output_folder, "pdf")
    output_basename = os.path.join(output_folder, output_base_filename)
    output_basename_pdf = os.path.join(output_folder_pdf, output_base_filename)
    output_basename_av_3aa = os.path.join(output_folder, output_base_filename + "_av_3aa")
    output_basename_av_3aa_pdf = os.path.join(output_folder_pdf, output_base_filename + "_av_3aa")

    # create output folders, if necessary
    list_paths = [output_folder, output_folder_pdf]
    for path in list_paths:
        if not os.path.exists(path):
            os.makedirs(path)

    # set some figure settings
    plt.rcParams["savefig.dpi"] = 120
    plt.rc('font', family='Arial')
    fontsize = 12
    # create custom list of colours for figures
    colour_lists = tools.create_colour_lists()

    #####################################################################################################
    #                                                                                                   #
    #    CALCULATE SINE FIT TO AVERAGE DISRUPTION, WITHOUT AVERAGING OVER A WINDOW                      #
    #                                                                                                   #
    #####################################################################################################

    # create new fig
    fig, ax = plt.subplots()
    # plot the raw disruption data against the amino acid number
    ax.plot(x, y, color = colour_lists['TUM_colours']['TUMBlue'], label = "experimental")
    # create a number of x data points for smooth sine curves
    x_smooth = np.linspace(x.min(), x.max(), 500)

    # guess the constants in the initial curve
    sine_constants_guess = [1.0,1.7,0.2,0]
    # fit a sine curve to the data using the leastsq method
    sine_constants, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,
                                                                      sine_constants_guess,
                                                                      args=(sine,x,y),
                                                                      full_output=1)

    # print the periodicity of the curve that best fits the data (hopefully ~3.6, corresponding to an alpha-helix!)
    periodicity_fitted_curve = 2 * np.pi / sine_constants[1]
    #print("periodicity of fitted curve = %s" % periodicity_fitted_curve)
    # plot the fitted curve to the data
    yvalues_fitted = sine(sine_constants, x_smooth)
    ax.plot(x_smooth, yvalues_fitted, color = colour_lists['TUM_colours']['TUM4'], label = "fit to sine (unbiased)")
    ax.annotate(s='periodicity = %0.1f' % periodicity_fitted_curve,
                xy = (0.04,0.9),color = colour_lists['TUM_colours']['TUM4'],
                fontsize=fontsize, xytext=None, xycoords='axes fraction',
                alpha=0.75)

    #RECALCULATE SINE CONSTANTS ASSUMING A PERFECT HELIX (a=0.2, b=1.7, periodicity = 3.6)
    sine_constants_guess_perfhelix = [1.0,0.5]
    # fit a sine curve to the data using the leastsq method
    sine_constants_perfhelix, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = scipy.optimize.leastsq(residuals,
                                                                         sine_constants_guess_perfhelix,
                                                                         args=(sine_perfect_helix,x,y),
                                                                         full_output=1)
    # plot the fitted curve to the data
    yvalues_fitted_perfhelix = sine_perfect_helix(sine_constants_perfhelix, x_smooth)
    ax.plot(x_smooth, yvalues_fitted_perfhelix, color = colour_lists['TUM_accents']['green'], label = "fit to sine (fixed periodicity)")
    ax.annotate(s='periodicity = 3.6',
                xy = (0.04,0.8),color = colour_lists['TUM_accents']['green'],
                fontsize=fontsize, xytext=None, xycoords='axes fraction',
                alpha=0.75)

    set_ticks_and_labels_and_save(fig, ax, x, y, output_basename, output_basename_pdf, fontsize, title)

    #####################################################################################################
    #                                                                                                   #
    #    CALCULATE AVERAGE DISRUPTION OVER A WINDOW of residues on same side of perfect helix           #
    #                                                                                                   #
    #####################################################################################################

    # convert the x and y arrays to a new dataframe
    orig_data_df = pd.Series(x).to_frame(name='aa_position')
    orig_data_df['orig_disruption'] = y

    #set the aa position as the index
    orig_data_df.set_index("aa_position", inplace=True)
    # calculate the average disrupt index for each across a window
    orig_data_df = tools.calc_av_over_3aa_window(df = orig_data_df, data_col = "orig_disruption", new_col = "av_disrupt_3aa", n_aa = n_aa)
    x = orig_data_df.index
    y = orig_data_df.av_disrupt_3aa
    # create new fig with the averaged 3aa window data

    plt.close("all")
    fig, ax = plt.subplots()
    # plot the raw disruption data against the amino acid number
    ax.plot(x, y, color = colour_lists['TUM_colours']['TUMBlue'], label = "experimental")
    # create a number of x data points for smooth sine curves
    x_smooth = np.linspace(x.min(), x.max(), 500)
    # cguess the constants in the initial curve
    sine_constants_guess = [1.0,1.7,0.2,0]
    # fit a sine curve to the data using the leastsq method
    sine_constants, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,
                                                                         sine_constants_guess,
                                                                         args=(sine,x,y),
                                                                         full_output=1)
    # calculate the periodicity of the curve that best fits the data (hopefully ~3.6, corresponding to an alpha-helix!)
    periodicity_fitted_curve = 2 * np.pi / sine_constants[1]
    #print("periodicity of fitted curve = %s" % periodicity_fitted_curve)
    # plot the fitted curve to the data
    yvalues_fitted = sine(sine_constants, x_smooth)
    ax.plot(x_smooth, yvalues_fitted, color = colour_lists['TUM_colours']['TUM4'], label = "fit to sine (unbiased)")
    ax.annotate(s='periodicity = %0.1f' % periodicity_fitted_curve,
                xy = (0.04,0.9),color = colour_lists['TUM_colours']['TUM4'],
                fontsize=fontsize, xytext=None, xycoords='axes fraction',
                alpha=0.75)

    #RECALCULATE SINE CONSTANTS ASSUMING A PERFECT HELIX (a=0.2, b=1.7, periodicity = 3.6)
    sine_constants_guess_perfhelix = [1.0,0.5]
    # fit a sine curve to the data using the leastsq method
    sine_constants_perfhelix, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = scipy.optimize.leastsq(residuals,
                                                                         sine_constants_guess_perfhelix,
                                                                         args=(sine_perfect_helix,x,y),
                                                                         full_output=1)
    # plot the fitted curve to the data
    yvalues_fitted_perfhelix = sine_perfect_helix(sine_constants_perfhelix, x_smooth)
    ax.plot(x_smooth, yvalues_fitted_perfhelix, color = colour_lists['TUM_accents']['green'], label = "fit to sine (fixed periodicity)")
    ax.annotate(s='periodicity = 3.6',
                xy = (0.04,0.8),color = colour_lists['TUM_accents']['green'],
                fontsize=fontsize, xytext=None, xycoords='axes fraction',
                alpha=0.75)

    title_av_3aa = title + " av 3aa"
    set_ticks_and_labels_and_save(fig, ax, x, y, output_basename_av_3aa, output_basename_av_3aa_pdf, fontsize, title_av_3aa)

    # save data to new excel file
    writer = pd.ExcelWriter(output_basename + ".xlsx")
    orig_data_df.to_excel(writer, sheet_name='orig_disruption')
    # save the constants for the sine curve fitted to the data
    sine_constants_series = pd.Series(sine_constants, index = ["a","b","c","d"])
    sine_constants_series["periodicity"] = periodicity_fitted_curve
    df_sine_constants = sine_constants_series.to_frame(name="sine_constants_unbiased")
    df_sine_constants.loc["c":"d","sine_constants_fixed_periodicity"] = sine_constants_perfhelix

    df_sine_constants.to_excel(writer, sheet_name = "sine_constants")
    # save the datapoints for the fitted curve, so that it can be recreated in excel, etc.
    fitted_curve_as_series = pd.Series(yvalues_fitted, index = x_smooth)
    df_fitted_curve = fitted_curve_as_series.to_frame(name='yvalues_fitted_unbiased')
    df_fitted_curve['yvalues_fitted_fixed_periodicity'] = yvalues_fitted_perfhelix
    df_fitted_curve.to_excel(writer, sheet_name='yvalues_fitted')
    writer.save()
    writer.close()
    return sine_constants, periodicity_fitted_curve

def set_ticks_and_labels_and_save(fig, ax, x, y, output_basename, output_basename_pdf, fontsize, title):
    """ Set the ticks and other parameters in the output figure, and save as png and pdf

    Parameters
    ----------
        fig : figure object
        ax : ax (plot) object
        x : array-like data of the x axis (e.g. aa position)
        y : array-like data of the y axis (e.g. disruption)
        output_basename : base filepath of output figures
        output_basename_pdf : base filepath of output pdf figures
        fontsize : figure fontsize
        title : figure title

    Saved Files and Figures
    -------
    output_basename + ".png" : output figure, png
    output_basename_pdf + ".pdf" : output figure, pdf
    """

    """ Adjusts some of the figure parameters, and saves as a png and pdf using the basepath indicated
    """
    ax.set_xticks(x)
    ax.set_xticklabels(x, fontsize  = fontsize, rotation = 90)
    ax.set_xlabel('residue', fontsize=fontsize)


    #### Yao Changed here from average disruption to ToxR activity####
    ax.set_ylabel('Predicted disruption', fontsize=fontsize)




    #ax.set_ylim(0,1.4)
    y_dist = max(y)-min(y)
    ylim_min = min(y) - y_dist*0.2 if min(y) - y_dist*0.2 > 0 else 0
    ylim_max = max(y) + y_dist*0.2
    #ax.set_ylim(min(y),max(y) + 0.4)
    ax.set_ylim(ylim_min, ylim_max + 0.4)
    ax.set_title(title, fontsize=fontsize)
    # show the legend(as previously labelled)
    ax.legend(fontsize=fontsize)
    # save figure
    fig.savefig(output_basename + ".png")
    fig.savefig(output_basename_pdf + ".pdf", format="pdf")

def extract_aa_pos_from_sample_names(df_with_sn_as_index, start_aa, end_aa, newcol_name):
    """Extract the amino acid position from sample names in the index of a dataframe.

    Note this cannot be used if the data contains double-mutants. Consider using a dictionary approach instead.

    Parameters
    ----------
    df_with_sn_as_index : dataframe with the sample name as the index
        Dataframe to be altered, will be returned with the new column
    start_aa : int, starting amino acid in range of amino acids mutated
        UniProt style, NOT python indexing style
    end_aa : int, end amino acid in range of amino acids mutated
        UniProt style where the number is the last aa to be included, NOT python indexing style
    newcol_name : name of the new column containing the index

    Returns
    -------
    df_with_sn_as_index : original dataframe, with the new column
    """

    for i in range(start_aa,end_aa+1):
        # for each sample name in the dataframe with unique samples
        for sn in df_with_sn_as_index.index:
            # if the amino acid number (233, 234, etc) is in the sample name
            if str(i) in sn:
                # add the amino acid number to a new column
                df_with_sn_as_index.loc[sn, newcol_name] = i