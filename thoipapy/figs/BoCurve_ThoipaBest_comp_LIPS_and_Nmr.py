import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from plotly import graph_objs as go, plotly
from plotly.tools import FigureFactory as FF
from scipy.interpolate import interp1d
from scipy.stats import ttest_ind
import plotly
import plotly.plotly as py
import tlabtools.tools as tools
import scipy.optimize
import scipy.stats
from tlabtools.mathfunctions import sine, residuals
from plotly.graph_objs import *
import thoipapy
plt.rcParams["font.family"] = "Verdana"
import matplotlib as mpl
from plotly import figure_factory as FF
import statistics
import glob
import os
import korbinian
plt.rcParams["font.family"] = "Verdana"
colour_dict = korbinian.utils.create_colour_lists()
colour_lists = tools.create_colour_lists()
color_thoipa = "k"
color_Lips = colour_dict["TUM_colours"]['TUM5']
blue1 = colour_dict["TUM_colours"]['TUM1']
blue5 = colour_dict["TUM_colours"]['TUM5']
black = "k"
TUMblue = colour_dict["TUM_colours"]['TUMBlue']

def thoipabest_vs_LIPS(bo_curve_handle,overlap_num,output_basename_LIPS,output_basename_LIPS_png,Fontsize,Linewidth):
    fig, ax = plt.subplots(figsize=(3.42, 3.42))
    ax1 = plt.subplot(311)
    # ax1 = plt.subplot(311)
    ax1.plot(overlap_num, bo_curve_handle['ThoipaTrainRatio68'].values, color=blue5, linestyle="-", linewidth=1.0,label="THOIPA best")
    ax1.plot(overlap_num, bo_curve_handle['ThoipaTrainLipsRatio68'].values, color=blue1, linestyle='--', linewidth=1.0,label="LIPS")
    ax1.fill_between(overlap_num, bo_curve_handle['ThoipaTrainRatio68'].values,bo_curve_handle['ThoipaTrainLipsRatio68'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax1.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth= Linewidth, label="random", alpha=0.6)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, ncol=1, fontsize=Fontsize, frameon=True)

    plt.setp(ax1.get_xticklabels(), fontsize=Fontsize, visible=False)

    ax1.set_ylabel("test crystal/NMR", fontsize=Fontsize)
    ax1.yaxis.set_label_position("right")

    plt.yticks(fontsize=Fontsize)

    # share x only
    ax2 = plt.subplot(312, sharex=ax1)
    # ax2 = plt.subplot(312, sharex=ax1)
    ax2.plot(overlap_num, bo_curve_handle['ThoipaTrain68Test13nmr'].values, color=blue5, linestyle="-", linewidth=1.0,label="THOIPA best")
    ax2.plot(overlap_num, bo_curve_handle['ThoipaTrain68LipsTest13nmr'].values, color=blue1, linestyle="--",linewidth=1.0, label="LIPS")
    ax2.fill_between(overlap_num, bo_curve_handle['ThoipaTrain68Test13nmr'].values,bo_curve_handle['ThoipaTrain68LipsTest13nmr'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax2.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth=Linewidth, label="random", alpha=0.6)
    # plt.plot(t, s2)
    # make these tick labels invisible
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax2.set_ylabel("test NMR", fontsize=Fontsize)
    ax2.yaxis.set_label_position("right")
    plt.yticks(fontsize=Fontsize)

    # share x and y
    ax3 = plt.subplot(313, sharex=ax1)
    # ax3 = plt.subplot(313, sharex=ax1)
    ax3.plot(overlap_num, bo_curve_handle['ThoipaEtraRatio68'].values, color=blue5, linestyle="-", linewidth=1.0,label="THOIPA best")
    ax3.plot(overlap_num, bo_curve_handle['ThoipaEtraLipsRatio68'].values, color=blue1, linestyle="--", linewidth=1.0,label="LIPS")
    ax3.fill_between(overlap_num, bo_curve_handle['ThoipaEtraRatio68'].values, bo_curve_handle['ThoipaEtraLipsRatio68'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax3.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth=Linewidth, label="random", alpha=0.6)


    ax3.set_ylabel("test ETRA", fontsize=Fontsize)
    ax3.yaxis.set_label_position("right")

    # ax3.xaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')
    ax1.xaxis.set_ticks_position('none')

    ax1.set_ylim(0.1, 6)
    ax2.set_ylim(0.1, 6)
    ax3.set_ylim(0.1, 6)

    x_list = bo_curve_handle.index.tolist()
    x_list = [i + 1 for i in x_list]

    ax3.set_xticks(x_list)
    fig.text(0.06, 0.5, "performance ratio", fontsize=Fontsize, ha='center', va='center', rotation='vertical')
    fig.text(0.5, 0.04, 'sample size', ha='center', va='center')
    fig.subplots_adjust(hspace=0.05)
    plt.xlim(0.9, 10)
    #plt.subplots_adjust(hspace=0.22, bottom=0.125)
    plt.xticks(fontsize=Fontsize)
    plt.yticks(fontsize=Fontsize)
    plt.rcParams['xtick.labelsize'] = Fontsize
    plt.rcParams['ytick.labelsize'] = Fontsize
    # plt.show()
    ax1.grid(False)
    ax2.grid(False)
    ax3.grid(False)
    plt.savefig(output_basename_LIPS)
    plt.savefig(output_basename_LIPS_png)
    plt.close()

################################################################################################################################################
###############################################################################################################################################
#####################################


def thoipabest_vs_thoipanmr(bo_curve_handle,overlap_num,output_basename_nmr,output_basename_nmr_png,Fontsize,Linewidth):
    fig, ax = plt.subplots(figsize=(3.42, 3.42))
    # plt.close()
    ax1 = plt.subplot(311)
    ax1.plot(overlap_num, bo_curve_handle['ThoipaTrainRatio68'].values, color=blue5, linestyle="-", linewidth=1.0,label="train crystal/NMR")
    ax1.plot(overlap_num, bo_curve_handle['ThoipaTrain68nmrLipsTestTrain68'].values, color=blue1, linestyle='--',linewidth=1.0, label="train NMR")
    ax1.fill_between(overlap_num, bo_curve_handle['ThoipaTrainRatio68'].values,
                     bo_curve_handle['ThoipaTrain68nmrLipsTestTrain68'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax1.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth=Linewidth, label="random", alpha=0.6)
    plt.setp(ax1.get_xticklabels(), fontsize=Fontsize, visible=False)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, ncol=1, fontsize=Fontsize, frameon=True)

    ax1.set_ylabel("test crystal/NMR", fontsize=Fontsize)
    ax1.yaxis.set_label_position("right")
    plt.yticks(fontsize=Fontsize)

    # share x only
    ax2 = plt.subplot(312, sharex=ax1)
    ax2.plot(overlap_num, bo_curve_handle['ThoipaTrain68Test13nmr'].values, color=blue5, linestyle="-", linewidth=1.0,label="train crystal/NMR")
    ax2.plot(overlap_num, bo_curve_handle['ThoipaTrainRatio68nmr'].values, color=blue1, linestyle="--", linewidth=1.0,label="train NMR")
    ax2.fill_between(overlap_num, bo_curve_handle['ThoipaTrain68Test13nmr'].values,
                     bo_curve_handle['ThoipaTrainRatio68nmr'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax2.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth=Linewidth, label="random", alpha=0.6)
    # plt.plot(t, s2)
    # make these tick labels invisible
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax2.set_ylabel("test NMR", fontsize=Fontsize)
    ax2.yaxis.set_label_position("right")

    plt.yticks(fontsize=Fontsize)

    # share x and y
    ax3 = plt.subplot(313, sharex=ax1)
    ax3.plot(overlap_num, bo_curve_handle['ThoipaEtraRatio68'].values, color=blue5, linestyle="-", linewidth=1.0,label="train crystal/NMR")
    ax3.plot(overlap_num, bo_curve_handle['ThoipaEtraRatio68nmr'].values, color=blue1, linestyle="--", linewidth=1.0,label="train NMR")
    ax3.fill_between(overlap_num, bo_curve_handle['ThoipaEtraRatio68'].values,bo_curve_handle['ThoipaEtraRatio68nmr'].values, color=color_thoipa, alpha=0.26,zorder=0)
    ax3.plot(overlap_num, [1] * 10, color=black, linestyle=':', linewidth=Linewidth, label="random", alpha=0.6)


    ax3.set_ylabel("test ETRA", fontsize=Fontsize)
    ax3.yaxis.set_label_position("right")

    # ax3.xaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')
    ax1.xaxis.set_ticks_position('none')

    ax1.set_ylim(0.1, 6)
    ax2.set_ylim(0.1, 6)
    ax3.set_ylim(0.1, 6)

    fig.text(0.06, 0.5, "performance ratio", fontsize=Fontsize, ha='center', va='center', rotation='vertical')
    fig.text(0.5, 0.04, 'sample size', ha='center', va='center')

    plt.xlim(0.9, 10)
    plt.ylim(0, 6)

    x_list = bo_curve_handle.index.tolist()
    x_list = [i+1 for i in x_list]
    ax3.set_xticks(x_list)

    #plt.subplots_adjust(hspace=0.22, bottom=0.125)
    plt.xticks(fontsize=Fontsize)
    plt.yticks(fontsize=Fontsize)
    plt.rcParams['xtick.labelsize'] = Fontsize
    plt.rcParams['ytick.labelsize'] = Fontsize
    fig.subplots_adjust(hspace=0.05)
    ax1.grid(False)
    ax2.grid(False)
    ax3.grid(False)

    plt.savefig(output_basename_nmr)
    plt.savefig(output_basename_nmr_png)

    plt.close()

def run_bocurve_comp( Fontsize, Width, Size, s,Linewidth):

    sys.stdout.write('\n~~~~~~~~~~~~                 starting run_bocurve_comp              ~~~~~~~~~~~~\n')
    sys.stdout.flush()

    Fontsize -= 7
    plt.rcParams["savefig.dpi"] = 240
    plt.rcParams.update({'font.size': 6})
    output_base_filename_nmr = "Thoipabest_vs_NMR.pdf"
    output_base_filename_nmr_png = "Thoipabest_vs_NMR.png"

    save_path_nmr = os.path.join("H:")
    settings_folder_nmr = os.path.dirname(save_path_nmr)
    output_folder_nmr = os.path.join(settings_folder_nmr,  "figs","FigBZ14_Bocurve_THOIPAbest_NMR")
    output_folder_pdf_nmr = os.path.join(output_folder_nmr, "pdf")
    output_basename_nmr = os.path.join(output_folder_nmr, output_base_filename_nmr)
    output_basename_nmr_png = os.path.join(output_folder_nmr, output_base_filename_nmr_png)

    output_base_filename_LIPS = "Thoipabest_vs_LIPS.pdf"
    output_base_filename_LIPS_png = "Thoipabest_vs_LIPS.png"
    save_path_LIPS = os.path.join("H:", "figs")
    settings_folder_LIPS = os.path.dirname(save_path_LIPS)
    output_folder_LIPS = os.path.join(settings_folder_LIPS, "figs","FigBZ15_Bocurve_THOIPAbest_LIPS")
    output_folder_pdf_LIPS = os.path.join(output_folder_LIPS, "pdf")

    output_basename_LIPS = os.path.join(output_folder_LIPS, output_base_filename_LIPS)
    output_basename_LIPS_png = os.path.join(output_folder_LIPS, output_base_filename_LIPS_png)

    list_paths = [output_folder_nmr, output_folder_LIPS, output_folder_pdf_nmr, output_folder_pdf_LIPS]
    for path in list_paths:
        if not os.path.exists(path):
            os.makedirs(path)

    bo_curve_file = r"H:\\figs\\FigBZ14_Bocurve_THOIPAbest_NMR\\zBo_Curve_save_for_python.csv"
    bo_curve_handle = pd.read_csv(bo_curve_file)
    bo_curve_handle.dropna(inplace=True)
    overlap_num = bo_curve_handle['OverlapNum'].values
    plt.close("all")
    thoipabest_vs_LIPS(bo_curve_handle, overlap_num, output_basename_LIPS,output_basename_LIPS_png,Fontsize,Linewidth)
    thoipabest_vs_thoipanmr(bo_curve_handle, overlap_num, output_basename_nmr,output_basename_nmr_png,Fontsize,Linewidth)

    sys.stdout.write('\n~~~~~~~~~~~~                 finished run_bocurve_comp              ~~~~~~~~~~~~')
    sys.stdout.flush()

