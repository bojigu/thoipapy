import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
plt.rcParams["font.family"] = "Verdana"
import statistics
import os
from thoipapy.utils import create_colour_lists
colour_dict = create_colour_lists()
color_thoipa = "k"
color_Lips = colour_dict["TUM_colours"]['TUM5']
blue1 = colour_dict["TUM_colours"]['TUM1']
blue5 = colour_dict["TUM_colours"]['TUM5']
black = "k"
TUMblue = colour_dict["TUM_colours"]['TUMBlue']

def FigZB_07_hardlinked(Fontsize, Width, Size, s):
    """ use bar chart to show the interface based on the cutoff and limit to max 7 min 3,
    use different to show pattern

    Parameters
    ----------
    s_path : str
        Path to settings excel file
    acc_list : list
        List of uniprot accession numbers.
    df_list: list
        List of

    Dataframesinterface_con_cov_lip
    ----------
    df_MU
        example of notes, showing what is in the code.

    Saved Files and Figures
    -----------------------
    figXX : png
        example of notes describing output files

    """
    sys.stdout.write('\n~~~~~~~~~~~~                 starting FigZB_07_hardlinked              ~~~~~~~~~~~~\n')
    sys.stdout.flush()
    Width=0.35

    output_base_filename = ""
    output_base_filename_excel = "Data.xlsx"

    save_path=r"H:\figs\FigBZ07_AverageDI_FractionHighDI"
    settings_folder = os.path.dirname(save_path)
    output_folder = os.path.join(settings_folder, "FigBZ07_AverageDI_FractionHighDI","Figures")
    output_folder_pdf = os.path.join(output_folder, "pdf")


    output_basename = os.path.join(output_folder, output_base_filename)
    output_basename_excel= os.path.join(output_folder, output_base_filename_excel)

    output_basename_pdf = os.path.join(output_folder_pdf, output_base_filename)

    list_paths = [output_folder,output_folder_pdf]
    for path in list_paths:
        if not os.path.exists(path):
            os.makedirs(path)

    # read name from name file
    #name_file = r"H:\data_xy\protein_names.xlsx"
    names_excel_path = os.path.join(s["dropbox_dir"], "ETRA_NMR_names.xlsx")
    df_name_file = pd.read_excel(names_excel_path, sheet_name='proteins', index_col=0)


    #read crustal raw data
    crystal_excel = r"H:\figs\FigBZ07_AverageDI_FractionHighDI\Crystal_DI_rawdata.csv"
    df_crystal = pd.read_csv(crystal_excel, index_col=0)
    df_crystal.reset_index(inplace=True)
    # read NMR raw data
    NMR_excel = r"H:\figs\FigBZ07_AverageDI_FractionHighDI\NMR_DI_rawdata.csv"
    df_NMR = pd.read_csv(NMR_excel, index_col=0)

    # only contain NMRdata
    df_name_file_NMR = df_name_file[df_name_file.database == "NMR"]

    for i_NMR in df_NMR.index:
        Detail_name = df_name_file_NMR.loc[i_NMR, "short_name"]
        df_NMR.loc[i_NMR,"Protein"] = Detail_name

    df_NMR.reset_index(inplace=True)

    # read ETRA raw data
    ETRA_excel = r"H:\figs\FigBZ07_AverageDI_FractionHighDI\ETRA_DI_rawdata.csv"
    df_ETRA = pd.read_csv(ETRA_excel, index_col=0)

    # only contain ETRA data
    df_name_file_ETRA = df_name_file[df_name_file.database == "ETRA"]

    for i_ETRA in df_ETRA.index:
        Detail_name = df_name_file_ETRA.loc[i_ETRA, "short_name"]
        df_ETRA.loc[i_ETRA, "Protein"] = Detail_name
    df_ETRA.reset_index(inplace=True)


    list_plot=[df_crystal,df_NMR,df_ETRA]
    list_name=["Crystal","NMR","ETRA"]
    list_size_x=[12,8,8]
    list_size_y = [6,6,6]
    list_fontsize=[Fontsize,Fontsize,Fontsize]
    list_space=[0.025,0.05,0.05]
    list_bottom=[0.17,0.35,0.35]

    #plot for Average DI
    for n, m ,sx,sy,fs,space,bottom in zip(list_plot,list_name,list_size_x,list_size_y,list_fontsize,list_space,list_bottom):
        x = n.index.tolist()
        x_1 = [int(i) + Width for i in x]
        x_label_location= [int(i) + Width/2 for i in x]
        df_plot_interface=n["AverageInter"].tolist()
        df_plot_noninterface = n["AverageNoninter"].tolist()
        protein_name=n["Protein"].tolist()

        fig, ax = plt.subplots(figsize=(sx,sy))
        ax.bar(x, df_plot_interface, width=Width, color='k', align='center', label="interface",alpha=1)
        ax.bar(x_1, df_plot_noninterface, width=Width, color='grey', align='center', label="non-interface",alpha=0.6)
        ax.set_ylabel("average coevolution (DI)", fontsize=fs)
        ax.set_xticks(x_label_location)
        plt.title(m, fontsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        ax.set_xticklabels(protein_name, fontsize=fs, rotation=90)

        plt.grid(False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, ncol=2, fontsize=fs,frameon=True)
        plt.gcf().subplots_adjust(bottom=bottom)

        plt.margins(space, 0.1)
        ax.set_ylim(0, )
        fig.tight_layout()
        plt.savefig(output_basename + m+"_average.png")
        plt.savefig(output_basename_pdf + m + "_average.pdf")
        plt.close()

    # plot for fraction DI
    for n, m, sx, sy, fs,space, bottom in zip(list_plot, list_name, list_size_x, list_size_y, list_fontsize,list_space,list_bottom):
        x = n.index.tolist()

        x_1 = [int(i) + Width for i in x]
        x_label_location = [int(i) + Width / 2 for i in x]

        df_plot_interface = n["FractionInter"].tolist()

        df_plot_noninterface = n["FractionNoninter"].tolist()
        protein_name = n["Protein"].tolist()

        fig, ax = plt.subplots(figsize=(sx, sy))

        ax.bar(x, df_plot_interface, width=Width, color='k', align='center', label="interface", alpha=1)
        ax.bar(x_1, df_plot_noninterface, width=Width, color='grey', align='center', label="non-interface",alpha=0.6)
        ax.set_ylabel("fraction coevolution (DI)", fontsize=fs)
        ax.set_xticks(x_label_location)
        plt.title(m, fontsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        ax.set_xticklabels(protein_name, fontsize=fs, rotation=90)
        ax.set_ylim(0, 1.3)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, ncol=2, fontsize=fs,frameon=True)
        plt.gcf().subplots_adjust(bottom=bottom)
        plt.margins(space, 0.1)
        plt.grid(False)
        fig.tight_layout()
        plt.savefig(output_basename + m + "_fraction.png")
        plt.savefig(output_basename_pdf + m + "_fraction.pdf")
        plt.close()


    #plot average DI for 3 dataset use barchart.
    list_mean=[]
    for n,m  in zip(list_plot, list_name):
        mean_AverageInter=np.mean(n["AverageInter"].tolist())
        SD_mean_AverageInter = statistics.stdev(n["AverageInter"].tolist())

        mean_AverageNoninter = np.mean(n["AverageNoninter"].tolist())
        SD_mean_AverageNoninter = statistics.stdev(n["AverageNoninter"].tolist())

        mean_FractionInter = np.mean(n["FractionInter"].tolist())
        SD_mean_FractionInter = statistics.stdev(n["FractionInter"].tolist())

        mean_FractionNoninter = np.mean(n["FractionNoninter"].tolist())
        SD_mean_FractionNoninter = statistics.stdev(n["FractionNoninter"].tolist())


        list_mean.append((m,mean_AverageInter,SD_mean_AverageInter,mean_AverageNoninter,SD_mean_AverageNoninter,
                   mean_FractionInter,SD_mean_FractionInter,mean_FractionNoninter,SD_mean_FractionNoninter))

    df_list_mean=pd.DataFrame(list_mean)
    df_list_mean.columns=["Dataset","mean_AverageInter","SD_mean_AverageInter","mean_AverageNoninter","SD_mean_AverageNoninter",
                         "mean_FractionInter","SD_mean_FractionInter","mean_FractionNoninter","SD_mean_FractionNoninter"]

    df_list_mean.reset_index(inplace=True)

    # save excel
    writer = pd.ExcelWriter(output_basename_excel)
    df_list_mean.to_excel(writer, 'Sheet1')

    #plot for mean average DI
    fig, ax = plt.subplots(figsize=(4,6))
    x = df_list_mean.index.tolist()
    x_1 = [int(i) + Width for i in x]

    plot_inter=df_list_mean["mean_AverageInter"].tolist()
    SD_inter=df_list_mean["SD_mean_AverageInter"].tolist()

    plot_noninter = df_list_mean["mean_AverageNoninter"].tolist()
    SD_noninter = df_list_mean["SD_mean_AverageNoninter"].tolist()

    name= df_list_mean["Dataset"].tolist()

    ax.bar(x, plot_inter, width=Width, color='k', align='center',yerr=SD_inter, label="interface", alpha=0.6)
    ax.bar(x_1, plot_noninter, width=Width, color='grey', align='center',yerr=SD_noninter, label="non-interface",hatch="///", alpha=0.6)
    ax.set_ylabel("mean average coevolution(DI)", fontsize=Fontsize)
    ax.set_xticks(x)
    ax.tick_params(axis='y', labelsize=Fontsize)
    ax.set_xticklabels(name, fontsize=Fontsize, rotation=0)
    mpl.rcParams['hatch.linewidth'] = 2.5

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, ncol=1, fontsize=Fontsize,frameon=True)
    plt.gcf().subplots_adjust(bottom=0.18)
    plt.axis('tight')
    fig.tight_layout()
    plt.grid(False)
    plt.savefig(output_basename + m + "_mean_average_DI.png")
    plt.savefig(output_basename_pdf + m + "_mean_average_DI.pdf")
    plt.close()

    #plot mean fraction DI
    fig, ax = plt.subplots(figsize=(4,6))
    x = df_list_mean.index.tolist()
    x_1 = [int(i) + Width for i in x]

    plot_inter = df_list_mean["mean_FractionInter"].tolist()
    SD_inter = df_list_mean["SD_mean_FractionInter"].tolist()

    plot_noninter = df_list_mean["mean_FractionNoninter"].tolist()
    SD_noninter = df_list_mean["SD_mean_FractionNoninter"].tolist()

    name = df_list_mean["Dataset"].tolist()

    ax.bar(x, plot_inter, width=Width, color='grey', align='center', yerr=SD_inter,label="interface", alpha=0.6)
    ax.bar(x_1, plot_noninter, width=Width, color='grey', align='center', yerr=SD_noninter, label="non-interface",
           hatch="///", alpha=0.6)
    ax.set_ylabel("mean fraction DI", fontsize=Fontsize)
    ax.set_xticks(x)
    ax.tick_params(axis='y', labelsize=Fontsize)
    ax.set_xticklabels(name, fontsize=Fontsize, rotation=0)
    mpl.rcParams['hatch.linewidth'] = 2.5

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, ncol=1, fontsize=Fontsize,frameon=True)
    plt.gcf().subplots_adjust(bottom=0.18)
    plt.axis('tight')
    plt.grid(False)
    fig.tight_layout()
    plt.savefig(output_basename + m + "_mean_fraction_DI.png")
    plt.savefig(output_basename_pdf + m + "_mean_fraction_DI.pdf")
    plt.close()

    sys.stdout.write('\n~~~~~~~~~~~~                 finished FigZB_07_hardlinked              ~~~~~~~~~~~~')
    sys.stdout.flush()

