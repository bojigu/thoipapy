import sys
import pandas as pd
from matplotlib import pyplot as plt
from plotly import graph_objs as go, plotly
from plotly.tools import FigureFactory as FF
from scipy.interpolate import interp1d
import plotly.plotly as py
import tlabtools.tools as tools
plt.rcParams["font.family"] = "Verdana"
from plotly import figure_factory as FF
import glob
import os
import korbinian
import thoipapy
plt.rcParams["font.family"] = "Verdana"
colour_dict = korbinian.utils.create_colour_lists()
colour_lists = tools.create_colour_lists()
color_thoipa = "k"
color_Lips = colour_dict["TUM_colours"]['TUM5']
blue1 = colour_dict["TUM_colours"]['TUM1']
blue5 = colour_dict["TUM_colours"]['TUM5']
black = "k"
TUMblue = colour_dict["TUM_colours"]['TUMBlue']


def FigZB_18(Fontsize,Width,Size):
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
    thoipapy.utils.setup_biopol_plotly(username=s["plotly_username"], api_key=s["plotly_api_key"])

    sys.stdout.write('\n~~~~~~~~~~~~                 starting FigZB_18              ~~~~~~~~~~~~\n')
    sys.stdout.flush()

    dataset_list = ["crystal","NMR","ETRA"]
    save_name_list=[50,46,47]

    for dataset,mmm in zip(dataset_list,save_name_list):
        Path = glob.glob(r"H:\figs\FigBZ18-PreddimerTmdockComparison\{}\*.csv" .format(dataset))

        for pathfile in Path:
            save_name = pathfile[mmm:-8]
            sys.stdout.write('{}, '.format(save_name))
            sys.stdout.flush()
            output_base_filename = save_name
            save_folder = r"H:\figs\FigBZ18-PreddimerTmdockComparison\NMR"
            settings_folder = os.path.dirname(save_folder)
            output_folder = os.path.join(settings_folder, "heatmap",dataset)
            output_basename = os.path.join(output_folder, output_base_filename)

            list_paths = [output_folder]
            for path in list_paths:
                if not os.path.exists(path):
                    os.makedirs(path)
            # read crustal raw data

            df_row_data = pd.read_csv(pathfile, index_col=0)
            df_row_data.reset_index(inplace=True)
            df_row_data.columns = df_row_data.columns.str.replace('Res_name', 'aa_position')
            df_row_data.columns = df_row_data.columns.str.replace('Res_num', 'orig_aa')
            df_row_data.columns = df_row_data.columns.str.replace('Thoipa', 'THOIPA')
            df_row_data.columns = df_row_data.columns.str.replace('Tmdock', 'TMDOCK')
            df_row_data.columns = df_row_data.columns.str.replace('Preddimer', 'PREDDIMER')
            df_row_data["Closedist"] = -df_row_data["Closedist"]
            df_row_data["TMDOCK"] = -df_row_data["TMDOCK"]
            df_row_data["PREDDIMER"] = -df_row_data["PREDDIMER"]

            # put orignal aa in each bo

            x = df_row_data["orig_aa"].tolist()
            # name the 5 rows
            y = ["THOIPA", "PREDDIMER", "TMDOCK", "-closedist","interface"]

            z = [[float(i) for i in df_row_data["THOIPA"].tolist()], [float(i) for i in df_row_data["PREDDIMER"].tolist()], \
                 [float(i) for i in df_row_data["TMDOCK"].tolist()],
                 [float(i) for i in df_row_data["Closedist"].tolist()], \
                 [float(i) for i in df_row_data["Interface"].tolist()]]

            x_label = df_row_data.index.tolist()
            m = interp1d([min(z[0]), max(z[0])], [0, 1])
            z[0] = m(z[0]).tolist()

            m = interp1d([min(z[1]), max(z[1])], [0, 1])
            z[1] = m(z[1]).tolist()

            m = interp1d([min(z[2]), max(z[2])], [0, 1])
            z[2] = m(z[2]).tolist()

            m = interp1d([min(z[3]), max(z[3])], [0, 1])
            z[3] = m(z[3]).tolist()

            m = interp1d([min(z[4]), max(z[4])], [0, 1])
            z[4] = m(z[4]).tolist()

            # set heat map color, blue
            colorscale = [[0, '#afcdd6'], [1, '#1083a3']]
            # set font color, white and black
            font_colors = ['#efecee', '#3c3636']
            # set annotation
            z_text = [x, x, x, x, x]

            columns = ['Colors']
            df_color = pd.DataFrame(index=x_label, columns=columns)
            df_color = df_color.fillna("black")

            # plot
            fig = FF.create_annotated_heatmap(z, x=x_label, y=y, annotation_text=z_text, colorscale=colorscale,
                                              font_colors=font_colors)

            fig['layout'].update(
                title=save_name,
                xaxis=dict(ticks='', side='top', color="black"),
                font=dict(family='Verdana', size=13),
                yaxis=dict(ticks='', ticksuffix=' ', categoryarray=y, autorange='reversed'),
                width=600,
                height=280,
                margin=go.Margin(
                    l=100,
                    r=30))
            # save figure
            py.image.save_as(fig, output_basename + ".png")
            # save figure as pdf
            #py.image.save_as(fig, output_basename_pdf + ".pdf")

    sys.stdout.write('\n~~~~~~~~~~~~                 finished FigZB_18              ~~~~~~~~~~~~')
    sys.stdout.flush()