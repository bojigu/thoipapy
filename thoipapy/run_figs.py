import matplotlib.pyplot as plt
import thoipapy
import sys
import eccpy
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

# read the command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "\Path\to\your\settingsfile.xlsx"')

if __name__ == "__main__":

    plt.rcParams["font.family"] = "Verdana"

    #s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\Setting_figure_Bo_win10.xlsx"
    #s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\thoipapy_run_settings_Bo_win10.xlsx"

    # get the command-line arguments
    args = parser.parse_args()
    s_path = args.s

    df_setting_control = pd.read_excel(s_path, sheetname='run_figs', index_col=0)
    df_setting_control.dropna(inplace = True, subset=["Run"])

    df_setting_control["Run"] = df_setting_control["Run"].apply(eccpy.tools.convert_truelike_to_bool, convert_nontrue=False)
    s = df_setting_control.Run.to_dict()
    #sys.stdout.write('{}, '.format(s))
    #sys.stdout.flush()

    thoipapy.utils.setup_biopol_plotly(username=s["plotly_username"], api_key=s["plotly_api_key"])

    Fontsize = s["Fontsize"]
    Filter = s["Filter"]
    Width= s["Width"]
    Size= s["Size"]
    Linewidth= s["Linewidth"]


    if s["FigZB_07"] == True:
        # barcharts of coevolution values for interface and non-interface
        thoipapy.figs.Average_Fraction_DI.FigZB_07(Fontsize, Width, Size, s)


    if s["FigZB_18"] == True:
        # heatmap of prediction from THOIPA, PREDDIMER, TMDOCK
        thoipapy.figs.Preddimer_TMdock_heatmap.FigZB_18(Fontsize,Width,Size)

    if s["pred_interf_single_prot_using_sel_train_datasets"] == True:
        #thoipapy.figs.Create_Bo_Curve_files.Test_Etra(s)
        thoipapy.figs.Create_Bo_Curve_files.pred_interf_single_prot_using_sel_train_datasets(s)

    if s["combine_BO_curve_files"] == True:
        thoipapy.figs.Combine_Bo_Curve_files.combine_BO_curve_files(s)

    if s["fig_plot_BO_curve_mult_train_datasets"] == True:
        thoipapy.figs.Combine_Bo_Curve_files.fig_plot_BO_curve_mult_train_datasets(s)

    if s["fig_plot_BO_curve_mult_predictors"] == True:
        thoipapy.figs.Combine_Bo_Curve_files.fig_plot_BO_curve_mult_predictors(s)

    if s["run_bocurve_comp"] == True:
        thoipapy.figs.BoCurve_ThoipaBest_comp_LIPS_and_Nmr.run_bocurve_comp(Fontsize,Width,Size,s,Linewidth)

    if s["calc_PREDDIMER_TMDOCK_closedist"] == True:
        thoipapy.figs.Calc_PREDDIMER_TMDOCK_Closedist.calc_closedist_from_PREDDIMER_TMDOCK_best_model(s)



