import matplotlib.pyplot as plt
import pandas as pd
import thoipapy
import sys
import eccpy


if __name__ == "__main__":

    plt.rcParams["font.family"] = "Verdana"

    #s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\Setting_figure_Bo_win10.xlsx"
    s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\thoipapy_run_settings_Bo_win10.xlsx"

    df_setting_control = pd.read_excel(s_path, sheetname='run', index_col=0)
    df_setting_control.dropna(inplace = True, subset=["Run"])

    df_setting_control["Run"] = df_setting_control["Run"].apply(eccpy.tools.convert_truelike_to_bool, convert_nontrue=False)
    s = df_setting_control.Run.to_dict()
    sys.stdout.write('{}, '.format(s))
    sys.stdout.flush()

    thoipapy.utils.setup_biopol_plotly(username=s["plotly_username"], api_key=s["plotly_api_key"])

    Fontsize = s["Fontsize"]
    Filter = s["Filter"]
    Width= s["Width"]
    Size= s["Size"]
    Linewidth= s["Linewidth"]


    if s["FigZB_05"] == True:
        thoipapy.figs.Average_Fraction_DI.FigZB_05(Fontsize,Width,Size,s)


    if s["FigZB_18"] == True:
        thoipapy.figs.Preddimer_TMdock_heatmap.FigZB_18(Fontsize,Width,Size)


    if s["run_bocurve_comp"] == True:
        thoipapy.figs.BoCurve_ThoipaBest_comp_LIPS_and_Nmr.run_bocurve_comp(Fontsize,Width,Size,s,Linewidth)

    if s["create_bocurve_files"] == True:
        thoipapy.figs.Create_Bo_Curve_files.Test_Etra(s)

