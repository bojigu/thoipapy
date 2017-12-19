print("FIGS.PY HAS BEEN DEPRECATED, and all functions moved to run.py")
#
# import matplotlib.pyplot as plt
# import thoipapy
# import sys
# import eccpy
# import argparse
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)
# import pandas as pd
#
# # read the command line arguments
# parser = argparse.ArgumentParser()
#
# parser.add_argument("-s",  # "-settingsfile",
#                     help=r'Full path to your excel settings file.'
#                          r'E.g. "\Path\to\your\settingsfile.xlsx"')
#
# if __name__ == "__main__":
#
#     plt.rcParams["font.family"] = "Verdana"
#
#     #s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\Setting_figure_Bo_win10.xlsx"
#     #s_path = r"C:\Users\ZENGBO\Dropbox\tm_homodimer_dropbox\thoipapy_run_settings_Bo_win10.xlsx"
#
#     # get the command-line arguments
#     args = parser.parse_args()
#     s_path = args.s
#
#     df_setting_control = pd.read_excel(s_path, sheetname='run_figs', index_col=0)
#     df_setting_control.dropna(inplace = True, subset=["Run"])
#
#     df_setting_control["Run"] = df_setting_control["Run"].apply(korbinian.utils.convert_truelike_to_bool, convert_nontrue=False)
#     s = df_setting_control.Run.to_dict()
#     #sys.stdout.write('{}, '.format(s))
#     #sys.stdout.flush()
#
#



