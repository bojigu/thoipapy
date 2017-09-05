#
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import sys
# from eccpy.tools import normalise_0_1
# from eccpy.settings import setup_analysed_data_folder, read_settings_file
# from tlabassays.sin_disrupt_fit import fit_sin_to_data, extract_aa_pos_from_sample_names
#
# #from tlabtools.sine_fit import fit_sin_to_data, extract_aa_pos_from_sample_names
#
# # aa_pos_csv = r"D:\DATABASES_spitfire\BlaTM\20151229\20151229_collected.csv"
# # csv containing the scanning mutagenesis data
# # aa_pos_csv = r"D:\DATABASES_spitfire\BlaTM\collected\20160103\20160103_collected.csv"
# aa_pos_csv = r"D:\Schweris\Projects\Programming\Python\files\20151217_blatm_module\analysed\20160105 all TWDT,excl 1010p2_1022p1_1104p2_1201\20160105_analysed.csv"
# dfc_uniq_summ = pd.read_csv(aa_pos_csv, index_col = 0)
#
# # excel file containing the dictionary linking the sample name to a certain amino acid mutation position
# # aa_pos_excel = r"D:\DATABASES_spitfire\BlaTM\SineDisrupt_aa_pos_dict.xlsx"
# aa_pos_excel = r"D:\Schweris\Projects\Programming\Python\files\20151217_blatm_module\SineDisrupt_aa_pos_dict.xlsx"
#
# # define the base-name for the summary figures, including scanning mutagenesis sine-fit output
# # scan_mut_sine_output_basename = r"D:\DATABASES_spitfire\BlaTM\sindisrupt\summary"
# # scan_mut_sine_output_folder = r"D:\DATABASES_spitfire\BlaTM\sindisrupt"
# scan_mut_sine_output_folder = r"D:\Schweris\Projects\Programming\Python\files\20151217_blatm_module\scan_mut_output"
# scan_mut_sine_output_file_basename = "scanmut_sine"
#
# # set up boolean values to determine the method used to extract the amino acid position from the sample names
# extract_aa_pos_from_sample_names = False
# extract_aa_pos_from_sample_name_dict = True
#
# # open the excel file containing the dictionary linking the sample name to a certain amino acid mutation position
# df_aa_pos = pd.read_excel(aa_pos_excel, sheet_name = "aa_pos_dict")
#
# # extract the amino acid position from the sample names
# if extract_aa_pos_from_sample_names:
#     extract_aa_pos_from_sample_names(dfc_uniq_summ, start_aa=233, end_aa=251, newcol_name="aa_pos")
# if extract_aa_pos_from_sample_name_dict:
#     # convert dataframe to a python dictionary
#     aa_pos_dict = dict(zip(df_aa_pos.full_sample_name, df_aa_pos.aa_pos))
#     # use the dictionary to extract the relevant amino acid positions from the sample names
#     dfc_uniq_summ["aa_pos"] = dfc_uniq_summ["longname"].apply(lambda x : aa_pos_dict[x] if x in list(aa_pos_dict.keys()) else np.nan)
#
# # drop any of the columns that do not have a particular amino acid position (for example, positive control samples)
# dfc_uniq_summ.dropna(subset=["aa_pos"], inplace=True)
# # confirm that all of the amino acid positions are integers
# dfc_uniq_summ["aa_pos"] = dfc_uniq_summ["aa_pos"].astype(int)
# # sort the samples by the amino acid position in the dataframe
# dfc_uniq_summ = dfc_uniq_summ.sort_values(by="aa_pos")
# # create a list of unique amino acid positions (for which there could be multiple mutations)
# list_unique_pos = dfc_uniq_summ.aa_pos.dropna().unique()
# # find the min and max of the amino acid positions
# min_aa_pos = list_unique_pos.min()
# max_aa_pos = list_unique_pos.max()
# # use min and max to create a contiguous list of positions, including any gaps
# list_contiguous_aa_pos = list(range(min_aa_pos,max_aa_pos+1))
#
# # set the index as the aa_pos
# dfc_uniq_summ.set_index("aa_pos", inplace=True)
# # create counter for the "gap" positions, where data is not contiguous for all residues
# gapcounter = 0
#
# # iterate through the database types
# datasets = ["_orig", "_ful"]
# for d in datasets:
#     # create a new dataframe to store the average values for all mutations
#     dfc_mult = pd.DataFrame()
#     for pos in list_contiguous_aa_pos:
#         if pos in dfc_uniq_summ.index:
#             # select the row (or rows) with mutations at that position
#             selected_data_at_aa_pos = dfc_uniq_summ.loc[pos,:]
#             # if selected_data_at_aa_pos is a series, there is only one mutation at that position.
#             if isinstance(selected_data_at_aa_pos,pd.Series):
#                 dfc_mult.loc[pos,"mean_mult{}".format(d)] = selected_data_at_aa_pos["mean{}".format(d)]
#                 dfc_mult.loc[pos,"std_mult{}".format(d)] = 0
#             # if selected_data_at_aa_pos is a dataframe, there are multiple mutations at that position. Take the mean.
#             elif isinstance(selected_data_at_aa_pos,pd.DataFrame):
#                 dfc_mult.loc[pos,"mean_mult{}".format(d)] = selected_data_at_aa_pos["mean{}".format(d)].mean()
#                 dfc_mult.loc[pos,"std_mult{}".format(d)] = selected_data_at_aa_pos["mean{}".format(d)].std()
#         else:
#             # the position is not in the index, a gap is counted, the values will automatically be np.nan
#             if gapcounter == 0:
#                 print("gaps without data:")
#             dfc_mult.loc[pos,"mean_mult{}".format(d)] = np.nan
#             dfc_mult.loc[pos,"std_mult{}".format(d)] = 0
#             # print a list of the positions that are not in the index
#             sys.stdout.write("{}, ".format(pos))
#             sys.stdout.flush()
#             gapcounter += 1
#
#     # define the x-axis for the sine calculation as the contiguous index integer
#     x_sin = dfc_mult.index
#     # the y-axis is the mean of all mutations at that position (note, NOT the mean of all replicates)
#     y_sin = normalise_0_1(dfc_mult["mean_mult{}".format(d)])[0]
#     title = "Scanning mutation of MHC II beta TMD, BlaTM heterodimer assay, {}".format(d[1:])
#     # HLA_sine_constants, HLA_periodicity = fit_sin_to_data(x_sin,y_sin, title, scan_mut_sine_output_folder,
#     #                                                       scan_mut_sine_output_file_basename + d)
#     output_tuple = fit_sin_to_data(x_sin, y_sin, title,
#                                    scan_mut_sine_output_folder,
#                                    scan_mut_sine_output_file_basename + d)
#     print(output_tuple)
#     print(scan_mut_sine_output_folder,)
#     print(scan_mut_sine_output_file_basename + d)
#
#     title_shuf = "Scan mut MHC II beta TMD, BlaTM heterodimer assay, RANDOMISED, {}".format(d[1:])
#     y_sin_shuffled = y_sin.copy()
#     np.random.shuffle(y_sin_shuffled)
#     output_tuple = fit_sin_to_data(x_sin, y_sin_shuffled, title_shuf,
#                                    scan_mut_sine_output_folder,
#                                    scan_mut_sine_output_file_basename + "_randomised" + d)
#
# #     list_long_names = ["N(P01920_B2_wt)_C(P01909_A1_wt)_sBla1.1","N(P02724_0_wt)_C(P02724_0_wt)_sBla1.1"]
# #     list_short_names = ["wt","GpA"]
# #     df_select = dfc_uniq_summ.loc[list_long_names,:]
# #     df_select.index = list_short_names
# #     df_select.columns = ["n{}".format(d), "mean_mult{}".format(d), "std_mult{}".format(d), 'pos', 'longname', 'shortname', 'pos_re']
# #     # ['n_ful', 'mean_ful', 'std_ful', 'n_orig', 'mean_orig', 'std_orig',
# #     #    'longname', 'shortname', 'pos', 'sample_name', 'pos_re']
# #     df_select = df_select.reindex(columns = ["mean_mult{}".format(d), "std_mult{}".format(d)])
# #     dfc_mult = pd.concat([dfc_mult,df_select])
# #
# #     #outfig = r"D:\Schweris\Projects\Programming\Python\files\20150920_sBla_96_well_multi_figures\20150920_MT_new_templates_96_well_sBla\quick_summary_mutations.png"
# #     yerr = dfc_mult["std_mult{}".format(d)]
# #     dfc_mult.mean_mult.plot(kind="bar", color = "#0489B1", yerr = [yerr,yerr])
# #     plt.savefig(aa_pos_csv[:-4] + d + "_scanmut.png", dpi=150)
# #
#
#
