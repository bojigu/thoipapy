import os
from collections import Counter
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["font.family"] = "Verdana"
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
plt.rcParams["font.family"] = "Verdana"
import matplotlib.ticker as ticker

dropbox_dir = r"D:\Dropbox\tm_homodimer_dropbox"
drive_dir = r"D:\drive\TMD_homodimer"
thoipa_dir = r"D:\data_thoipapy"

out_dir = os.path.join(drive_dir, "figs\\FigXY25_residues_composition")
# excel ="I:\sets\set05_ETRA_NMR_crystal_nr.xlsx"
excel = os.path.join(dropbox_dir, "sets\\set05_ETRA_NMR_crystal_nr.xlsx")

list_aa = []
list_all = []
Size = 10
Fontsize = 8
Width = 0.3

purpose = "paper"

if purpose == "thesis":
    color_list = ['grey', 'grey', 'k']

else:
    color_list = ['grey', '#559ed5', '#0065BD']

#excel = "I:\sets\set05_ETRA_NMR_crystal_nr.xlsx"
df_set = pd.read_excel(excel, index_col=0)

df_set.reset_index(inplace=True)

for nn in df_set.index:
    sour = df_set.loc[nn, "database"]
    uniprot_acc = df_set.loc[nn, "acc"]

    sys.stdout.write('{}, '.format(uniprot_acc))
    sys.stdout.flush()

    # read for crystal  dataset.

    # read combined data which include conservation, covaritation, lipophilicity and the interface from structure.
    Prediction_path = os.path.join(thoipa_dir, "Features\\combined\\{}\\{}.surr20.gaps5.combined_features.csv".format(sour, uniprot_acc))
    prediction = pd.read_csv(Prediction_path, index_col=0)

    prediction.columns = prediction.columns.str.replace('residue_num', 'aa_position')
    prediction.columns = prediction.columns.str.replace('residue_name', 'orig_aa')
    prediction['orig_aa'].replace('Q', 'B', inplace=True)
    prediction['orig_aa'].replace('N', 'B', inplace=True)
    prediction['orig_aa'].replace('H', 'B', inplace=True)
    prediction['orig_aa'].replace('D', 'B', inplace=True)
    prediction['orig_aa'].replace('E', 'B', inplace=True)
    prediction['orig_aa'].replace('K', 'B', inplace=True)
    prediction['orig_aa'].replace('R', 'B', inplace=True)

    # collect interfac eresidue
    list_aaa = []
    for n in prediction.index:
        if prediction.loc[n, "interface"] == 1:
            Orig_aa = prediction.loc[n, "orig_aa"]

            list_aaa.extend(Orig_aa)
    list_aa.extend(list_aaa)

    # collect all residue
    list_all.extend(prediction["orig_aa"].tolist())

# for all residues in 3 datasets.


# calculate total residue numbers for 3 datast seperate
# ETRA
total_number_interface = len(list_aa)
total_number_non_interface = len(list_all) - len(list_aa)

b = Counter(list_all)
c = Counter(list_aa)

Index_all = np.arange(len(b.items()))
Index = np.arange(len(c.items()))

Columns_all = ['Residue_all', 'Frequence_all']
Columns = ['Residue', 'Frequence']

df_count = pd.DataFrame(index=Index, columns=Columns)
n = 0
for key, count in c.items():
    df_count['Residue'][n] = key
    df_count['Frequence'][n] = count
    n += 1

# all residues
df_count_all = pd.DataFrame(index=Index_all, columns=Columns_all)
n = 0
for key, count in b.items():
    df_count_all['Residue_all'][n] = key
    df_count_all['Frequence_all'][n] = count
    n += 1

df_count.sort_values(by="Residue", ascending=True, inplace=True)
df_count_all.sort_values(by="Residue_all", ascending=True, inplace=True)

df_count.reset_index(drop=True, inplace=True)
df_count_all.reset_index(drop=True, inplace=True)

# calculate propensity for non interface
result = pd.merge(df_count_all, df_count, how='left', left_on='Residue_all', right_on='Residue')

interface_noninterface = result.copy()
interface_noninterface = interface_noninterface.fillna(0)

for n in interface_noninterface.index:
    all = interface_noninterface.loc[n, "Frequence_all"]
    interface = interface_noninterface.loc[n, "Frequence"]
    interface_noninterface.loc[n, "Frequence_noninterface"] = all - interface

# calculate use devide
for n in result.index:
    all = result.loc[n, "Frequence_all"]
    interface = result.loc[n, "Frequence"]
    non = all - interface
    if non != 0:
        result.loc[n, "Frequene_nocn"] = non
    else:
        result.loc[n, "Frequene_nocn"] = np.nan

    result.loc[n, "enrichment"] = (interface / total_number_interface) / (all / (total_number_interface + total_number_non_interface))

result.sort_values(by="Residue_all", inplace=True, ascending=True)
result.reset_index(drop=True, inplace=True)

#####
# for 3 dataset
#############################################################################################################
##                                                                                                     ######
##                                       #plot for enrichment                                          ######
##                                                                                                     ######
#############################################################################################################

Columns_plot_enrichment = ['Residue', 'enrichment_str', 'enrichment_NMR', 'enrichment_ETRA', 'total_number_str', 'total_number_NMR', 'total_number_ETRA']

index_all_resides = (result["Residue_all"].tolist())

unique_index_all_resides = np.unique(index_all_resides)

df_plot_enrichment = pd.DataFrame(columns=Columns_plot_enrichment)
df_plot_enrichment["Residue"] = unique_index_all_resides

for i in df_plot_enrichment.index:
    AA = df_plot_enrichment.loc[i, "Residue"]
    for m in result.index:
        AA_ERTA = result.loc[m, "Residue_all"]
        if AA == AA_ERTA:
            df_plot_enrichment.loc[i, "enrichment_ETRA"] = result.loc[m, "enrichment"]
            df_plot_enrichment.loc[i, "total_number_ETRA"] = result.loc[m, "Frequence_all"]

# reorder by enrichment of ETRA
df_plot_enrichment.sort_values(by='total_number_ETRA', ascending=False, inplace=True)
df_plot_enrichment.reset_index(inplace=True, drop=True)

# df_plot_enrichment['Residue'].replace('B', 'strongly polar\n (EDKRQNH)', inplace=True)
df_plot_enrichment['Residue'].replace('B', 'sp', inplace=True)

# reorder by total_number_ETRA of ETRA
df_plot_enrichment.sort_values(by='total_number_ETRA', ascending=False, inplace=True)
df_plot_enrichment.reset_index(inplace=True, drop=True)
out_csv = os.path.join(out_dir, "data\\aa_enrichment_at_interface.csv")
df_plot_enrichment.to_csv(out_csv)

y_etra_enrichment = df_plot_enrichment["enrichment_ETRA"]

# plot for frequence bar chart
fig, ax = plt.subplots(figsize=(3.42, 3.42))
x = df_plot_enrichment.index.tolist()
x_1 = [i + Width for i in x]
x_label = [i + Width * 0.5 for i in x]

y_etra = df_plot_enrichment["total_number_ETRA"].tolist()

df_xticklabel = df_plot_enrichment["Residue"].tolist()

len_ETRA = len(list_all)
# devide by the total residue number

y_etra = [i / len_ETRA for i in y_etra]

ax.bar(x, y_etra, width=Width, color=color_list[0], align='center', label='residue propensity in TMD')

ax2 = ax.twinx()

ax2.bar(x_1, y_etra_enrichment, width=Width, color=color_list[2], align='center', label='residue enrichment at interface')

# 0065BD
plt.rcParams['xtick.labelsize'] = Fontsize
plt.rcParams['ytick.labelsize'] = Fontsize
ax.tick_params(axis='y', labelsize=Fontsize)
ax.tick_params(axis='x', labelsize=Fontsize)
ax.set_ylabel('overall residue propensity in TMD', fontsize=Fontsize, color=color_list[0])
ax2.set_ylabel('residue enrichment at interface', fontsize=Fontsize, color=color_list[2])
ax.set_xticks(x_label)
ax.set_xticklabels(x, fontsize=Fontsize, rotation=0)
plt.margins(0.035, 0.1)
plt.tight_layout()
plt.grid(False)
ax.set_facecolor('white')
fig.patch.set_visible(True)
plt.rcParams['xtick.labelsize'] = Fontsize
plt.rcParams['ytick.labelsize'] = Fontsize
ax.tick_params(axis='y', labelsize=Fontsize, pad=2)
ax.tick_params(axis='x', labelsize=Fontsize, pad=2.2)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

ax2.tick_params(axis='y', labelsize=Fontsize, pad=2)
ax2.tick_params(axis='x', labelsize=Fontsize, pad=2)
ax2.tick_params(direction='out', length=0, width=1, colors='k')

ax.plot([-1, 15], [0.1352, 0.1352], linestyle="--", color=color_list[2], linewidth=0.8, zorder=0)
# ax2.plot([-1,15],[1, 1],linestyle="--",color='red', linewidth=0.8, zorder=0)

plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.linewidth"] = 0.4
plt.grid(False)
ax.tick_params(direction='out', length=0, width=1, colors='k')

ax.set_xticklabels(df_xticklabel, fontsize=Fontsize, rotation=0)

ax.tick_params("y", colors=color_list[0], labelsize=Fontsize)
ax2.tick_params("y", colors=color_list[2], labelsize=Fontsize)

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines + lines2, labels + labels2, loc=0,frameon=True,facecolor='white', fontsize=Fontsize)

plt.text(16.64, 0.4, "under represented", fontsize=Fontsize - 1, rotation=90, verticalalignment='center',
         horizontalalignment='center', color=color_list[2])

plt.text(16.64, 1.24, "over represented", fontsize=Fontsize - 1, rotation=90, verticalalignment='center',
         horizontalalignment='center', color=color_list[2])
# ax.legend(handle, label, ncol=1, loc=1, fontsize=Fontsize,frameon=True,facecolor='white')
# ax2.legend(handle, label, ncol=1, loc=1, fontsize=Fontsize,frameon=True,facecolor='white')

plt.ylim(0, 1.59)
plt.xlim(-0.6, 14)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "figures\\frequence_bar_chart_complete_plot_in_one_figure_{}.pdf".format(purpose)))
plt.savefig(os.path.join(out_dir, "figures\\frequence_bar_chart_complete_plot_in_one_figure_{}.png".format(purpose)), dpi=300)
plt.show()
