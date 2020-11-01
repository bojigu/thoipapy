from pathlib import Path
from random import shuffle
from typing import Union

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import precision_recall_curve, auc

from thoipapy import utils
from thoipapy.validation.gather import create_df_with_all_predictions_for_all_residues_in_set


def create_precision_recall_all_residues(s, df_set, logging):
    """Combine all residue predictions, so precision recall can be calculated from a single array.

    Effectively stacks the CSVs on top of each other.
    Code is directly copied and modified from create_ROC_all_residues.

    This is the recommended method for the pr of a dataset, and complements methods where the precision-recall
    is calculated for each protein separately.

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    predictions_csv : csv
        csv file with stacked predictions data for multiple proteins
        index = range(0, ..)
        columns =
    """
    logging.info('Starting combine_all_residue_predictions.')

    # output file with all predictions
    pred_all_res_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_pred_all_res.csv"
    # all_res_precision_recall_data_dict_pkl = os.path.join(s["thoipapy_data_folder"], "results", s["setname"], "precision_recall", "{}_all_res_precision_recall_data_dict.pickle".format(s["setname"]))
    all_res_precision_recall_data_csv: Path = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data.csv"
    all_res_precision_recall_png: Path = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall.png"

    utils.make_sure_path_exists(pred_all_res_csv, isfile=True)

    df_set_nonred = utils.drop_redundant_proteins_from_list(df_set, logging)

    df_all = create_df_with_all_predictions_for_all_residues_in_set(s, df_set_nonred, pred_all_res_csv, logging)

    save_fig_precision_recall_all_residues(s, df_all, all_res_precision_recall_png, all_res_precision_recall_data_csv, logging)

    df_all["subset"] = df_all.acc_db.str.split("-").str[1]

    subsets = ["ETRA", "NMR", "crystal"]
    for subset in subsets:
        df_subset = df_all.loc[df_all.subset == subset]
        if df_subset.empty:
            continue
        precision_recall_png: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data_{subset}_subset.png"
        precision_recall_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data_{subset}_subset.csv"
        save_fig_precision_recall_all_residues(s, df_subset, precision_recall_png, precision_recall_data_csv, logging)

    # with open(all_res_precision_recall_data_pkl, "wb") as f:
    #     pickle.dump(output_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Finished combine_all_residue_predictions.')


def save_fig_precision_recall_all_residues(s, df, all_res_precision_recall_png, all_res_precision_recall_data_csv, logging):
    """Save figure for precision recall plot of all residues joined together.

    Code is directly copied and modified from save_fig_ROC_all_residues

    """
    fontsize = 8
    fig, ax = plt.subplots(figsize=(5, 5))
    THOIPA_predictor = "THOIPA_{}_LOO".format(s["set_number"])
    predictors = [THOIPA_predictor, "TMDOCK", "LIPS_surface_ranked", "PREDDIMER", "random"]

    testsetname, trainsetname = utils.get_testsetname_trainsetname_from_run_settings(s)

    if s["setname"] == testsetname:
        predictors.append(f"thoipa.train{trainsetname}")

    output_dict = {}
    interface_random = df.interface_score.tolist()
    shuffle(interface_random)
    df["random"] = interface_random

    for predictor in predictors:
        df_sel = df[["interface", predictor]].dropna()
        if predictor in ["TMDOCK", "PREDDIMER"]:
            pred = - df_sel[predictor]
            # pred = normalise_between_2_values(df_sel[predictor], 2.5, 8, invert=True)
        else:
            pred = df_sel[predictor]
        precision, recall, thresholds_PRC = precision_recall_curve(df_sel.interface, pred)

        pred_auc = auc(recall, precision)
        # sys.stdout.write("{} AUC : {:.03f}\n".format(predictor, pred_auc))
        label = "{}. AUC : {:.03f}".format(predictor, pred_auc)
        ax.plot(recall, precision, label=label, linewidth=1)

        output_dict[predictor] = {"precision": list(precision), "recall": list(recall), "pred_auc": pred_auc}
    ax.grid(False)

    ax.set_xlabel("recall", fontsize=fontsize)
    ax.set_ylabel("precision", fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(all_res_precision_recall_png, dpi=240)
    fig.savefig(str(all_res_precision_recall_png)[:-4] + ".pdf")

    df_precision_recall_data = pd.DataFrame(output_dict).T
    df_precision_recall_data.to_csv(all_res_precision_recall_data_csv)

    logging.info("save_fig_precision_recall_all_residues finished ({})".format(all_res_precision_recall_data_csv))
