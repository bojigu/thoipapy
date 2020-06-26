import os
from pathlib import Path
from typing import Union

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import roc_curve, auc

import thoipapy.utils


def create_ROC_all_residues(s, df_set, logging):
    """Combine all residue predictions, so AUC can be calculated from a single array.

    Effectively stacks the CSVs on top of each other.

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
    pred_all_res_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_pred_all_res.csv"
    #all_res_ROC_data_dict_pkl: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data_dict.pickle"
    all_res_ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data.csv"
    all_res_ROC_png: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC.png"

    thoipapy.utils.make_sure_path_exists(pred_all_res_csv, isfile=True)

    df_set_nonred = thoipapy.utils.drop_redundant_proteins_from_list(df_set, logging)

    # set up a dataframe to hold the features for all proteins
    df_all = pd.DataFrame()
    for i in df_set_nonred.index:
        acc = df_set_nonred.loc[i, "acc"]
        database = df_set_nonred.loc[i, "database"]
        merged_data_csv_path = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "predictions", database, "{}.merged.csv".format(acc))

        df_merged_new_protein = pd.read_csv(merged_data_csv_path, index_col=0)
        df_merged_new_protein["acc_db"] = "{}-{}".format(acc, database)
#
        # for the first protein, replace the empty dataframe
        if df_all.empty:
            df_all = df_merged_new_protein
        else:
            # concatenate the growing dataframe of combined proteins and new dataframe
            df_all = pd.concat([df_all, df_merged_new_protein])

    # drop any positions where there is no interface_score (e.g. no mutations, or hetero contacts?)
    if "interface_score" in df_all.columns:
        df_all.dropna(subset=["interface_score"], inplace=True)
    else:
        logging.warning("No experimental data has been added to this dataset. Hope you're not trying to train with it!!!")

    # reset the index to be a range (0,...).
    df_all.index = range(df_all.shape[0])

    # reorder the columns
    column_list = ['acc_db', 'interface', 'interface_score', 'residue_num', 'residue_name']
    df_all = thoipapy.utils.reorder_dataframe_columns(df_all, column_list)

    df_all.to_csv(pred_all_res_csv)

    save_fig_ROC_all_residues(s, df_all, all_res_ROC_png, all_res_ROC_data_csv, logging)

    df_all["subset"] = df_all.acc_db.str.split("-").str[1]

    subsets = ["ETRA", "NMR", "crystal"]
    for subset in subsets:
        df_subset = df_all.loc[df_all.subset == subset]
        if df_subset.empty:
            continue
        ROC_png: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data_{subset}_subset.png"
        ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data_{subset}_subset.csv"
        save_fig_ROC_all_residues(s, df_subset, ROC_png, ROC_data_csv, logging)

    # with open(all_res_ROC_data_pkl, "wb") as f:
    #     pickle.dump(output_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Finished combine_all_residue_predictions.')


def save_fig_ROC_all_residues(s, df, all_res_ROC_png, all_res_ROC_data_csv, logging):
    fontsize=8
    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot([0, 1], [0, 1], color="0.5", linestyle="--", label="random", linewidth=1)
    THOIPA_predictor = "THOIPA_{}_LOO".format(s["set_number"])
    predictors = [THOIPA_predictor, "TMDOCK", "LIPS_surface_ranked", "PREDDIMER"]
    output_dict = {}
    for predictor in predictors:
        df_sel = df[["interface", predictor]].dropna()
        if predictor in ["TMDOCK", "PREDDIMER"]:
            pred = - df_sel[predictor]
            # pred = normalise_between_2_values(df_sel[predictor], 2.5, 8, invert=True)
        else:
            pred = df_sel[predictor]
        fpr, tpr, thresholds = roc_curve(df_sel.interface, pred, drop_intermediate=False)
        pred_auc = auc(fpr, tpr)
        #sys.stdout.write("{} AUC : {:.03f}\n".format(predictor, pred_auc))
        label = "{}. AUC : {:.03f}".format(predictor, pred_auc)
        ax.plot(fpr, tpr, label=label, linewidth=1)
        # output_dict["fpr_{}".format(predictor)] = fpr
        # output_dict["tpr_{}".format(predictor)] = tpr
        # output_dict["auc_{}".format(predictor)] = auc

        output_dict[predictor] = {"fpr" : list(fpr), "tpr" : list(tpr), "pred_auc" : pred_auc}
    ax.grid(False)
    ax.set_xlabel("false positive rate", fontsize=fontsize)
    ax.set_ylabel("true positive rate", fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(all_res_ROC_png, dpi=240)
    fig.savefig(str(all_res_ROC_png)[:-4] + ".pdf")

    df_ROC_data = pd.DataFrame(output_dict).T
    df_ROC_data.to_csv(all_res_ROC_data_csv)

    logging.info("save_fig_ROC_all_residues finished ({})".format(all_res_ROC_data_csv))