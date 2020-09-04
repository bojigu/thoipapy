import os
from pathlib import Path
from typing import Union

import pandas as pd
import thoipapy.utils as utils


def gather_validation_data_for_figs(s, df_set, logging):
    """Gather output files together to make validation figures easier
    """
    logging.info('Starting gather_validation_data.')

    ##########################################  Gather the cross-validation data for the full training dataset ##########################################

    # inputs
    indiv_validation_data_xlsx = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/indiv_validation/bocurve/indiv_validation_data.xlsx"
    all_res_precision_recall_data_csv: Path = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data.csv"
    all_res_ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data.csv"
    perc_interf_vs_PR_cutoff_linechart_data_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/indiv_validation/bocurve/perc_interf_vs_PR_cutoff_linechart_data.csv"
    # outputs
    validation_summary_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/validation_summary/validation_summary.xlsx"
    utils.make_sure_path_exists(validation_summary_xlsx, isfile=True)

    df_bocurve = pd.read_excel(indiv_validation_data_xlsx, index_col=0, sheet_name = "BO_o_minus_r")
    df_pr = pd.read_csv(all_res_precision_recall_data_csv, index_col=0)
    df_roc = pd.read_csv(all_res_ROC_data_csv, index_col=0)
    df_perc_interf_vs_pr = pd.read_csv(perc_interf_vs_PR_cutoff_linechart_data_csv, index_col=0)

    writer = pd.ExcelWriter(validation_summary_xlsx)
    df_bocurve.to_excel(writer, sheet_name="bocurve")
    df_pr.to_excel(writer, sheet_name="precision_recall")
    df_roc.to_excel(writer, sheet_name="roc")
    df_perc_interf_vs_pr.to_excel(writer, sheet_name="perc_interf_vs_pr")
    writer.save()


    ##########################################  Gather the cross-validation data for each subset of the training dataset ##########################################


    subsets = df_set["database"].unique()
    # for each dataset(e.g. ETRA) separately. Saved in "by_subset" subfolder
    for subset in subsets:

        BOCURVE_linechart_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/compare_selected_predictors/data/{s['setname']}.{subset}.4predictors_BOCURVE_linechart.csv"
        precision_recall_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data_{subset}_subset.csv"
        ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data_{subset}_subset.csv"
        #perc_interf_vs_PR_cutoff_linechart_single_database_data_csv = indiv_validation_dir / f"by_subset/{subset}_perc_interf_vs_PR_cutoff_linechart_data.csv"
        perc_interf_vs_PR_cutoff_linechart_single_database_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/crossvalidation/indiv_validation/bocurve/{subset}_perc_interf_vs_PR_cutoff_linechart_data.csv"

        df_bocurve_subset = pd.read_csv(BOCURVE_linechart_csv, index_col=0)
        df_pr_subset = pd.read_csv(precision_recall_data_csv, index_col=0)
        df_roc_subset = pd.read_csv(ROC_data_csv, index_col=0)
        df_perc_interf_vs_p_subset = pd.read_csv(perc_interf_vs_PR_cutoff_linechart_single_database_data_csv, index_col=0)

        crossvalidation_summary_subset_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/validation_summary/validation_summary_subset_{subset}.xlsx"

        writer = pd.ExcelWriter(crossvalidation_summary_subset_xlsx)
        df_bocurve_subset.to_excel(writer, sheet_name="bocurve")
        df_pr_subset.to_excel(writer, sheet_name="precision_recall")
        df_roc_subset.to_excel(writer, sheet_name="roc")
        df_perc_interf_vs_p_subset.to_excel(writer, sheet_name="perc_interf_vs_pr")
        writer.save()


def create_df_with_all_predictions_for_all_residues_in_set(s, df_set_nonred, pred_all_res_csv, logging):
    # set up a dataframe to hold the features for all proteins
    df_all = pd.DataFrame()
    for i in df_set_nonred.index:
        acc = df_set_nonred.loc[i, "acc"]
        database = df_set_nonred.loc[i, "database"]
        merged_data_csv_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/predictions/merged/{database}.{acc}.merged.csv"

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
    df_all = utils.reorder_dataframe_columns(df_all, column_list)
    df_all.to_csv(pred_all_res_csv)
    return df_all