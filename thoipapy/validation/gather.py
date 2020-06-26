import os
from pathlib import Path
from typing import Union

import pandas as pd


def gather_validation_data_for_figs(s, df_set, logging):
    """Gather output files together to make validation figures easier


    Where's the data?
        FULL TRAINING SET
            BO curve
            PR
            perc_interf_vs_PR_cutoff


    MOVE FILES:
    Bo curve data mult predictors
    file_location = r"/Users/yaoxiao/Dropbox/tm_homodimer_dropbox/THOIPA_data/Results/set08/indiv_validation/indiv_validation_data.xlsx"
    read_file = pd.read_excel(file_location, index_col=0,sheet_name = "BO_o_minus_r")

    PR mult pred
     file_location = r"/Users/yaoxiao/Dropbox/tm_homodimer_dropbox/THOIPA_data/Results/set08/precision_recall/set08_all_res_precision_recall_data.csv"

    perc interface PR
    csv_pathstr = r"/Users/yaoxiao/Dropbox/tm_homodimer_dropbox/THOIPA_data/Results/{}/indiv_validation/perc_interf_vs_PR_cutoff_linechart_data.csv".format(setname)


    /media/sindy/m_data/THOIPA_data/Results/set07/validation_summary/make 2 excel files here - one for blind and one for crossval


    CHANGE TO SET 8, MOVE TO RESULTS/set08 folder "crossvalidation"
    r"/Users/yaoxiao/Dropbox/tm_homodimer_dropbox/THOIPA_data/Results/compare_predictors/set05.{}.4predictors_BOCURVE_linechart.csv".
    """
    logging.info('Starting gather_validation_data.')

    ##########################################  Gather the cross-validation data for the full training dataset ##########################################

    # inputs
    indiv_validation_dir: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/indiv_validation"
    indiv_validation_data_xlsx = indiv_validation_dir / "indiv_validation_data.xlsx"
    all_res_precision_recall_data_csv: Path = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data.csv"
    all_res_ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data.csv"
    perc_interf_vs_PR_cutoff_linechart_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/indiv_validation/bocurve/perc_interf_vs_PR_cutoff_linechart_data.csv"
    # outputs
    crossvalidation_summary_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/cross_blind_validation_summary/crossvalidation_summary.xlsx"
    blindvalidation_summary_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/cross_blind_validation_summary/blindvalidation_summary.xlsx"

    df_bocurve = pd.read_excel(indiv_validation_data_xlsx, index_col=0, sheet_name = "BO_o_minus_r")
    df_pr = pd.read_csv(all_res_precision_recall_data_csv, index_col=0)
    df_roc = pd.read_csv(all_res_ROC_data_csv, index_col=0)
    df_perc_interf_vs_pr = pd.read_csv(perc_interf_vs_PR_cutoff_linechart_data_csv, index_col=0)

    writer = pd.ExcelWriter(crossvalidation_summary_xlsx)
    df_bocurve.to_excel(writer, sheet_name="bocurve")
    df_pr.to_excel(writer, sheet_name="precision_recall")
    df_roc.to_excel(writer, sheet_name="roc")
    df_perc_interf_vs_pr.to_excel(writer, sheet_name="perc_interf_vs_pr")
    writer.save()


    ##########################################  Gather the cross-validation data for each subset of the training dataset ##########################################


    subsets = df_set["database"].unique()
    # for each dataset(e.g. ETRA) separately. Saved in "by_subset" subfolder
    for subset in subsets:

        BOCURVE_linechart_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/compare_predictors/{s['setname']}.{subset}.4predictors_BOCURVE_linechart.csv"
        precision_recall_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/precision_recall/{s['setname']}_all_res_precision_recall_data_{subset}_subset.csv"
        ROC_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/ROC/{s['setname']}_all_res_ROC_data_{subset}_subset.csv"
        #perc_interf_vs_PR_cutoff_linechart_single_database_data_csv = indiv_validation_dir / f"by_subset/{subset}_perc_interf_vs_PR_cutoff_linechart_data.csv"
        perc_interf_vs_PR_cutoff_linechart_single_database_data_csv: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/indiv_validation/bocurve/{subset}_perc_interf_vs_PR_cutoff_linechart_data.csv"

        df_bocurve_subset = pd.read_csv(BOCURVE_linechart_csv, index_col=0)
        df_pr_subset = pd.read_csv(precision_recall_data_csv, index_col=0)
        df_roc_subset = pd.read_csv(ROC_data_csv, index_col=0)
        df_perc_interf_vs_p_subset = pd.read_csv(perc_interf_vs_PR_cutoff_linechart_single_database_data_csv, index_col=0)

        crossvalidation_summary_subset_xlsx: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/cross_blind_validation_summary/crossvalidation_summary_subset_{subset}.xlsx"

        writer = pd.ExcelWriter(crossvalidation_summary_subset_xlsx)
        df_bocurve_subset.to_excel(writer, sheet_name="bocurve")
        df_pr_subset.to_excel(writer, sheet_name="precision_recall")
        df_roc_subset.to_excel(writer, sheet_name="roc")
        df_perc_interf_vs_p_subset.to_excel(writer, sheet_name="perc_interf_vs_pr")
        writer.save()
