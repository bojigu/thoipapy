import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

import thoipapy.utils
from thoipapy.validation.auc import calc_PRAUC_ROCAUC_using_10F_validation
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters
from thoipapy.validation.bocurve import calc_best_overlap_from_selected_column_in_df, calc_best_overlap, parse_BO_data_csv_to_excel


def calc_feat_import_from_mean_decrease_accuracy(s, logging):
    """Calculate feature importances using mean decrease in accuracy.

    This method differs from calc_feat_import_from_mean_decrease_impurity.
    It's much slower, and involves the use of 10-fold cross-validation for each variable separately.

     - a feature (or group of features) is selected for randomisation
     - in theory, randomising important features will cause a drop in prediction accuracy
     - The feature (or group of features) is shuffled
     - precision-recall AUC and ROC-AUC is measured
     - the difference between the original AUC and the AUC with shuffled variable is measured
     - higher values suggest more important features

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    feat_imp_MDA_xlsx : xlsx
        Comma separated values, showing decrease in AUC for each feature or group of features.
    """
    logging.info('------------ starting calc_feat_import_from_mean_decrease_accuracy ------------')
    # input
    train_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
    tuned_ensemble_parameters_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"
    # output
    feat_imp_MDA_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "feat_imp_mean_decrease_accuracy.xlsx")
    feat_imp_temp_THOIPA_BO_curve_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/feat_imp/feat_imp_temp_THOIPA.best_overlap_data.csv"
    feat_imp_temp_bocurve_data_xlsx = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/feat_imp/feat_imp_temp_bocurve_data.xlsx"

    thoipapy.utils.make_sure_path_exists(feat_imp_MDA_xlsx, isfile=True)

    df_data = pd.read_csv(train_data_csv, index_col=0)
    df_data = df_data.dropna()

    # drop training data (full protein) that don't have enough homologues
    if s["min_n_homol_training"] != 0:
        df_data = df_data.loc[df_data.n_homologues >= s["min_n_homol_training"]]

    cols_excluding_y = [c for c in df_data.columns if c != s['bind_column']]
    X = df_data[cols_excluding_y]
    y = df_data["interface"]

    polarity_features = ["test_dropping_of_features_not_included", "polarity", "relative_polarity", "polarity4mean", "polarity3Nmean", "polarity3Cmean", "polarity1mean"]
    pssm_features = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "CS", "DE", "KR", "QN", "LIV"]
    coev_features = ["DImax", "MImax", "DItop4mean", "MItop4mean", "DItop8mean", "MItop8mean", "DI4max", "MI4max", "DI1mean", "MI1mean", "DI3mean", "MI3mean", "DI4mean", "MI4mean", "DI4cum", "MI4cum"]
    DI_features = ["DImax", "DItop4mean", "DItop8mean", "DI4max", "DI1mean", "DI3mean", "DI4mean", "DI4cum"]
    MI_features = ["MImax", "MItop4mean", "MItop8mean", "MI4max", "MI1mean", "MI3mean", "MI4mean", "MI4cum"]
    cons_features = ["entropy", "cons4mean", "rate4site"]
    motif_features =  ["GxxxG", "SmxxxSm"]
    physical_features = ["branched", "mass"]
    TMD_features = ["residue_depth", "n_TMDs", "n_homologues"]
    polarity_and_pssm_features = polarity_features + pssm_features
    features_nested_list = [polarity_and_pssm_features, coev_features, DI_features, MI_features, cons_features, motif_features, physical_features, TMD_features]
    features_nested_namelist = ["polarity_and_pssm_features", "coev_features",  "DI_features", "MI_features", "cons_features", "motif_features", "physical_features", "TMD_features"]

    for i in range(len(features_nested_list)):
        sys.stdout.write("\n{} : {}".format(features_nested_namelist[i], features_nested_list[i]))

    forest = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv)

    pr_auc_orig, roc_auc_orig = calc_PRAUC_ROCAUC_using_10F_validation(X, y, forest)
    auboc_orig = calc_AUBOC_for_feat_imp(y, X, forest, feat_imp_temp_THOIPA_BO_curve_data_csv, feat_imp_temp_bocurve_data_xlsx, logging)

    start = time.clock()

    sys.stdout.write("\nmean : {:.03f}\n".format(pr_auc_orig)), sys.stdout.flush()


    ################### grouped features ###################

    grouped_feat_decrease_PR_AUC_dict = {}
    grouped_feat_decrease_ROC_AUC_dict = {}
    grouped_feat_decrease_AUBOC_dict = {}

    for feature_type, feature_list in zip(features_nested_namelist, features_nested_list):
        feature_list = list(set(feature_list).intersection(set(X.columns.tolist())))
        logging.info("{} : {}".format(feature_type, feature_list))
        X_t = X.copy()
        for feature in feature_list:
            # shuffle the data for that feature
            row_to_shuffle = X_t[feature].to_numpy()
            np.random.shuffle(row_to_shuffle)
            X_t[feature] = row_to_shuffle
        # calculate prediction performance after shuffling
        PR_AUC, ROC_AUC = calc_PRAUC_ROCAUC_using_10F_validation(X_t, y, forest)

        decrease_PR_AUC = pr_auc_orig - PR_AUC
        grouped_feat_decrease_PR_AUC_dict[feature_type] = decrease_PR_AUC

        decrease_ROC_AUC = roc_auc_orig - ROC_AUC
        grouped_feat_decrease_ROC_AUC_dict[feature_type] = decrease_ROC_AUC

        auboc = calc_AUBOC_for_feat_imp(y, X_t, forest, feat_imp_temp_THOIPA_BO_curve_data_csv, feat_imp_temp_bocurve_data_xlsx, logging)

        decrease_auboc = auboc_orig - auboc
        grouped_feat_decrease_AUBOC_dict[feature_type] = decrease_auboc

        logging.info(f"  {feature_type} {decrease_auboc:.03f} | {decrease_PR_AUC:.03f} | {decrease_ROC_AUC:.03f}")


    # remove temp bocurve output files
    feat_imp_temp_THOIPA_BO_curve_data_csv.unlink()
    feat_imp_temp_bocurve_data_xlsx.unlink()

    ################### single features ###################

    single_feat_decrease_PR_AUC_dict = {}
    single_feat_decrease_ROC_AUC_dict = {}
    single_feat_decrease_AUBOC_dict = {}

    for feature in X.columns:
        X_t = X.copy()
        # shuffle the data for that feature
        row_to_shuffle = X_t[feature].to_numpy()
        np.random.shuffle(row_to_shuffle)
        X_t[feature] = row_to_shuffle
        # calculate prediction performance after shuffling
        PR_AUC, ROC_AUC = calc_PRAUC_ROCAUC_using_10F_validation(X_t, y, forest)

        decrease_PR_AUC = pr_auc_orig - PR_AUC
        single_feat_decrease_PR_AUC_dict[feature] = decrease_PR_AUC

        decrease_ROC_AUC = roc_auc_orig - ROC_AUC
        single_feat_decrease_ROC_AUC_dict[feature] = decrease_ROC_AUC

        auboc = calc_AUBOC_for_feat_imp(y, X_t, forest, feat_imp_temp_THOIPA_BO_curve_data_csv, feat_imp_temp_bocurve_data_xlsx, logging)

        decrease_auboc = auboc_orig - auboc
        single_feat_decrease_AUBOC_dict[feature] = decrease_auboc

        logging.info(f"  {feature} {decrease_auboc:.03f} | {decrease_PR_AUC:.03f} | {decrease_ROC_AUC:.03f}")

    df_grouped_feat = pd.DataFrame()
    df_grouped_feat["PR_AUC"] = pd.Series(grouped_feat_decrease_PR_AUC_dict)
    df_grouped_feat["ROC_AUC"] = pd.Series(grouped_feat_decrease_ROC_AUC_dict)
    df_grouped_feat["AUBOC"] = pd.Series(grouped_feat_decrease_AUBOC_dict)
    df_grouped_feat.sort_values("AUBOC", ascending=False, inplace=True)

    df_single_feat = pd.DataFrame()
    df_single_feat["PR_AUC"] = pd.Series(single_feat_decrease_PR_AUC_dict)
    df_single_feat["ROC_AUC"] = pd.Series(single_feat_decrease_ROC_AUC_dict)
    df_single_feat["AUBOC"] = pd.Series(single_feat_decrease_AUBOC_dict)
    df_single_feat.sort_values("AUBOC", ascending=False, inplace=True)

    writer = pd.ExcelWriter(feat_imp_MDA_xlsx)

    df_grouped_feat.to_excel(writer, sheet_name="grouped_feat")
    df_single_feat.to_excel(writer, sheet_name="single_feat")

    writer.save()
    writer.close()

    duration = time.clock() - start

    logging.info('{} calc_feat_import_from_mean_decrease_accuracy. PR_AUC({:.3f}). Time taken = {:.2f}.\nFeatures: {}'.format(s["setname"], pr_auc_orig, duration, X.columns.tolist()))
    logging.info(f'output: ({feat_imp_MDA_xlsx})')
    logging.info('------------ finished calc_feat_import_from_mean_decrease_accuracy ------------')


def calc_AUBOC_for_feat_imp(y, X_t, forest, feat_imp_temp_THOIPA_BO_curve_data_csv, feat_imp_temp_bocurve_data_xlsx, logging):
    THOIPA_BO_data_df = pd.DataFrame()
    acc_db_list = pd.Series(X_t.index).str.split("_").str[0].unique().tolist()
    for acc_db in acc_db_list:
        rows_including_test_tmd = pd.Series(X_t.index).str.contains(acc_db).to_list()
        rows_excluding_test_tmd = [not i for i in rows_including_test_tmd]
        y_test_tmd = y.loc[rows_including_test_tmd]
        y_excluding_test_tmd = y.loc[rows_excluding_test_tmd]
        X_t_test_tmd = X_t.loc[rows_including_test_tmd]
        X_t_excluding_test_tmd = X_t.loc[rows_excluding_test_tmd]
        assert acc_db not in "".join(X_t_excluding_test_tmd.index.to_list())

        probas_ = forest.fit(X_t_excluding_test_tmd, y_excluding_test_tmd).predict_proba(X_t_test_tmd)

        experiment_data = y_test_tmd
        prediction_data = probas_[:, 1]

        THOIPA_BO_single_prot_df = calc_best_overlap(acc_db, experiment_data, prediction_data)

        if THOIPA_BO_data_df.empty:
            THOIPA_BO_data_df = THOIPA_BO_single_prot_df
        else:
            THOIPA_BO_data_df = pd.concat([THOIPA_BO_data_df, THOIPA_BO_single_prot_df], axis=1, join="outer")
    THOIPA_BO_data_df.to_csv(feat_imp_temp_THOIPA_BO_curve_data_csv)
    # THOIPA_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(THOIPA_BO_curve_data_csv, BO_curve_folder, names_excel_path)
    parse_BO_data_csv_to_excel(feat_imp_temp_THOIPA_BO_curve_data_csv, feat_imp_temp_bocurve_data_xlsx, logging)
    df_bocurve = pd.read_excel(feat_imp_temp_bocurve_data_xlsx, sheet_name="mean_o_minus_r", index_col=0)
    df_bocurve = df_bocurve.iloc[:5]
    AUBOC = np.trapz(y=df_bocurve["mean_o_minus_r"], x=df_bocurve.index)
    return AUBOC


def fig_feat_import_from_mean_decrease_accuracy(s, logging):

    feat_imp_MDA_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "feat_imp_mean_decrease_accuracy.xlsx")
    feat_imp_MDA_png = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "feat_imp", "mean_decrease_accuracy_filtered_features.png")

    #df_grouped_feat = pd.read_excel(feat_imp_MDA_xlsx, index_col=0, sheet_name="grouped_feat")
    #df_single_feat = pd.read_excel(feat_imp_MDA_xlsx, index_col=0, sheet_name="single_feat")

    logging.warning("fig_feat_import_from_mean_decrease_accuracy is not yet implemented")

    pass