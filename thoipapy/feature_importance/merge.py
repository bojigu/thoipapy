from ast import literal_eval
from pathlib import Path
from typing import List, Set

import pandas as pd

from thoipapy.utils import reorder_dataframe_columns, make_sure_path_exists

def merge_top_features_anova_ensemble(s, logging):
    logging.info('starting merge_top_features_anova_ensemble')
    # inputs
    train_data_excl_duplicates_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data_excl_duplicates.csv"
    top_features_anova_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/feat_imp/top_features_anova.csv"
    top_features_rfe_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/feat_imp/top_features_rfe.csv"
    # outputs
    train_data_excl_dup_top_feat = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/train_data_excl_dup_top_feat.csv"

    make_sure_path_exists(train_data_excl_dup_top_feat, isfile=True)

    anova_ser = pd.read_csv(top_features_anova_csv, index_col=0).iloc[:, 0]
    df_rfe = pd.read_csv(top_features_rfe_csv, index_col=0)

    anova_top_features: List[str] = literal_eval(anova_ser["top_features"])
    rfe_top_features: List[str] = df_rfe.loc[df_rfe["ranking"] == 1]["features"].to_list()

    combined_top_features: List[str] = list(set(anova_top_features + rfe_top_features))
    n_combined_top_features: int = len(combined_top_features)
    features_in_anova_but_not_ensemble_rfe: Set[str] = set(anova_top_features) - set(rfe_top_features)
    features_in_ensemble_rfe_but_not_anova: Set[str] = set(rfe_top_features) - set(anova_top_features)

    df_train_data_excl_duplicates = pd.read_csv(train_data_excl_duplicates_csv)
    combined_top_features_incl_y: List[str] = [s["bind_column"]] + combined_top_features
    n_dropped_features: int = df_train_data_excl_duplicates.shape[1] - n_combined_top_features

    logging.info(f"n_combined_top_features : {n_combined_top_features}")
    logging.info(f"combined_top_features : {combined_top_features}")
    logging.info(f"features_in_anova_but_not_ensemble_rfe : {features_in_anova_but_not_ensemble_rfe}")
    logging.info(f"features_in_ensemble_rfe_but_not_anova : {features_in_ensemble_rfe_but_not_anova}")
    logging.info(f"n_dropped_features : {n_dropped_features}")

    df_train_data_excl_dup_top_feat = df_train_data_excl_duplicates.reindex(columns=combined_top_features_incl_y, index=df_train_data_excl_duplicates.index)
    df_train_data_excl_dup_top_feat.to_csv(train_data_excl_dup_top_feat)

    logging.info('finished merge_top_features_anova_ensemble')
