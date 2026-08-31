from ast import literal_eval

import pandas as pd

from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import make_sure_path_exists


def merge_top_features_anova_ensemble(
    paths: ArtefactPaths, bind_column: str, features_to_be_retained_during_selection, logging
):
    logging.info("starting merge_top_features_anova_ensemble")
    # inputs
    train_data_excl_duplicates_csv = paths.train_data_excl_duplicates_csv()
    top_features_anova_csv = paths.top_features_anova_csv()
    top_features_rfe_csv = paths.top_features_rfe_csv()
    # outputs
    train_data_after_first_feature_seln_csv = paths.train_data_after_first_feature_seln_csv()

    make_sure_path_exists(train_data_after_first_feature_seln_csv, isfile=True)

    anova_ser = pd.read_csv(top_features_anova_csv, index_col=0).iloc[:, 0]
    df_rfe = pd.read_csv(top_features_rfe_csv, index_col=0)
    features_to_be_retained_during_selection = features_to_be_retained_during_selection.split(",")

    anova_top_features: list[str] = literal_eval(anova_ser["top_features"])
    rfe_top_features: list[str] = df_rfe.loc[df_rfe["ranking"] == 1]["features"].to_list()

    # dict.fromkeys deduplicates while preserving first-seen order; list(set(...)) did not, which
    # made the selected feature set itself vary between runs.
    combined_top_features: list[str] = list(
        dict.fromkeys(anova_top_features + rfe_top_features + features_to_be_retained_during_selection)
    )
    n_combined_top_features: int = len(combined_top_features)
    features_in_anova_but_not_ensemble_rfe: set[str] = set(anova_top_features) - set(rfe_top_features)
    features_in_ensemble_rfe_but_not_anova: set[str] = set(rfe_top_features) - set(anova_top_features)

    df_train_data_excl_duplicates = pd.read_csv(train_data_excl_duplicates_csv, index_col=0)
    combined_top_features_incl_y: list[str] = [bind_column] + combined_top_features
    n_dropped_features: int = df_train_data_excl_duplicates.shape[1] - n_combined_top_features

    logging.info(f"n_combined_top_features : {n_combined_top_features}")
    logging.info(f"combined_top_features : {combined_top_features}")
    logging.info(f"features_in_anova_but_not_ensemble_rfe : {features_in_anova_but_not_ensemble_rfe}")
    logging.info(f"features_in_ensemble_rfe_but_not_anova : {features_in_ensemble_rfe_but_not_anova}")
    logging.info(f"n_dropped_features : {n_dropped_features}")
    logging.info(f"total number of retained features : {len(combined_top_features)}")

    for column_name in combined_top_features_incl_y:
        if column_name not in df_train_data_excl_duplicates.columns:
            raise Exception(f"df_train_data_excl_duplicates does not contain {column_name}")

    df_train_data_excl_dup_top_feat = df_train_data_excl_duplicates.reindex(
        columns=combined_top_features_incl_y, index=df_train_data_excl_duplicates.index
    )
    df_train_data_excl_dup_top_feat.to_csv(train_data_after_first_feature_seln_csv)

    logging.info("finished merge_top_features_anova_ensemble")
