import os
import sys
from pathlib import Path
from typing import List

import pandas as pd

from thoipapy.utils import reorder_dataframe_columns, make_sure_path_exists
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML


def remove_duplicate_features_with_lower_MDI(s, logging):
    logging.info('starting remove_duplicate_features_with_lower_MDI')
    # inputs
    train_data_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_train_data_orig.csv"
    max_similarity_duplicate_features = s["max_similarity_duplicate_features"]
    # outputs
    mean_decrease_impurity_all_features_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_feat_imp_MDI_before_feature_seln.csv"
    results_remove_dup_feat_with_low_MDI_csv = Path(s["data_dir"]) / "results" / s["setname"] / "feat_imp/results_remove_dup_feat_with_low_MDI.csv"
    train_data_excl_duplicates_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/02_train_data_excl_duplicates.csv"
    make_sure_path_exists(results_remove_dup_feat_with_low_MDI_csv, isfile=True)

    df_MDI = pd.read_csv(mean_decrease_impurity_all_features_csv, index_col=0)
    df_data = pd.read_csv(train_data_csv, index_col=0)

    if True in df_data.columns.str.contains("Unnamed").tolist():
        raise ValueError(f"unnamed column found when reading {train_data_csv}")
    df_X = drop_cols_not_used_in_ML(logging, df_data, s["settings_path"])

    correlated_feature_pairs = set()
    correlation_matrix = df_X.corr()

    duplicate_features_to_be_removed: List[str] = []

    features_to_be_retained_during_selection = s['features_to_be_retained_during_selection'].split(",")
    logging.info(f"features_to_be_retained_during_selection : {features_to_be_retained_during_selection}")

    for i in range(len(correlation_matrix.columns)):
        for j in range(i):
            if abs(correlation_matrix.iloc[i, j]) > max_similarity_duplicate_features:

                featurename1 = correlation_matrix.columns[i]
                featurename2 = correlation_matrix.columns[j]

                feature1_mean_decrease_impurity = df_MDI.at[featurename1, "mean_decrease_impurity"]
                feature2_mean_decrease_impurity = df_MDI.at[featurename2, "mean_decrease_impurity"]

                feature_with_lower_MDI = featurename1 if feature1_mean_decrease_impurity <= feature2_mean_decrease_impurity else featurename2

                # OVERRIDE SELECTION IF ANY FEATURES IN SETTINGS LIST TO BE RETAINED (FOR USE IN PLOTS AND COMPARISONS IN THE ARTICLE)
                if featurename1 in features_to_be_retained_during_selection and featurename2 in features_to_be_retained_during_selection:
                    logging.warn(f"features_to_be_retained_during_selection contains two correlated features ({featurename1} and {featurename2})")
                elif featurename1 in features_to_be_retained_during_selection:
                    feature_to_be_removed = featurename2
                elif featurename2 in features_to_be_retained_during_selection:
                    feature_to_be_removed = featurename1
                else:
                    feature_to_be_removed = feature_with_lower_MDI

                duplicate_features_to_be_removed.append(feature_to_be_removed)
                feature_to_be_retained = featurename1 if feature_to_be_removed == featurename2 else featurename2

                MDI_of_feature_with_lower_MDI = "{:.03f}".format(df_MDI.at[feature_with_lower_MDI, "mean_decrease_impurity"])
                MDI_of_feature_with_higher_MDI = "{:.03f}".format(df_MDI.at[feature_to_be_retained, "mean_decrease_impurity"])

                correlated_feature_tuple = (feature_to_be_retained, MDI_of_feature_with_higher_MDI, feature_with_lower_MDI, MDI_of_feature_with_lower_MDI)
                correlated_feature_pairs.add(correlated_feature_tuple)
                logging.info(f"duplicate feature removed: keep '{feature_to_be_retained}' remove '{feature_with_lower_MDI}' R2={correlation_matrix.iloc[i, j]}")

    duplicate_ser = pd.Series()
    duplicate_ser["correlated_feature_pairs"] = correlated_feature_pairs
    duplicate_ser["duplicate_features_to_be_removed"] = duplicate_features_to_be_removed
    num_duplicate_features_with_low_MDI_to_be_removed = len(duplicate_features_to_be_removed)
    duplicate_ser["num_duplicate_features_with_low_MDI_to_be_removed"] = num_duplicate_features_with_low_MDI_to_be_removed

    duplicate_ser.to_csv(results_remove_dup_feat_with_low_MDI_csv)

    logging.info(f"correlated_feature_pairs: {correlated_feature_pairs}")
    logging.info(f"duplicate_features_to_be_removed: {duplicate_features_to_be_removed}")

    cols_to_keep = [c for c in df_X.columns if c not in duplicate_features_to_be_removed]
    df_X_excl_duplicates = df_X.reindex(columns=cols_to_keep, index=df_X.index)
    df_X_excl_duplicates[s["bind_column"]] = df_data[s["bind_column"]]
    if df_X_excl_duplicates.isnull().values.any():
        raise Exception("df_X_excl_duplicates contains nan values")

    df_X_excl_duplicates = reorder_dataframe_columns(df_X_excl_duplicates, ["interface"])
    df_X_excl_duplicates.to_csv(train_data_excl_duplicates_csv)

    sys.stdout.write("----- training data filtered to remove duplicate features ----")
    sys.stdout.write(str(df_X_excl_duplicates.head()))
    sys.stdout.write("--------------------------------------------------------------")

    logging.info(f"num_duplicate_features_with_low_MDI_to_be_removed: {num_duplicate_features_with_low_MDI_to_be_removed}")
    logging.info(f"{df_data.shape[1] - 4} features collected in total, {df_X.shape[1]} included in initial analysis, {df_X_excl_duplicates.shape[1] - 1} remaining after removing duplicates.")
    logging.info(f"output saved to {train_data_excl_duplicates_csv}")
    logging.info("finished remove_duplicate_features_with_lower_MDI")
