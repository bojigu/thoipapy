from pathlib import Path

import pandas as pd
from thoipapy.utils import make_sure_path_exists
from sklearn.feature_selection import SelectKBest, f_classif


def select_best_features_with_anova(s, logging):
    """ The f_classif function is an ANOVA implementation in python.

    This function selects the top features according to the ANOVA analysis.
    """
    logging.info('starting select_best_features_with_ANOVA')
    # inputs
    train_data_excl_duplicates_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/02_train_data_excl_duplicates.csv"
    # outputs
    top_features_anova_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/feat_imp/top_features_anova.csv"

    make_sure_path_exists(top_features_anova_csv, isfile=True)

    df_data = pd.read_csv(train_data_excl_duplicates_csv, index_col=0)
    if True in df_data.columns.str.contains("Unnamed").tolist():
        raise ValueError(f"unnamed column found when reading {train_data_excl_duplicates_csv}")

    X = df_data.copy()
    del X[s["bind_column"]]
    assert "interface" not in X.columns

    y = df_data[s["bind_column"]]

    cls = SelectKBest(score_func=f_classif, k=s["n_top_features_to_keep"])
    fit: SelectKBest = cls.fit(X, y)

    X_selected = fit.transform(X)

    top_features_anova = []
    for orig_featurename in X.columns:
        for selected_col in X_selected.T:
            if list(selected_col) == X[orig_featurename].to_list():
                top_features_anova.append(orig_featurename)
                break

    # make sure that there are no duplicate columns that ruin the mapping to column names
    if len(top_features_anova) != len(set(top_features_anova)):
        raise Exception(f"top_features_anova contains duplicate values")

    top_features_anova_ser = pd.Series()
    top_features_anova_ser["top_features"] = top_features_anova
    top_features_anova_ser["n_top_features_anova"] = len(top_features_anova)
    top_features_anova_ser["n_features_excluded"] = X.shape[1] - len(top_features_anova)
    top_features_anova_ser.to_csv(top_features_anova_csv)

    logging.info(f"output saved to {top_features_anova_csv}")
    logging.info('finished select_best_features_with_ANOVA')
