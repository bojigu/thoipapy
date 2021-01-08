from pathlib import Path
from typing import List

import pandas as pd
from thoipapy.utils import make_sure_path_exists
from sklearn.feature_selection import RFE
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters


def select_best_features_with_ensemble_rfe(s, logging):
    """ The f_classif function is an rfe implementation in python.

    This function selects the top features according to the rfe analysis.
    """

    logging.info('starting select_best_features_with_rfe')
    # inputs
    tuned_ensemble_parameters_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_tuned_ensemble_parameters_before_feature_seln.csv"
    train_data_excl_duplicates_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/02_train_data_excl_duplicates.csv"
    # outputs
    top_features_rfe_csv = Path(s["data_dir"]) / f"results/{s['setname']}/feat_imp/top_features_rfe.csv"

    make_sure_path_exists(top_features_rfe_csv, isfile=True)

    df_data = pd.read_csv(train_data_excl_duplicates_csv, index_col=0)
    if True in df_data.columns.str.contains("Unnamed").tolist():
        raise ValueError(f"unnamed column found when reading {train_data_excl_duplicates_csv}")

    X = df_data.copy()
    del X[s["bind_column"]]
    assert "interface" not in X.columns

    y = df_data[s["bind_column"]]

    forest = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv)
    rfe = RFE(forest, s["n_top_features_to_keep"])
    fit = rfe.fit(X, y)
    logging.info("Num Features: %d" % fit.n_features_)
    logging.info("Selected Features: %s" % fit.support_)
    logging.info("Feature Ranking: %s" % fit.ranking_)

    df_rfe = pd.DataFrame()
    df_rfe["features"] = X.columns
    df_rfe["ranking"] = fit.ranking_
    df_rfe.sort_values("ranking", ascending=True, inplace=True)
    df_rfe.to_csv(top_features_rfe_csv)

    top_features_ensemble_rfe: List[str] = df_rfe.loc[df_rfe["ranking"] == 1]["features"].to_list()
    logging.info(f"best {s['n_top_features_to_keep']} according to ensemble RFE : {top_features_ensemble_rfe}")
    logging.info(f"output saved to {top_features_rfe_csv}")
    logging.info('finished select_best_features_with_ensemble_rfe')

    # logging.info("----------------------------------------------------------------------")
    # logging.info("starting linear RFE")
    #
    # model = LogisticRegression(solver="lbfgs")
    # rfe = RFE(model, 5)
    # fit = rfe.fit(X, y)
    # logging.info("Num Features: %d" % fit.n_features_)
    # logging.info("Selected Features: %s" % fit.support_)
    # logging.info("Feature Ranking: %s" % fit.ranking_)
    # df_lin_ranking = pd.DataFrame()
    # df_lin_ranking["features"] = X.columns
    # df_lin_ranking["ranking"] = fit.ranking_
    #
    # df_lin_ranking.sort_values("ranking", ascending=True, inplace=True)
    # logging.info(df_lin_ranking)
    #
    #
    # logging.info("----------------------------------------------------------------------")
    # logging.info("starting PCA RFE")
    #
    # pca = PCA(n_components=5)
    # fit = pca.fit(X, y)
    #
    # logging.info("Explained Variance: %s" % fit.explained_variance_ratio_)
    # logging.info(fit.components_)
    # X_selected = fit.components_
    #
    # for selected_col in X_selected.T:
    #     for orig_featurename in X.columns:
    #         if list(selected_col) == X[orig_featurename].to_list():
    #             logging.info(orig_featurename, "is a PCA hit")
    #
