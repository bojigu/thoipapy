"""10-fold cross-validation of the trained model.

NOT MAINTAINED. Reachable only from the run_validation pipeline stage, which is switched off in both
shipped settings files and disabled in the functional test. Last exercised for the 2020
publication. Not covered by any test, and excluded from the refactoring applied to the maintained
part of the package.
"""

import os
import pickle
import sys
import time

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import StratifiedGroupKFold

import thoipapy.utils
from thoipapy.artefacts import ArtefactPaths
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters
from thoipapy.utils import dropna_with_report
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML


def run_10fold_cross_validation(
    paths: ArtefactPaths, min_n_homol_training: int, cross_validation_number_of_splits: int, bootstrap: bool, logging
):
    """Run 10-fold cross-validation for a particular set of TMDs (e.g. set04).

    The SAME SET is used for both training and cross-validation.

    The number of folds is determined by "cross_validation_number_of_splits" in the settings file.

    IMPORTANT. CURRENTLY THERE IS NO AUTOMATIC REDUNDANCY CHECK.
     - homologues of the tested protein could be in the training dataset

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Saved Files
    -----------
    crossvalidation_pkl : pickle
        Pickled dictionary (xv_dict) containing the results for each fold of validation.
        Also contains the mean ROC curve, and the mean AUC.
    """
    sys.stdout.write("\n--------------- starting run_10fold_cross_validation ---------------\n")
    train_data_after_first_feature_seln_csv = (
        paths.results_dir / "train_data/03_train_data_after_first_feature_seln.csv"
    )
    tuned_ensemble_parameters_csv = paths.results_dir / "train_data/04_tuned_ensemble_parameters.csv"
    crossvalidation_pkl = os.path.join(paths.crossvalidation_dir, "data", f"{paths.setname}_10F_data.pkl")
    features_csv = paths.results_dir / "feat_imp/test_features.csv"

    thoipapy.utils.make_sure_path_exists(crossvalidation_pkl, isfile=True)
    thoipapy.utils.make_sure_path_exists(features_csv, isfile=True)

    df_data = pd.read_csv(train_data_after_first_feature_seln_csv, index_col=0)
    df_data = dropna_with_report(df_data, "run_10fold_cross_validation", logging)

    # drop training data (full protein) that don't have enough homologues
    if min_n_homol_training != 0:
        df_data = df_data.loc[df_data.n_homologues >= min_n_homol_training]

    X = drop_cols_not_used_in_ML(logging, df_data)
    y = df_data["interface"]

    # Group by protein. StratifiedKFold splits residues, so residues of the same TMD landed in
    # both train and test. The shipped rows happen to be ordered protein-by-protein, which made
    # unshuffled folds roughly contiguous protein blocks and kept the published number honest --
    # but adding shuffle=True, or simply reordering the rows, silently bought about +0.03 AUC.
    # StratifiedGroupKFold makes the grouping explicit instead of accidental.
    groups = df_data["acc_db"] if "acc_db" in df_data.columns else df_data.index.str.split("_").str[0]
    sgkf = StratifiedGroupKFold(n_splits=cross_validation_number_of_splits)
    cv = list(sgkf.split(X, y, groups=groups))

    X.shape[1]
    forest = return_classifier_with_loaded_ensemble_parameters(tuned_ensemble_parameters_csv, bootstrap)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    # save all outputs to a cross-validation dictionary, to be saved as a pickle file
    xv_dict = {}

    start = time.perf_counter()

    for i, (train, test) in enumerate(cv):
        sys.stdout.write(f"f{i + 1}."), sys.stdout.flush()
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1], drop_intermediate=False)
        xv_dict[f"fpr{i}"] = fpr
        xv_dict[f"tpr{i}"] = tpr
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    sys.stdout.write("\n"), sys.stdout.flush()

    logging.info(f"tree depths : {[estimator.tree_.max_depth for estimator in forest.estimators_]}")

    duration = time.perf_counter() - start

    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0

    ROC_AUC = auc(mean_fpr, mean_tpr)

    xv_dict["true_positive_rate_mean"] = mean_tpr
    xv_dict["false_positive_rate_mean"] = mean_fpr
    xv_dict["ROC_AUC"] = ROC_AUC

    # save dict as pickle
    with open(crossvalidation_pkl, "wb") as f:
        pickle.dump(xv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    features_ser = pd.Series(X.columns)
    features_ser.to_csv(features_csv)
    logging.info(
        f"{paths.setname} 10-fold validation. AUC({ROC_AUC:.3f}). Time taken = {duration:.2f}.\nFeatures: {X.columns.tolist()}"
    )
    sys.stdout.write("\n--------------- finished run_10fold_cross_validation ---------------\n")


def create_10fold_cross_validation_fig(paths: ArtefactPaths, cross_validation_number_of_splits: int, logging):
    """Create figure showing ROC curve for each fold in a 10-fold validation.

    The underlying data is created by run_10fold_cross_validation. If this has not been run,
    it will return a file-not-found error.

    Parameters
    ----------
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    sys.stdout.write("\n--------------- starting create_10fold_cross_validation_fig ---------------\n")
    # plt.rcParams.update({'font.size': 7})
    crossvalidation_png = os.path.join(paths.crossvalidation_dir, f"{paths.setname}_10F_ROC.png")
    crossvalidation_pkl = os.path.join(paths.crossvalidation_dir, "data", f"{paths.setname}_10F_data.pkl")

    # open pickle file
    with open(crossvalidation_pkl, "rb") as f:
        xv_dict = pickle.load(f)

    figsize = np.array([3.42, 3.42]) * 2  # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for i in range(cross_validation_number_of_splits):
        roc_auc = auc(xv_dict[f"fpr{i}"], xv_dict[f"tpr{i}"])
        ax.plot(xv_dict[f"fpr{i}"], xv_dict[f"tpr{i}"], lw=1, label=f"fold {i} (area = {roc_auc:0.2f})", alpha=0.8)

    ROC_AUC = xv_dict["ROC_AUC"]

    ax.plot(
        xv_dict["false_positive_rate_mean"],
        xv_dict["true_positive_rate_mean"],
        color="k",
        label=f"mean (area = {ROC_AUC:0.2f})",
        lw=1.5,
    )
    ax.plot([0, 1], [0, 1], "--", color=(0.6, 0.6, 0.6), label="random")
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(crossvalidation_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(crossvalidation_png))
    sys.stdout.write("\n--------------- finished create_10fold_cross_validation_fig ---------------\n")
