import sys

import numpy as np
from sklearn.metrics import precision_recall_curve, auc, roc_curve
from sklearn.model_selection import StratifiedKFold


def calc_PRAUC_ROCAUC_using_10F_validation(X, y, forest):
    """Calculate mean precision-recall and ROC AUC using 10-fold cross-validation.

    Parameters
    ----------
    X : pd.DataFrame
    y : pd.Series
    forest : sklearn.ensemble.ExtraTreesClassifier

    Returns
    -------
    (PR_AUC, ROC_AUC) : tuple
        PR_AUC - float of precision-recall AUC
        ROC_AUC - float of ROC-AUC
    """

    skf = StratifiedKFold(n_splits=10)
    cv = list(skf.split(X, y))
    # save precision-recall AUC of each fold
    pr_auc_list = []
    roc_auc_list = []
    for i, (train, test) in enumerate(cv):
        sys.stdout.write("f{}.".format(i + 1)), sys.stdout.flush()
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])

        precision, recall, thresholds_PRC = precision_recall_curve(y.iloc[test], probas_[:, 1])
        pred_auc = auc(recall, precision)
        pr_auc_list.append(pred_auc)

        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1], drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        roc_auc_list.append(roc_auc)

    PR_AUC = np.mean(pr_auc_list)
    ROC_AUC = np.mean(roc_auc_list)
    return PR_AUC, ROC_AUC
