import pandas as pd
import numpy as np
import scipy.stats as stats
import os
from pathlib import Path

from thoipapy.utils import make_sure_path_exists, get_testsetname_trainsetname_from_run_settings


def generate_boot_matrix(z, B):
    """Returns all bootstrap samples in a matrix
    (original source and licence uncertain)
    """

    sample_size = len(z)  # sample size
    idz = np.random.randint(0, sample_size, size=(B, sample_size))  # indices to pick for all boostrap samples
    return z[idz]


def calc_ttest_pvalue_from_bootstrapped_data(x, y, equal_var=False, B=10000, plot=False):
    """Calculate Student's two-sample t-test from bootstrapped data

    Returns bootstrap p-value, test statistics and parametric p-value

    (original source and licence uncertain)
    """

    # Original t-test statistic
    orig = stats.ttest_ind(x, y, equal_var=equal_var)

    # Generate bootstrap distribution of t statistic
    xboot = generate_boot_matrix(x - np.mean(x), B=B)  # important centering step to get sampling distribution under the null
    yboot = generate_boot_matrix(y - np.mean(y), B=B)
    sampling_distribution = stats.ttest_ind(xboot, yboot, axis=1, equal_var=equal_var)[0]

    # Calculate proportion of bootstrap samples with at least as strong evidence against null
    p = np.mean(sampling_distribution >= orig[0])

    return 2 * min(p, 1 - p)

def conduct_ttest_for_selected_features_used_in_model(s, logging):
    logging.info('starting conduct_ttest_for_selected_features_used_in_model')
    testsetname, trainsetname = get_testsetname_trainsetname_from_run_settings(s)
    # inputs
    train_data_csv = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/01_train_data_orig.csv"
    # IMPORTANT: the features used in the model are taken from the trainset as defined in "train_datasets" in the excel file, not from the actual set being investigated
    feat_imp_MDA_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", trainsetname, "feat_imp", "feat_imp_mean_decrease_accuracy.xlsx")

    # outputs
    "/media/sindy/m_data/THOIPA_data/Results/set08/feat_imp/feat_imp_mean_decrease_accuracy.xlsx"
    ttest_pvalues_bootstrapped_data_using_traindata_selected_features_xlsx = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/ttest/ttest_pvalues_bootstrapped_data_using_traindata_selected_features.xlsx"

    # parameters
    bootstrap_repetitions = 100000

    make_sure_path_exists(ttest_pvalues_bootstrapped_data_using_traindata_selected_features_xlsx, isfile=True)

    df = pd.read_csv(train_data_csv, index_col=0)

    df_feat_imp_MDA_trainset = pd.read_excel(feat_imp_MDA_xlsx, sheet_name="single_feat", index_col=0)
    cols = list(df_feat_imp_MDA_trainset.index)

    dfs0 = df[df.interface == 0]
    dfs1 = df[df.interface == 1]
    dft = pd.DataFrame()

    for col in cols:
        x = dfs0[col].to_numpy()
        y = dfs1[col].to_numpy()
        if dfs1[col].mean() - dfs0[col].mean() > 0:
            Rbool = "True"
        else:
            Rbool = "False"
        dft.loc[col, "higher for interface residues"] = Rbool

        p_boot = calc_ttest_pvalue_from_bootstrapped_data(x, y, equal_var=True, B=bootstrap_repetitions)
        print(p_boot)

        dft.loc[col, "p_boot"] = p_boot

        p_ttest = stats.ttest_ind(x, y)[1]
        dft.loc[col, "p_ttest"] = p_ttest

        # print(p)
        # if p < 0.05:
        # print(col)

    dft.sort_values("p_ttest", ascending=True, inplace=True)

    cutoff = 0.3
    dfcorr = df[cols].corr()
    for col in cols:
        correlated_features = dfcorr[col].loc [dfcorr[col] > cutoff].index.tolist()
        correlated_features.remove(col)
        if len(correlated_features) > 0:
            dft.loc[col, "correlated_features"] = ", ".join(correlated_features)
        else:
            dft.loc[col, "correlated_features"] = ""

    dft.index.name = "feature"

    dft_sign = dft[dft.p_ttest < 0.05].copy()
    dft_sign.drop("p_boot", axis=1, inplace=True)
    sign_coevolution_feats = [x for x in dft_sign.index.tolist() if "MI" in x or "DI" in x]
    sign_coevolution_feats = sorted(sign_coevolution_feats)
    list_sign_bootstrapped = ", ".join(sign_coevolution_feats)
    print("{} coevolution features significantly different between int and non-interface, REGULAR TTEST\n".format(len(sign_coevolution_feats)))
    print(list_sign_bootstrapped)

    dft_sign_boot = dft[dft.p_boot < 0.05].copy()
    dft_sign_boot.drop("p_ttest", axis=1, inplace=True)
    sign_coevolution_feats = [x for x in dft_sign_boot.index.tolist() if "MI" in x or "DI" in x]
    sign_coevolution_feats = sorted(sign_coevolution_feats)
    list_sign_bootstrapped = ", ".join(sign_coevolution_feats)
    print("{} coevolution features significantly different between int and non-interface, BOOTSTRAPPED TTEST\n".format(len(sign_coevolution_feats)))
    print(list_sign_bootstrapped)

    dft["p_boot"] = dft["p_boot"].replace(0.0, f"<{1/bootstrap_repetitions}")

    with pd.ExcelWriter(ttest_pvalues_bootstrapped_data_using_traindata_selected_features_xlsx) as writer:
        dft.to_excel(writer, sheet_name="all")
        dft_sign.to_excel(writer, sheet_name="sign_TTEST")
        dft_sign_boot.to_excel(writer, sheet_name="sign_boot_TTEST")

    logging.info('finished conduct_ttest_for_selected_features_used_in_model')