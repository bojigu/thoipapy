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


def calc_ttest_pvalue_from_bootstrapped_data(x, y, equal_var=False, B=100000, plot=False):
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
    ttest_pvalues_bootstrapped_data_using_traindata_selected_features_xlsx = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/ttest/ttest_pvalues_bootstrapped_data_using_traindata_selected_features(train{trainsetname}).xlsx"

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

        p_bootstrapped_ttest = calc_ttest_pvalue_from_bootstrapped_data(x, y, equal_var=True)
        print(p_bootstrapped_ttest)

        dft.loc[col, "p_bootstrapped_ttest"] = p_bootstrapped_ttest

        p_ttest = stats.ttest_ind(x, y)[1]
        dft.loc[col, "p_ttest"] = p_ttest

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
    dft_sign.drop("p_bootstrapped_ttest", axis=1, inplace=True)
    sign_coevolution_feats = [x for x in dft_sign.index.tolist() if "MI" in x or "DI" in x]
    sign_coevolution_feats = sorted(sign_coevolution_feats)
    list_sign_bootstrapped = ", ".join(sign_coevolution_feats)
    print("{} coevolution features significantly different between int and non-interface, REGULAR TTEST\n".format(len(sign_coevolution_feats)))
    print(list_sign_bootstrapped)

    dft_sign_boot = dft[dft.p_bootstrapped_ttest < 0.05].copy()
    dft_sign_boot.drop("p_ttest", axis=1, inplace=True)
    sign_coevolution_feats = [x for x in dft_sign_boot.index.tolist() if "MI" in x or "DI" in x]
    sign_coevolution_feats = sorted(sign_coevolution_feats)
    list_sign_bootstrapped = ", ".join(sign_coevolution_feats)
    print("{} coevolution features significantly different between int and non-interface, BOOTSTRAPPED TTEST\n".format(len(sign_coevolution_feats)))
    print(list_sign_bootstrapped)

    dft["p_bootstrapped_ttest_formatted"] = dft["p_bootstrapped_ttest"].apply(convert_pvalue_to_text)

    with pd.ExcelWriter(ttest_pvalues_bootstrapped_data_using_traindata_selected_features_xlsx) as writer:
        dft.to_excel(writer, sheet_name="all")
        dft_sign.to_excel(writer, sheet_name="significant")
        dft_sign_boot.to_excel(writer, sheet_name="significant_bootstrapped")

    logging.info('finished conduct_ttest_for_selected_features_used_in_model')


def convert_pvalue_to_text(p, bootstrap_replicates=100000):

    n_significant_figures = len(str(bootstrap_replicates)) - 1

    formatter = "{:0.%if}" % n_significant_figures

    assert isinstance(p, float)
    if p == 0.0:
        lowest_possible_pvalue_based_on_n_bootstrap_replicates = 1 / bootstrap_replicates
        min_value_tested = f"{float(1 / bootstrap_replicates):}"
        return "<" + formatter.format(lowest_possible_pvalue_based_on_n_bootstrap_replicates)

    return formatter.format(p)
