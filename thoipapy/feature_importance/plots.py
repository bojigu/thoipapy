from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import thoipapy.utils
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters
from thoipapy.feature_importance.mean_decrease_impurity import calculate_mean_decrease_impurity_for_dataset
from thoipapy.utils import create_colour_lists


def plot_feature_importance(s, logging):
    """Create figures showing ML feature importance.
    """
    plt.style.use('seaborn-whitegrid')
    plt.rcParams['errorbar.capsize'] = 1
    plt.rcParams.update({'font.size': 4})
    colour_dict = create_colour_lists()

    full_dataset = "set05"
    trainset = "set08"

    # input
    mean_decrease_impurity_all_features_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/01_feat_imp_MDI_before_feature_seln.csv"
    feat_imp_mean_decrease_accuracy_xlsx = Path(s["thoipapy_data_folder"]) / f"results/{trainset}/feat_imp/feat_imp_mean_decrease_accuracy.xlsx"
    tuned_ensemble_parameters_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"
    # output
    variable_importance_png = Path(s["thoipapy_data_folder"]) / "results" / s["setname"] / "feat_imp/var_import.png"

    train_data_after_first_feature_seln_csv = Path(s["thoipapy_data_folder"]) / f"results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
    df_data = pd.read_csv(train_data_after_first_feature_seln_csv, index_col=0)
    y = df_data[s["bind_column"]]
    X = df_data[[c for c in df_data.columns if c != s["bind_column"]]]
    assert "interface" not in X.columns

    forest = return_classifier_with_loaded_ensemble_parameters(s, tuned_ensemble_parameters_csv)
    df_MDI = calculate_mean_decrease_impurity_for_dataset(X, y, forest, "", logging)
    df_MDI.set_index("feature", inplace=True)

    df_MDA = pd.read_excel(feat_imp_mean_decrease_accuracy_xlsx, sheet_name="single_feat", index_col=0)

    df = pd.DataFrame()
    df["mean_decrease_impurity"] = df_MDI["mean_decrease_impurity"]
    df["MDA_AUBOC"] = df_MDA["AUBOC"]
    df["MDA_PR_AUC"] = df_MDA["PR_AUC"]
    df.sort_values("MDA_AUBOC", ascending=True, inplace=True)

    n_features_in_plot = df.shape[0]
    # determine the plot height by the number of features
    # currently set for 30
    plot_height = 4 * n_features_in_plot / 30
    figsize = np.array([4.42, plot_height])
    fig, ax = plt.subplots(figsize=figsize)
    ax2 = ax.twiny()

    TUMblue = colour_dict["TUM_colours"]['TUMBlue']
    colour_MDI = TUMblue
    colour_MDA_AUBOC = "#17a8a5"
    df["mean_decrease_impurity"].plot(kind="barh", color=colour_MDI, ax=ax)  # xerr=df_sel["std"]
    #ax.errorbar(df_sel["mean_decrease_impurity"], range(len(df_sel.index)), xerr=df_sel["std"], fmt="none", ecolor="k", ls="none", capthick=0.5, elinewidth=0.5, capsize=1, label=None)
    df["MDA_AUBOC"].plot(kind="barh", color=colour_MDA_AUBOC, ax=ax2)  # xerr=df_sel["std"]

    #ax.set_xlim(0)

    ax.set_ylabel("")
    ax.set_xlabel("variable importance\n(mean decrease impurity)")
    #ax.grid(False)
    fig.tight_layout()
    fig.savefig(variable_importance_png, dpi=240)
    # fig.savefig(thoipapy.utils.pdf_subpath(variable_importance_png))
