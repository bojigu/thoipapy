from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from thoipapy.utils import normalise_0_1

import thoipapy.utils
from thoipapy.ML_model.train_model import return_classifier_with_loaded_ensemble_parameters
from thoipapy.feature_importance.mean_decrease_impurity import calculate_mean_decrease_impurity_for_dataset
from thoipapy.utils import create_colour_lists, normalise_between_2_values


def plot_feature_importance(s, logging):
    """Create figures showing ML feature importance.
    """
    plt.style.use('seaborn-whitegrid')
    plt.rcParams['errorbar.capsize'] = 1
    plt.rcParams.update({'font.size': 4})
    colour_dict = create_colour_lists()

    trainset = "set08"

    # input
    mean_decrease_impurity_all_features_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/01_feat_imp_MDI_before_feature_seln.csv"
    feat_imp_mean_decrease_accuracy_xlsx = Path(s["data_dir"]) / f"results/{s['setname']}/feat_imp/feat_imp_mean_decrease_accuracy.xlsx"
    tuned_ensemble_parameters_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/04_tuned_ensemble_parameters.csv"
    # output
    variable_importance_png = Path(s["data_dir"]) / "results" / s["setname"] / "feat_imp/FigS17_BZ13_feature_importance.png"
    variable_importance_xlsx = Path(s["data_dir"]) / "results" / s["setname"] / "feat_imp/FigS17_BZ13_feature_importance.xlsx"

    train_data_after_first_feature_seln_csv = Path(s["data_dir"]) / f"results/{s['setname']}/train_data/03_train_data_after_first_feature_seln.csv"
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
    df.sort_values("MDA_AUBOC", ascending=False, inplace=True)

    with pd.ExcelWriter(variable_importance_xlsx) as writer:
        df.to_excel(writer, sheet_name="var_import")

        df_norm = pd.DataFrame()
        for col in df.columns:
            df_norm[col] = normalise_0_1(df[col])[0]
        df_norm.index = df.index
        df_norm.to_excel(writer, sheet_name="var_import_norm")
        writer.save()

    # min_max_scaler = MinMaxScaler()
    # x_scaled = min_max_scaler.fit_transform(df.values)
    # df_norm = pd.DataFrame(x_scaled)
    # df_norm.columns = df.columns
    # df_norm.index = df.index
    # df_norm.sort_values("MDA_AUBOC", ascending=False, inplace=True)

    n_features_in_plot = df.shape[0]
    # determine the plot height by the number of features
    # currently set for 30
    plot_height = 4 * n_features_in_plot / 30
    figsize = np.array([4.42, plot_height])
    fig, ax = plt.subplots()

    TUMblue = colour_dict["TUM_colours"]['TUMBlue']

    df["MDA_AUBOC"].plot(kind="bar", ax=ax, color=TUMblue)
    # df["MDA_PR_AUC"].plot(kind="bar", ax=ax, position=1, color="#17a8a5")

    ax.set_xlabel("")
    ax.set_ylabel("variable importance\n(mean decrease AUBOC5)")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(variable_importance_png, dpi=240)
    fig.savefig(thoipapy.utils.pdf_subpath(variable_importance_png))
