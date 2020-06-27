import numpy as np
from matplotlib import pyplot as plt

import thoipapy.utils


def create_var_imp_plot(df_imp, colour_dict, variable_importance_png, n_features_in_plot):
    """Plot function for fig_feat_import_from_mean_decrease_impurity, allowing a variable number of features.
    """
    df_sel = df_imp.iloc[:n_features_in_plot, :].copy()

    # regular trees, or Totally Randomized Trees
    model_types = ["", "_TRT"]

    for model_type in model_types:

        # add suffix for the totally randomised trees
        variable_importance_png = str(variable_importance_png)[:-4] + model_type + ".png"

        df_sel.sort_values("mean_decrease_impurity{}".format(model_type), ascending=True, inplace=True)
        #min_ = df_sel.mean_decrease_impurity.min()

        # determine the plot height by the number of features
        # currently set for 30
        plot_height = 4 * n_features_in_plot / 30
        figsize = np.array([4.42, plot_height])
        fig, ax = plt.subplots(figsize=figsize)

        TUMblue = colour_dict["TUM_colours"]['TUMBlue']
        df_sel["mean_decrease_impurity{}".format(model_type)].plot(kind="barh", color="#17a8a5", ax=ax)# xerr=df_sel["std"]
        ax.errorbar(df_sel["mean_decrease_impurity{}".format(model_type)], range(len(df_sel.index)), xerr=df_sel["std{}".format(model_type)], fmt="none", ecolor="k", ls="none", capthick=0.5, elinewidth=0.5, capsize=1, label=None)

        ax.set_xlim(0)

        ax.set_ylabel("")
        ax.set_xlabel("variable importance\n(mean decrease impurity)")
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(variable_importance_png, dpi=240)
        #fig.savefig(thoipapy.utils.pdf_subpath(variable_importance_png))


