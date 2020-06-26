import os
import pickle
from pathlib import Path
from typing import Union

import joblib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy import interp
from sklearn.metrics import roc_curve, auc

import thoipapy.common
import thoipapy.figs
import thoipapy.utils
import thoipapy.validation
from thoipapy.figs.create_BOcurve_files import parse_BO_data_csv_to_excel, save_BO_linegraph_and_barchart, save_extra_BO_figs


def run_testset_trainset_validation(s, logging):
    if not os.path.exists(os.path.join(s["thoipapy_data_folder"], "Results", "Bo_curve")):
        os.makedirs(os.path.join(s["thoipapy_data_folder"], "Results", "Bo_curve"))

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset
    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    validate_LIPS_for_testset(s, logging)
    validate_LIPS_for_testset(s, logging, LIPS_name="LIPS_surface_ranked", pred_col="LIPS_surface_ranked")

    validate_THOIPA_for_testset_trainset_combination(s, test_set_list, train_set_list, logging)


def validate_THOIPA_for_testset_trainset_combination(s, test_set_list, train_set_list, logging):
    """ Creates ROC and BO-curve for a particular testset-trainset combination.

    Parameters
    ----------
	s : dict
        Settings dictionary for figures.
    test_set_list : list
        List of test datasets in selection
        E.g. ["set03", "set31"]
    train_set_list : list
        List of training datasets in selection
        E.g. ["set02", "set04"]

    Saved Files
    -----------
    THOIPA_pred_csv : csv
        THOIPA result for this testset-trainset combination
        Columns = "residue_num", "residue_name", "THOIPA"
        Index = range index of residues
    combined_incl_THOIPA_csv : csv
        The combined file with all features. THOIPA prediction is added as a new column
    THOIPA_ROC_pkl : pickle
        Pickled output dictionary with ROC curves
        keys = accessions
        values = dictionary with fpr, tpr etc for each protein
        Could not be saved easily as a dataframe, because the number of residues is different for each protein

    """
    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    for n, train_set in enumerate(train_set_list):
        trainsetname = "set{:02d}".format(int(train_set))
        model_pkl = os.path.join(s["Results_folder"], trainsetname, "{}_ML_model.lpkl".format(trainsetname))

        for test_set in test_set_list:
            testsetname = "set{:02d}".format(int(test_set))

            #BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname))
            BO_curve_folder = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/THOIPA.train{trainsetname}"
            #THOIPA_BO_curve_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "data", "Test{}_Train{}.THOIPA.best_overlap_data.csv".format(testsetname, trainsetname))
            THOIPA_BO_curve_data_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/THOIPA.train{trainsetname}/THOIPA.best_overlap_data.csv"
            #THOIPA_ROC_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}_Train{}.THOIPA".format(testsetname, trainsetname), "data", "Test{}_Train{}.THOIPA.ROC_data.pkl".format(testsetname, trainsetname))
            THOIPA_ROC_pkl = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/THOIPA.train{trainsetname}/ROC_data.pkl"


            BO_data_excel = os.path.join(BO_curve_folder, "data", "BO_curve_data.xlsx")
            BO_linechart_png = os.path.join(BO_curve_folder, "BO_linechart.png")
            BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")

            thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

            testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["sets_folder"])

            testdataset_df = pd.read_excel(testset_path)
            THOIPA_BO_data_df = pd.DataFrame()
            #LIPS_BO_data_df = pd.DataFrame()

            # save all outputs to a cross-validation dictionary, to be saved as a pickle file
            xv_dict_THOIPA = {}
            mean_tpr = 0.0
            mean_fpr = np.linspace(0, 1, 100)

            for i in testdataset_df.index:
                acc = testdataset_df.loc[i, "acc"]
                database = testdataset_df.loc[i, "database"]
                acc_db = acc + "-" + database
                testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                      "{}.surr20.gaps5.combined_features.csv".format(acc))
                #THOIPA_pred_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database,
                #                                      "{}.THOIPA.train{}.csv".format(acc, trainsetname))
                THOIPA_pred_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/data/{database}/{acc}.THOIPA.train{trainsetname}.csv"
                #combined_incl_THOIPA_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database,
                #                                      "{}.THOIPA_incl_combined.train{}.csv".format(acc, trainsetname))
                combined_incl_THOIPA_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/data/{database}/{acc}.THOIPA_incl_combined.train{trainsetname}.csv"


                thoipapy.utils.make_sure_path_exists(combined_incl_THOIPA_csv, isfile=True)

                combined_incl_THOIPA_df = save_THOIPA_pred_indiv_prot(s, model_pkl, testdata_combined_file, THOIPA_pred_csv, combined_incl_THOIPA_csv, trainsetname, logging)

                #######################################################################################################
                #                                                                                                     #
                #                           Processing BO curve data for each single protein                          #
                #                                                                                                     #
                #######################################################################################################

                combined_incl_THOIPA_df["LIPS_L*E"] = -1 * combined_incl_THOIPA_df["LIPS_L*E"]

                if database == "crystal" or database == "NMR":
                    # (it is closest distance and low value means high propencity of interfacial)
                    combined_incl_THOIPA_df["interface_score"] = -1 * combined_incl_THOIPA_df["interface_score"]

                THOIPA_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, combined_incl_THOIPA_df, pred_col=f"THOIPA_train{trainsetname}")

                if THOIPA_BO_data_df.empty:
                    THOIPA_BO_data_df = THOIPA_BO_single_prot_df
                else:
                    THOIPA_BO_data_df = pd.concat([THOIPA_BO_data_df, THOIPA_BO_single_prot_df], axis=1, join="outer")

                #######################################################################################################
                #                                                                                                     #
                #                     Processing ROC data for each single protein, saving to nested dict              #
                #                                                                                                     #
                #######################################################################################################

                df_for_roc = combined_incl_THOIPA_df.dropna(subset=["interface_score"])

                predictor_name = f"THOIPA_train{trainsetname}"

                fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[predictor_name], drop_intermediate=False)
                roc_auc = auc(fpr, tpr)
                mean_tpr += interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0

                xv_dict_THOIPA[acc_db] = {"fpr" : fpr, "tpr" : tpr, "roc_auc" : roc_auc}

            #######################################################################################################
            #                                                                                                     #
            #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
            #                                                                                                     #
            #######################################################################################################

            THOIPA_BO_data_df.to_csv(THOIPA_BO_curve_data_csv)

            #THOIPA_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(THOIPA_BO_curve_data_csv, BO_curve_folder, names_excel_path)
            parse_BO_data_csv_to_excel(THOIPA_BO_curve_data_csv, BO_data_excel, logging)
            AUC_ser = pd.Series(xv_dict_THOIPA[acc_db]["roc_auc"])
            AUBOC10 = save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser)

            if "you_want_more_details" == "TRUE":
                other_figs_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/other_figs"
                save_extra_BO_figs(BO_data_excel, other_figs_path)

            #######################################################################################################
            #                                                                                                     #
            #                     Processing dictionary with ROC data, saving to pickle                           #
            #                                                                                                     #
            #######################################################################################################
            mean_tpr /= testdataset_df.shape[0]
            mean_tpr[-1] = 1.0

            mean_roc_auc = auc(mean_fpr, mean_tpr)

            ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_THOIPA}
            ROC_out_dict["true_positive_rate_mean"] = mean_tpr
            ROC_out_dict["false_positive_rate_mean"] = mean_fpr
            ROC_out_dict["mean_roc_auc"] = mean_roc_auc

            # save dict as pickle
            with open(THOIPA_ROC_pkl, "wb") as f:
                pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

            create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl)

            logging.info("test{}_train{} AUC({:.03f}), AUBOC10({:.2f}). ({})".format(testsetname, trainsetname, mean_roc_auc, AUBOC10, BO_barchart_png))


def validate_LIPS_for_testset(s, logging, LIPS_name="LIPS_LE", pred_col="LIPS_L*E"):

    names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")
    namedict = thoipapy.utils.create_namedict(names_excel_path)

    # create list of test and train datasets
    # if only one is given, make a list with only one dataset
    test_set_list, train_set_list = thoipapy.figs.fig_utils.get_test_and_train_set_lists(s)

    for test_set in test_set_list:
        testsetname = "set{:02d}".format(int(test_set))
        LIPS_BO_curve_data_csv = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/{LIPS_name}/{LIPS_name}.best_overlap_data.csv.csv"
        #LIPS_BO_curve_data_csv = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name), "Test{}.{}.best_overlap_data.csv".format(testsetname, LIPS_name))
        #BO_curve_folder = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name))
        BO_curve_folder = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/{LIPS_name}"
        LIPS_ROC_pkl = Path(s["thoipapy_data_folder"]) / "Results" / testsetname / f"blindvalidation/{LIPS_name}/ROC_data.pkl"
        #LIPS_ROC_pkl = os.path.join(s["thoipapy_data_folder"], "Results", "compare_testset_trainset", "data", "Test{}.{}".format(testsetname, LIPS_name), "data", "Test{}.{}.ROC_data.pkl".format(testsetname, LIPS_name))
        thoipapy.utils.make_sure_path_exists(LIPS_BO_curve_data_csv, isfile=True)
        BO_data_excel = os.path.join(BO_curve_folder, "data", "BO_curve_data.xlsx")
        BO_linechart_png = os.path.join(BO_curve_folder, "BO_linechart.png")
        BO_barchart_png = os.path.join(BO_curve_folder, "AUBOC10_barchart.png")
        thoipapy.utils.make_sure_path_exists(BO_data_excel, isfile=True)

        testset_path = thoipapy.common.get_path_of_protein_set(testsetname, s["sets_folder"])

        testdataset_df = pd.read_excel(testset_path)
        LIPS_BO_data_df = pd.DataFrame()

        # save all outputs to a cross-validation dictionary, to be saved as a pickle file
        xv_dict_LIPS = {}
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        for i in testdataset_df.index:
            acc = testdataset_df.loc[i, "acc"]
            database = testdataset_df.loc[i, "database"]
            acc_db = acc + "-" + database

            testdata_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database,
                                                  "{}.surr20.gaps5.combined_features.csv".format(acc))

            combined_df = pd.read_csv(testdata_combined_file, index_col=0)

            #######################################################################################################
            #                                                                                                     #
            #                           Processing BO curve data for each single protein                          #
            #                                                                                                     #
            #######################################################################################################
            # SAVE LIPS PREDICTION DATA
            # this is somewhat inefficient, as it is conducted for every test dataset
            LIPS_pred_csv = os.path.join(os.path.dirname(s["thoipapy_data_folder"]), "Features", "Predictions", "testset_trainset", database, "{}.LIPS_pred.csv".format(acc, testsetname))
            LIPS_pred_df = combined_df[["residue_name", "residue_num", "LIPS_polarity", "LIPS_entropy", "LIPS_L*E", "LIPS_surface", "LIPS_surface_ranked"]]
            thoipapy.utils.make_sure_path_exists(LIPS_pred_csv, isfile=True)
            LIPS_pred_df.to_csv(LIPS_pred_csv)

            if pred_col == "LIPS_L*E":
                combined_df[pred_col] = -1 * combined_df[pred_col]

            if database == "crystal" or database == "NMR":
                # (it is closest distance and low value means high propencity of interfacial)
                combined_df["interface_score"] = -1 * combined_df["interface_score"]

            LIPS_BO_single_prot_df = thoipapy.figs.fig_utils.calc_best_overlap(acc_db, combined_df, experiment_col="interface_score", pred_col=pred_col)

            if LIPS_BO_single_prot_df.empty:
                LIPS_BO_data_df = LIPS_BO_single_prot_df
            else:
                LIPS_BO_data_df = pd.concat([LIPS_BO_data_df, LIPS_BO_single_prot_df], axis=1, join="outer")

            #######################################################################################################
            #                                                                                                     #
            #                     Processing ROC data for each single protein, saving to nested dict              #
            #                                                                                                     #
            #######################################################################################################

            df_for_roc = combined_df.dropna(subset=["interface_score"])

            fpr, tpr, thresholds = roc_curve(df_for_roc.interface, df_for_roc[pred_col], drop_intermediate=False)
            roc_auc = auc(fpr, tpr)
            mean_tpr += interp(mean_fpr, fpr, tpr)
            mean_tpr[0] = 0.0

            xv_dict_LIPS[acc_db] = {"fpr" : fpr, "tpr" : tpr, "roc_auc" : roc_auc}

        #######################################################################################################
        #                                                                                                     #
        #      Processing BO CURVE data, saving to csv and running the BO curve analysis script               #
        #                                                                                                     #
        #######################################################################################################

        LIPS_BO_data_df.to_csv(LIPS_BO_curve_data_csv)
        names_excel_path = os.path.join(s["dropbox_dir"], "protein_names.xlsx")

        #LIPS_linechart_mean_obs_and_rand = analyse_bo_curve_underlying_data(LIPS_BO_curve_data_csv, BO_curve_folder, names_excel_path)

        #parse_BO_data_csv_to_excel(LIPS_BO_curve_data_csv, BO_curve_folder, names_excel_path)
        parse_BO_data_csv_to_excel(LIPS_BO_curve_data_csv, BO_data_excel, logging)
        AUC_ser = pd.Series(xv_dict_LIPS[acc_db]["roc_auc"])
        AUBOC10 = save_BO_linegraph_and_barchart(s, BO_data_excel, BO_linechart_png, BO_barchart_png, namedict, logging, AUC_ser)

        if "you_want_more_details" == "TRUE":
            other_figs_path: Union[Path, str] = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/crossvalidation/other_figs"
            save_extra_BO_figs(BO_data_excel, other_figs_path)

        #######################################################################################################
        #                                                                                                     #
        #                     Processing dictionary with ROC data, saving to pickle                           #
        #                                                                                                     #
        #######################################################################################################
        mean_tpr /= testdataset_df.shape[0]
        mean_tpr[-1] = 1.0

        mean_roc_auc = auc(mean_fpr, mean_tpr)

        ROC_out_dict = {"xv_dict_THOIPA" : xv_dict_LIPS}
        ROC_out_dict["true_positive_rate_mean"] = mean_tpr
        ROC_out_dict["false_positive_rate_mean"] = mean_fpr
        ROC_out_dict["mean_roc_auc"] = mean_roc_auc

        # save dict as pickle
        with open(LIPS_ROC_pkl, "wb") as f:
            pickle.dump(ROC_out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

        create_ROC_fig_for_testset_trainset_combination(LIPS_ROC_pkl)

        logging.info("test{}.{} AUC({:.03f}), AUBOC10({:.2f}). ({})".format(testsetname, LIPS_name, mean_roc_auc, AUBOC10, BO_barchart_png))


def create_ROC_fig_for_testset_trainset_combination(THOIPA_ROC_pkl):

    #plt.rcParams.update({'font.size': 6})
    ROC_pkl_basename = os.path.basename(THOIPA_ROC_pkl)[:-4]
    ROC_pkl_dir = os.path.dirname(THOIPA_ROC_pkl)

    ROC_png = os.path.join(ROC_pkl_dir, "{}.ROC.png".format(ROC_pkl_basename))
    thoipapy.utils.make_sure_path_exists(ROC_png, isfile=True)

    # open pickle file
    with open(THOIPA_ROC_pkl, "rb") as f:
        ROC_out_dict = pickle.load(f)

    xv_dict_THOIPA = ROC_out_dict["xv_dict_THOIPA"]

    figsize = np.array([3.42, 3.42]) * 2 # DOUBLE the real size, due to problems on Bo computer with fontsizes
    fig, ax = plt.subplots(figsize=figsize)

    for acc_db in xv_dict_THOIPA:
        roc_auc = xv_dict_THOIPA[acc_db]["roc_auc"]
        ax.plot(xv_dict_THOIPA[acc_db]["fpr"], xv_dict_THOIPA[acc_db]["tpr"], lw=1, label='{} ({:0.2f})'.format(acc_db, roc_auc), alpha=0.8)

    #mean_roc_auc = auc(df_xv["false_positive_rate"], df_xv["true_positive_rate"])
    mean_roc_auc = ROC_out_dict["mean_roc_auc"]

    ax.plot(ROC_out_dict["false_positive_rate_mean"], ROC_out_dict["true_positive_rate_mean"], color="k", label='mean (area = %0.2f)' % mean_roc_auc, lw=1.5)
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='random')
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(ROC_png, dpi=240)
    fig.savefig(thoipapy.utils.pdf_subpath(ROC_png))
    #fig.savefig(ROC_png[:-4] + ".pdf")


def save_THOIPA_pred_indiv_prot(s, model_pkl, testdata_combined_file, THOIPA_pred_csv, test_combined_incl_pred, trainsetname, logging):

    combined_incl_THOIPA_df = pd.read_csv(testdata_combined_file, sep=',', engine='python', index_col=0)
    train_data_filtered = Path(s["thoipapy_data_folder"]) / f"Results/{s['setname']}/train_data/train_data_filtered.csv"
    #combined_incl_THOIPA_df = combined_incl_THOIPA_df.dropna()

    #drop_cols_not_used_in_ML
    #X=train_df.drop(train_features_del,axis=1)
    # X = thoipapy.features.RF_Train_Test.drop_cols_not_used_in_ML(logging, train_df, s["excel_file_with_settings"])
    # y = train_df["interface"]
    # clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    # fit = clf.fit(X,y)

    fit = joblib.load(model_pkl)

    # if database == "ETRA":
    #     dirupt_path = os.path.join(s["toxr_disruption_folder"], "{}_mul_scan_average_data.xlsx".format(acc))
    #     ddf = pd.read_excel(dirupt_path, index_col=0)
    #     interface_score = ddf.Disruption
    #     interface_score = interface_score      #(it is closest experimental disruption and high value means high propencity of interfacial)

    #Lips_score = test_df.LIPS_polarity * test_df.LIPS_entropy

    #tX=test_df.drop(test_features_del,axis=1)
    test_X = thoipapy.validation.feature_selection.drop_cols_not_used_in_ML(logging, combined_incl_THOIPA_df, s["excel_file_with_settings"])

    try:
        prob_arr = fit.predict_proba(test_X)[:, 1]
    except ValueError:
        df_train_data_used_for_model_csv = pd.read_csv(train_data_filtered)
        del df_train_data_used_for_model_csv[s["bind_column"]]
        train_cols = set(df_train_data_used_for_model_csv.columns)
        test_cols = set(test_X.columns)
        cols_not_in_train = test_cols - train_cols
        cols_not_in_test = train_cols - test_cols
        logging.warning("cols_not_in_train : {}\ncols_not_in_test : {}".format(cols_not_in_train, cols_not_in_test))
        prob_arr = fit.predict_proba(test_X)[:, 1]

    # if hasattr(clf,'predict_proba'):
    #     prob_arr = fit.predict_proba(tX)[:,1]
    # else:
    #     prob_arr = fit.decision_function(tX)
    #     prob_arr = (prob_arr - prob_arr.min())/(prob_arr.max() - prob_arr.min())

    combined_incl_THOIPA_df[f"THOIPA_train{trainsetname}"] = prob_arr
    # save full combined file with THOIPA prediction
    combined_incl_THOIPA_df.to_csv(test_combined_incl_pred)
    # save only THOIPA prediction
    THOIPA_pred_df = combined_incl_THOIPA_df[["residue_num", "residue_name", f"THOIPA_train{trainsetname}"]]
    THOIPA_pred_df.to_csv(THOIPA_pred_csv)
    return combined_incl_THOIPA_df