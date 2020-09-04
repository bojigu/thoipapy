#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         BO ZENG
Created:        Monday November 20 12:33:08 2017
Operation system required: Linux (currently not available for windows)
Dependencies:   Python 3.5
                numpy
                Bio
                freecontact (currently only availble in linux)
                pandas
Purpose:        Self-interacting single-pass membrane protein interface residues prediction

"""
# I'm getting sick of the warnings that occur due to imported seaborn and statsmodels.stats.api modules, and have nothing to do with your code.
# you should turn on warnings once a month to check if there is anything related to your code
import warnings

import thoipapy.experimental_data.add_experimental_data_to_train_set
import thoipapy.feature_importance.mean_decrease_accuracy
import thoipapy.feature_importance.mean_decrease_impurity
import thoipapy.features.physical_parameters
import thoipapy.validation.gather
from thoipapy.clustering.pairwise_aln_similarity_matrix import create_identity_matrix_from_protein_set
from thoipapy.utils import get_testsetname_trainsetname_from_run_settings

warnings.filterwarnings("ignore")

import argparse
import os
import platform
import sys
import pandas as pd
import thoipapy

parser = argparse.ArgumentParser()

parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "\Path\to\your\settingsfile.xlsx"')

if __name__ == "__main__":
    sys.stdout.write('\nRun thoipapy as follows:')
    sys.stdout.write(r'python \Path\to\run.py -s \Path\to\your\settingsfile.xlsx')
    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    s = thoipapy.common.create_settingdict(args.s)

    ##############################################################################################
    #                                                                                            #
    #                               setname, logging, results folder                             #
    #                                                                                            #
    ##############################################################################################
    sets_folder = os.path.join(s["dropbox_dir"], "sets")

    # if multiple sets need to be run, split them by comma
    if isinstance(s["set_number"], str) and "," in s["set_number"]:
        list_protein_sets = [int(n) for n in s["set_number"].split(",")]
    else:
        list_protein_sets = [s["set_number"]]

    for set_number in list_protein_sets:
        s["set_number"] = set_number
        # define set name, which should be in the excel file name
        setname = "set{:02d}".format(s["set_number"])
        # add to the dictionary itself
        s["setname"] = setname
        # create a results folder for that set
        if not os.path.isdir(os.path.join(s["thoipapy_data_folder"], "Results", setname)):
            os.makedirs(os.path.join(s["thoipapy_data_folder"], "Results", setname))

        logging = thoipapy.common.setup_keyboard_interrupt_and_error_logging(s, setname)
        logging.info("STARTING PROCESSING OF {}.".format(setname))

        set_path = thoipapy.common.get_path_of_protein_set(setname, sets_folder)

        ##############################################################################################
        #                                                                                            #
        #                     open and process a set of protein sequences                            #
        #                                                                                            #
        ##############################################################################################
        # load the protein set (e.g. set01.xlsx) as a dataframe
        df_set = pd.read_excel(set_path, sheet_name='proteins')

        # create list of uniprot accessions to run
        acc_list = df_set.acc.tolist()
        sys.stdout.write("settings file : {}\nsettings : {}\nprotein set number {}, acc_list : {}\n".format(os.path.basename(args.s), s, s["set_number"], acc_list))
        sys.stdout.flush()

        dfset = thoipapy.common.process_set_protein_seqs(s, setname, df_set, set_path)

        # create a database label. Either crystal, NMR, ETRA or "mixed"
        subsets = df_set["database"].unique()
        if len(subsets.shape) == 1:
            database_for_full_set = subsets[0]
        else:
            database_for_full_set = "mixed"

        if "create_identity_matrix_from_set_seqs" in s:
            if s["create_identity_matrix_from_set_seqs"]:
                create_identity_matrix_from_protein_set(s, logging)


        ###################################################################################################
        #                                                                                                 #
        #                  calculate closedistance from NMR and crystal structures                        #
        #                                                                                                 #
        ###################################################################################################

        # DEPRECATED. Use atom_dist module instead to get closest heavy-atom distances
        # if s["Get_Tmd_Homodimers"] :
        #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.download_xml_get_alphahelix_get_homo_pair(s, logging)
        #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.Download_trpdb_Calc_inter_rr_pairs(s, logging)
        #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.create_redundant_interact_homodimer_rm_shorttm(s, logging)
        #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.extract_crystal_resolv035_interact_pairs_and_create_fasta_file(s, logging)
        #     thoipapy.structures.deprecated.get_tmd_nr_homodimer.create_multiple_bind_closedist_file(s, logging)
        #     pass

        if s["retrospective_coevolution"]:
            #thoipapy.figs.retrospective.calc_retrospective_coev_from_list_interf_res(s, dfset, logging)
            thoipapy.figs.retrospective.calc_retrospective_coev_from_struct_contacts(s, dfset, logging)

        #if s["calc_NMR_closedist"] :
        #    thoipapy.structures.deprecated.NMR_data.calc_closedist_from_NMR_best_model(s)

        if s["Atom_Close_Dist"]:
            infor = thoipapy.experimental_data.closest_heavy_atom_dist.homodimer_residue_closedist_calculate_from_complex(thoipapy, s, logging)
            sys.stdout.write(infor)

        ###################################################################################################
        #                                                                                                 #
        #                   homologues download from NCBI. parse, filter and save                         #
        #                                                                                                 #
        ###################################################################################################

        if s["run_retrieve_NCBI_homologues_with_blastp"]:
            thoipapy.homologues.NCBI_download.download_homologues_from_ncbi_mult_prot(s, df_set, logging)

        if s["run_parse_homologues_xml_into_csv"]:
            thoipapy.homologues.NCBI_parser.parse_NCBI_xml_to_csv_mult_prot(s, df_set, logging)

        if s["parse_csv_homologues_to_alignment"]:
            thoipapy.homologues.NCBI_parser.extract_filtered_csv_homologues_to_alignments_mult_prot(s, df_set, logging)


        ###################################################################################################
        #                                                                                                 #
        #                   machine learning feature calculation                                             #
        #                                                                                                 #
        ###################################################################################################

        if s["pssm_calculation"]:
            thoipapy.features.pssm.create_PSSM_from_MSA_mult_prot(s, df_set, logging)

        if s["entropy_calculation"]:
            thoipapy.features.entropy.entropy_calculation_mult_prot(s, df_set, logging)

        if s["rate4site_calculation"]:
            thoipapy.features.rate4site.rate4site_calculation(s, df_set, logging)

        if s["coevolution_calculation"]:
            if "Windows" in platform.system():
                sys.stdout.write("\n Freecontact cannot be run in Windows! Skipping coevolution_calculation_with_freecontact_mult_prot.")
                thoipapy.features.freecontact.parse_freecontact_coevolution_mult_prot(s, df_set, logging)
            else:
                thoipapy.features.freecontact.coevolution_calculation_with_freecontact_mult_prot(s, df_set, logging)
                thoipapy.features.freecontact.parse_freecontact_coevolution_mult_prot(s, df_set, logging)

        if s["clac_relative_position"]:
            thoipapy.features.relative_position.calc_relative_position_mult_prot(s, df_set, logging)

        if s["calc_lipo_from_pssm"]:
            thoipapy.features.lipophilicity.lipo_from_pssm_mult_prot(s, df_set, logging)

        if s["lips_score_calculation"]:
            thoipapy.features.lips.LIPS_score_calculation_mult_prot(s, df_set, logging)
            thoipapy.features.lips.parse_LIPS_score_mult_prot(s, df_set, logging)

        if s["motifs_from_seq"]:
            thoipapy.features.motifs.motifs_from_seq_mult_protein(s, df_set, logging)

        if s["combine_feature_into_train_data"]:
            thoipapy.features.combine_features.combine_all_features_mult_prot(s, df_set, logging)
            thoipapy.features.physical_parameters.add_physical_parameters_to_features_mult_prot(s, df_set, logging)
            thoipapy.experimental_data.add_experimental_data_to_train_set.add_experimental_data_to_combined_features_mult_prot(s, df_set, logging)
            if s["generate_randomised_interfaces"]:
                thoipapy.validation.random_interface.add_random_interface_to_combined_features_mult_prot(s, df_set, logging)
            if "add_PREDDIMER_TMDOCK_to_combined_features" in s:
                if s["add_PREDDIMER_TMDOCK_to_combined_features"]:
                    thoipapy.features.preddimer_tmdock.add_PREDDIMER_TMDOCK_to_combined_features_mult_prot(s, df_set, logging)
            if s["remove_crystal_hetero"]:
                thoipapy.experimental_data.remove_hetero_contacts.remove_crystal_hetero_contact_residues_mult_prot(s, df_set, logging)
            thoipapy.features.combine_features.combine_all_train_data_for_machine_learning(s, df_set, logging)

        ###################################################################################################
        #                                                                                                 #
        #                                    model validation                                             #
        #                                                                                                 #
        ###################################################################################################
        if s["run_feature_selection"]:
            thoipapy.feature_importance.mean_decrease_impurity.get_initial_ensemble_parameters_before_feature_selection(s, logging)
            thoipapy.feature_importance.mean_decrease_impurity.calc_feat_import_using_MDI_before_feature_seln(s, logging)
            thoipapy.feature_importance.mean_decrease_impurity.fig_feat_import_from_mean_decrease_impurity(s, logging)
            thoipapy.feature_importance.remove_duplicates.remove_duplicate_features_with_lower_MDI(s, logging)
            thoipapy.feature_importance.anova.select_best_features_with_anova(s, logging)
            thoipapy.feature_importance.ensemble_rfe.select_best_features_with_ensemble_rfe(s, logging)
            thoipapy.feature_importance.merge.merge_top_features_anova_ensemble(s, logging)

        if s["tune_ensemble_parameters"]:
            thoipapy.ML_model.tune.tune_ensemble_parameters_after_feature_seln(s, logging)

        if s["calc_feature_importances"]:
            thoipapy.feature_importance.mean_decrease_accuracy.calc_feat_import_from_mean_decrease_accuracy(s, logging)
            thoipapy.feature_importance.mean_decrease_accuracy.fig_feat_import_from_mean_decrease_accuracy(s, logging)

        if s["conduct_ttest"]:
            thoipapy.experimental_data.ttest_features.conduct_ttest_for_selected_features_used_in_model(s, logging)

        if s["train_machine_learning_model"]:
            thoipapy.ML_model.train_model.train_machine_learning_model(s, logging)

        if s["run_testset_trainset_validation"] == True:
            thoipapy.validation.testset_trainset.run_testset_trainset_validation(s, logging)

        ###################################################################################################
        #                                                                                                 #
        #                                               figures                                           #
        #                                                                                                 #
        ###################################################################################################

        Fontsize = s["Fontsize"]
        Filter = s["Filter"]
        Width= s["Width"]
        Size= s["Size"]
        Linewidth= s["Linewidth"]

        if s["compare_selected_predictors"] == True:
            thoipapy.figs.combine_BOcurve_files.compare_selected_predictors(s, logging)

        if s["calc_PREDDIMER_TMDOCK_closedist"] == True:
            thoipapy.figs.calc_PREDDIMER_TMDOCK_closedist.calc_closedist_from_PREDDIMER_TMDOCK_best_model(s, df_set, logging)

        if s["run_validation"] == True:
            sys.stdout.write("\n--------------- starting run_validation ---------------\n")
            namedict = thoipapy.utils.create_namedict(os.path.join(s["dropbox_dir"], "protein_names.xlsx"))
            THOIPA_predictor_name = "THOIPA_{}_LOO".format(s["set_number"])
            predictors = [THOIPA_predictor_name, "PREDDIMER", "TMDOCK", "LIPS_surface_ranked", "random"]
            testsetname, trainsetname = get_testsetname_trainsetname_from_run_settings(s)
            thoipa_trainsetname = f"thoipa.train{trainsetname}"
            if s["setname"] == testsetname:
                predictors.append(f"thoipa.train{trainsetname}")

            thoipapy.validation.tenfold.run_10fold_cross_validation(s, logging)
            thoipapy.validation.tenfold.create_10fold_cross_validation_fig(s, logging)

            thoipapy.validation.leave_one_out.run_LOO_validation(s, df_set, logging)
            thoipapy.validation.leave_one_out.create_LOO_validation_fig(s, df_set, logging)

            thoipapy.validation.combine_mult_predictors.merge_predictions(s, df_set, logging)

            thoipapy.validation.indiv_validation.collect_indiv_validation_data(s, df_set, logging, namedict, predictors, THOIPA_predictor_name, subsets)
            thoipapy.validation.indiv_validation.create_indiv_validation_figs(s, logging, namedict, predictors, THOIPA_predictor_name, subsets)

            thoipapy.validation.multiple_predictors.validate_multiple_predictors_and_subsets_auboc10(s, df_set, logging)
            thoipapy.validation.multiple_predictors.validate_multiple_predictors_and_subsets_auc(s, df_set, logging)

            thoipapy.validation.roc.create_ROC_all_residues(s, df_set, logging)
            thoipapy.validation.precision_recall.create_precision_recall_all_residues(s, df_set, logging)

            thoipapy.validation.gather.gather_validation_data_for_figs(s, df_set, logging)
            sys.stdout.write("\n--------------- finished run_validation ---------------\n")


        if s["create_merged_heatmap"] == True:
            thoipapy.figs.create_heatmap_from_merge_file.create_merged_heatmap(s, df_set, logging)

        if "download_10_homologues_from_ncbi" in s:
            if s["download_10_homologues_from_ncbi"] == True:
                thoipapy.homologues.NCBI_download.download_10_homologues_from_ncbi(s, df_set, logging)

        if s["plot_coev_vs_res_dist"] == True:
            thoipapy.figs.retrospective.calc_coev_vs_res_dist(s, dfset, logging)
            thoipapy.figs.retrospective.plot_coev_vs_res_dist(s, logging)

        thoipapy.setting.deployment_helper.docker_deployment_had_better_work_now()
        thoipapy.ML_model.deployment_helper2.docker_deployment_had_better_work_now2()

        # close the logger. A new one will be made for the next protein list.
        logging.info("FINISHED PROCESSING OF {}.".format(setname))
        logging.shutdown()
