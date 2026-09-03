#!/usr/bin/env python
"""

Author:   BO ZENG
Created:  Monday November 20 12:33:08 2017
Purpose:  Self-interacting single-pass membrane protein interface residues prediction

"""

# Avoid warnings due to imported seaborn and statsmodels.stats.api packages, etc. This has to run
# before those packages are imported, so it sits above every other import in the module.
import warnings

warnings.filterwarnings("ignore")

import os
import platform
import sys
from pathlib import Path

import pandas as pd

import thoipapy
import thoipapy.common
import thoipapy.experimental_data.add_experimental_data_to_train_set
import thoipapy.experimental_data.closest_heavy_atom_dist
import thoipapy.experimental_data.remove_hetero_contacts
import thoipapy.experimental_data.ttest_features
import thoipapy.feature_importance.anova
import thoipapy.feature_importance.ensemble_rfe
import thoipapy.feature_importance.mean_decrease_accuracy
import thoipapy.feature_importance.mean_decrease_impurity
import thoipapy.feature_importance.merge
import thoipapy.feature_importance.plots
import thoipapy.feature_importance.remove_duplicates
import thoipapy.features.combine_features
import thoipapy.features.entropy
import thoipapy.features.freecontact
import thoipapy.features.lipophilicity
import thoipapy.features.lips
import thoipapy.features.motifs
import thoipapy.features.physical_parameters
import thoipapy.features.preddimer_tmdock
import thoipapy.features.pssm
import thoipapy.features.rate4site
import thoipapy.features.relative_position
import thoipapy.homologues.colabfold_download
import thoipapy.homologues.colabfold_parser
import thoipapy.homologues.NCBI_download
import thoipapy.homologues.NCBI_parser
import thoipapy.ML_model.train_model
import thoipapy.ML_model.tune
import thoipapy.paper_figures.calc_PREDDIMER_TMDOCK_closedist
import thoipapy.paper_figures.combine_BOcurve_files
import thoipapy.paper_figures.create_heatmap_from_merge_file
import thoipapy.paper_figures.retrospective
import thoipapy.validation.combine_mult_predictors
import thoipapy.validation.gather
import thoipapy.validation.indiv_validation
import thoipapy.validation.leave_one_out
import thoipapy.validation.multiple_predictors
import thoipapy.validation.precision_recall
import thoipapy.validation.random_interface
import thoipapy.validation.roc
import thoipapy.validation.tenfold
import thoipapy.validation.testset_trainset
from thoipapy.artefacts import ArtefactPaths
from thoipapy.clustering.pairwise_aln_similarity_matrix import create_identity_matrix_from_protein_set
from thoipapy.paths import RUN_SETTINGS_CSV
from thoipapy.run_settings import RunSettings
from thoipapy.utils import get_test_and_train_set_lists, get_testsetname_trainsetname_from_run_settings


def run(s: dict):
    """Run the full THOIPA pipeline for every protein set named in the settings.

    Parameters
    ----------
    s : dict
        Settings dictionary, as returned by thoipapy.common.create_settingdict.
    """
    # if multiple sets need to be run, split them by comma
    if isinstance(s["set_number"], str) and "," in s["set_number"]:
        list_protein_sets = [int(n) for n in s["set_number"].split(",")]
    else:
        list_protein_sets = [s["set_number"]]

    for set_number in list_protein_sets:
        run_one_set(s, set_number)


def run_one_set(s: dict, set_number: int):
    """Run the pipeline for a single protein set.

    Split out of run() so that the settings are not rewritten mid-run. The loop used to assign
    s["set_number"] and s["setname"] back into the shared dictionary on every iteration, which
    made s both configuration and scratch space and meant it could never be frozen or reused.

    Parameters
    ----------
    s : dict
        Settings dictionary.
    set_number : int
        Number of the protein set to process, e.g. 8 for set08.
    """
    setname = f"set{set_number:02d}"
    # Work on a per-set copy. The loop previously assigned set_number and setname back into the
    # caller's dictionary, so the settings object a caller passed in came back modified and the
    # value of s["setname"] depended on how far the loop had got.
    s = {**s, "set_number": set_number, "setname": setname}
    paths = ArtefactPaths.from_settings(s)
    stages = RunSettings.from_settings(s)
    sets_dir = os.path.join(s["sets_dir"])
    # create a results folder for that set
    if not os.path.isdir(os.path.join(s["data_dir"], "results", setname)):
        os.makedirs(os.path.join(s["data_dir"], "results", setname))

    logging = thoipapy.common.setup_keyboard_interrupt_and_error_logging(s, setname)
    logging.info(f"STARTING PROCESSING OF {setname}.")

    set_path = thoipapy.common.get_path_of_protein_set(setname, sets_dir)

    ##############################################################################################
    #                                                                                            #
    #                     open and process a set of protein sequences                            #
    #                                                                                            #
    ##############################################################################################
    # load the protein set (e.g. set08_train.csv) as a dataframe
    df_set = pd.read_csv(set_path)

    # create list of uniprot accessions to run
    acc_list = df_set.acc.tolist()
    logging.info(f"settings file: {Path(s['settings_path']).name}")
    logging.info(f"stages enabled: {', '.join(stages.enabled_stages())}")
    logging.info(f"protein set {set_number}, {len(acc_list)} proteins: {acc_list}")

    dfset = thoipapy.common.process_set_protein_seqs(s, setname, df_set, set_path)

    # create a database label. Either crystal, NMR, ETRA or "mixed"
    subsets = df_set["database"].unique()
    if len(subsets.shape) == 1:
        subsets[0]
    else:
        pass

    if stages.create_identity_matrix_from_set_seqs:
        if True:
            create_identity_matrix_from_protein_set(paths, logging)

    ###################################################################################################
    #                                                                                                 #
    #                  calculate closedistance from NMR and crystal structures                        #
    #                                                                                                 #
    ###################################################################################################

    # DEPRECATED. Use atom_dist module instead to get closest heavy-atom distances
    # if stages.Get_Tmd_Homodimers :
    #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.download_xml_get_alphahelix_get_homo_pair(s, logging)
    #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.Download_trpdb_Calc_inter_rr_pairs(s, logging)
    #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.create_redundant_interact_homodimer_rm_shorttm(s, logging)
    #     #thoipapy.structures.deprecated.get_tmd_nr_homodimer.extract_crystal_resolv035_interact_pairs_and_create_fasta_file(s, logging)
    #     thoipapy.structures.deprecated.get_tmd_nr_homodimer.create_multiple_bind_closedist_file(s, logging)
    #     pass

    if stages.retrospective_coevolution:
        # thoipapy.paper_figures.retrospective.calc_retrospective_coev_from_list_interf_res(s, dfset, logging)
        thoipapy.paper_figures.retrospective.calc_retrospective_coev_from_struct_contacts(s, dfset, logging)

    # if stages.calc_NMR_closedist :
    #    thoipapy.structures.deprecated.NMR_data.calc_closedist_from_NMR_best_model(s)

    if stages.Atom_Close_Dist:
        # This stage needs a helix-pair file listing the PDB chains and their TM boundaries. No
        # setting supplies one, and the call here passed the thoipapy module and the settings dict
        # into a function expecting ArtefactPaths and a file path, so it could never have run.
        raise NotImplementedError(
            "Atom_Close_Dist needs a homodimer helix-pair file, which no setting provides. "
            "Supply the path and call homodimer_residue_closedist_calculate_from_complex directly."
        )

    ###################################################################################################
    #                                                                                                 #
    #                   homologues download from NCBI. parse, filter and save                         #
    #                                                                                                 #
    ###################################################################################################

    if stages.run_retrieve_NCBI_homologues_with_blastp:
        thoipapy.homologues.NCBI_download.download_homologues_from_ncbi_mult_prot(
            paths, df_set, s["expect_value"], s["hit_list_size"], s["rerun_existing_blast_results"], logging
        )

    if stages.run_parse_homologues_xml_into_csv:
        thoipapy.homologues.NCBI_parser.parse_NCBI_xml_to_csv_mult_prot(
            paths, df_set, s["surres"], s["e_value_cutoff"], logging
        )

    if stages.run_retrieve_homologues_from_colabfold:
        thoipapy.homologues.colabfold_download.download_homologues_from_colabfold_mult_prot(
            paths, df_set, s["colabfold_mode"], s["rerun_existing_blast_results"], logging
        )

    if stages.run_parse_colabfold_a3m_into_csv:
        thoipapy.homologues.colabfold_parser.parse_a3m_to_csv_mult_prot(
            paths,
            df_set,
            s["e_value_cutoff"],
            thoipapy.homologues.colabfold_download.mode_searches_env_db(s["colabfold_mode"]),
            logging,
        )

    if stages.parse_csv_homologues_to_alignment:
        thoipapy.homologues.NCBI_parser.extract_filtered_csv_homologues_to_alignments_mult_prot(
            paths, df_set, s["max_n_gaps_in_TMD_query_seq"], s["min_identity_of_TMD_seq"], logging
        )

    ###################################################################################################
    #                                                                                                 #
    #                   machine learning feature calculation                                             #
    #                                                                                                 #
    ###################################################################################################

    if stages.pssm_calculation:
        thoipapy.features.pssm.create_PSSM_from_MSA_mult_prot(paths, df_set, logging)

    if stages.entropy_calculation:
        thoipapy.features.entropy.entropy_calculation_mult_prot(paths, df_set, s["surres"], logging)

    if stages.rate4site_calculation:
        thoipapy.features.rate4site.rate4site_calculation_mult_prot(paths, df_set, logging)

    if stages.coevolution_calculation:
        if "Windows" in platform.system():
            sys.stdout.write(
                "\n Freecontact cannot be run in Windows! Skipping coevolution_calculation_with_freecontact_mult_prot."
            )
            thoipapy.features.freecontact.parse_freecontact_coevolution_mult_prot(paths, df_set, logging)
        else:
            thoipapy.features.freecontact.coevolution_calculation_with_freecontact_mult_prot(paths, df_set, logging)
            thoipapy.features.freecontact.parse_freecontact_coevolution_mult_prot(paths, df_set, logging)

    if stages.clac_relative_position:
        thoipapy.features.relative_position.calc_relative_position_mult_prot(paths, df_set, s["surres"], logging)

    if stages.calc_lipo_from_pssm:
        thoipapy.features.lipophilicity.lipo_from_pssm_mult_prot(paths, df_set, s["lipophilicity_scale"], logging)

    if stages.lips_score_calculation:
        thoipapy.features.lips.LIPS_score_calculation_mult_prot(paths, df_set, logging)
        thoipapy.features.lips.parse_LIPS_score_mult_prot(paths, df_set, logging)

    if stages.motifs_from_seq:
        thoipapy.features.motifs.motifs_from_seq_mult_protein(paths, df_set, logging)

    if stages.combine_feature_into_train_data:
        thoipapy.features.combine_features.combine_all_features_mult_prot(
            paths, df_set, s["lipophilicity_scale"], s["surres"], logging
        )
        thoipapy.features.physical_parameters.add_physical_parameters_to_features_mult_prot(paths, df_set, logging)
        thoipapy.experimental_data.add_experimental_data_to_train_set.add_experimental_data_to_combined_features_mult_prot(
            paths, df_set, s["inter_pair_max"], logging
        )
        if stages.generate_randomised_interfaces:
            thoipapy.validation.random_interface.add_random_interface_to_combined_features_mult_prot(
                paths, df_set, s["inter_pair_max"], logging
            )
        if stages.add_PREDDIMER_TMDOCK_to_combined_features:
            if True:
                thoipapy.features.preddimer_tmdock.add_PREDDIMER_TMDOCK_to_combined_features_mult_prot(
                    paths, df_set, logging
                )
        if stages.remove_crystal_hetero:
            thoipapy.experimental_data.remove_hetero_contacts.remove_crystal_hetero_contact_residues_mult_prot(
                paths, df_set, logging
            )
        thoipapy.features.combine_features.combine_all_train_data_for_machine_learning(
            paths, df_set, stages.remove_crystal_hetero, logging
        )

    ###################################################################################################
    #                                                                                                 #
    #                                    model validation                                             #
    #                                                                                                 #
    ###################################################################################################
    if stages.run_feature_selection:
        thoipapy.feature_importance.mean_decrease_impurity.get_initial_ensemble_parameters_before_feature_selection(
            paths, s["bind_column"], s["bootstrap"], logging
        )
        thoipapy.feature_importance.mean_decrease_impurity.calc_feat_import_using_MDI_before_feature_seln(
            paths, s["bootstrap"], logging
        )
        thoipapy.feature_importance.remove_duplicates.remove_duplicate_features_with_lower_MDI(
            paths,
            s["bind_column"],
            s["features_to_be_retained_during_selection"],
            s["max_similarity_duplicate_features"],
            logging,
        )
        thoipapy.feature_importance.anova.select_best_features_with_anova(
            paths, s["bind_column"], s["n_top_features_to_keep"], logging
        )
        thoipapy.feature_importance.ensemble_rfe.select_best_features_with_ensemble_rfe(
            paths, s["bind_column"], s["n_top_features_to_keep"], s["bootstrap"], logging
        )
        thoipapy.feature_importance.merge.merge_top_features_anova_ensemble(
            paths, s["bind_column"], s["features_to_be_retained_during_selection"], logging
        )
        thoipapy.ML_model.tune.tune_ensemble_parameters_after_feature_seln(paths, s["bind_column"], logging)

    if stages.calc_feature_importances:
        thoipapy.feature_importance.mean_decrease_accuracy.calc_feat_import_from_mean_decrease_accuracy(
            paths,
            s["bind_column"],
            s["min_n_homol_training"],
            s["bootstrap"],
            s["n_residues_AUBOC_validation"],
            logging,
        )
        thoipapy.feature_importance.plots.plot_feature_importance(paths, s["bind_column"], s["bootstrap"], logging)

    if stages.conduct_ttest:
        _, ttest_trainsetname = get_testsetname_trainsetname_from_run_settings(s)
        thoipapy.experimental_data.ttest_features.conduct_ttest_for_all_features(
            paths, s["bind_column"], ttest_trainsetname, logging
        )

    if stages.train_machine_learning_model:
        thoipapy.ML_model.train_model.train_machine_learning_model(
            paths, s["bind_column"], s["min_n_homol_training"], s["bootstrap"], logging
        )

    if stages.run_testset_trainset_validation:
        thoipapy.validation.testset_trainset.run_testset_trainset_validation(
            paths, s["n_residues_AUBOC_validation"], *get_test_and_train_set_lists(s), logging
        )

    ###################################################################################################
    #                                                                                                 #
    #                                               figures                                           #
    #                                                                                                 #
    ###################################################################################################

    if stages.compare_selected_predictors:
        thoipapy.paper_figures.combine_BOcurve_files.compare_selected_predictors(s, logging)

    if stages.calc_PREDDIMER_TMDOCK_closedist:
        thoipapy.paper_figures.calc_PREDDIMER_TMDOCK_closedist.calc_closedist_from_PREDDIMER_TMDOCK_best_model(
            s, df_set, logging
        )

    if stages.run_validation:
        sys.stdout.write("\n--------------- starting run_validation ---------------\n")
        namedict = thoipapy.utils.create_namedict(paths.protein_names_csv)
        THOIPA_predictor_name = "THOIPA_{}_LOO".format(s["set_number"])
        predictors = [THOIPA_predictor_name, "PREDDIMER", "TMDOCK", "LIPS_surface_ranked", "random"]
        testsetname, trainsetname = get_testsetname_trainsetname_from_run_settings(s)
        if s["setname"] == testsetname:
            predictors.append(f"thoipa.train{trainsetname}")

        thoipapy.validation.tenfold.run_10fold_cross_validation(
            paths, s["min_n_homol_training"], s["cross_validation_number_of_splits"], s["bootstrap"], logging
        )
        thoipapy.validation.tenfold.create_10fold_cross_validation_fig(
            paths, s["cross_validation_number_of_splits"], logging
        )

        thoipapy.validation.leave_one_out.run_LOO_validation(
            paths,
            df_set,
            s["bind_column"],
            s["min_n_homol_training"],
            s["bootstrap"],
            s["n_residues_AUBOC_validation"],
            s["use_multiprocessing"],
            s["multiple_tmp_simultaneous"],
            logging,
        )
        thoipapy.validation.leave_one_out.create_LOO_validation_fig(
            paths, df_set, s["n_residues_AUBOC_validation"], logging
        )

        thoipapy.validation.combine_mult_predictors.merge_predictions(paths, df_set, logging, testsetname, trainsetname)

        thoipapy.validation.indiv_validation.collect_indiv_validation_data(
            paths,
            df_set,
            s["n_residues_AUBOC_validation"],
            logging,
            namedict,
            predictors,
            THOIPA_predictor_name,
            subsets,
        )
        thoipapy.validation.indiv_validation.create_indiv_validation_figs(
            paths, logging, namedict, predictors, THOIPA_predictor_name, subsets
        )

        thoipapy.validation.multiple_predictors.validate_multiple_predictors_and_subsets_auboc(
            paths, df_set, s["n_residues_AUBOC_validation"], logging
        )
        thoipapy.validation.multiple_predictors.validate_multiple_predictors_and_subsets_auc(paths, df_set, logging)

        thoipapy.validation.roc.create_ROC_all_residues(paths, df_set, testsetname, trainsetname, logging)
        thoipapy.validation.precision_recall.create_precision_recall_all_residues(
            paths, df_set, testsetname, trainsetname, logging
        )

        thoipapy.validation.gather.gather_validation_data_for_figs(paths, df_set, logging)
        sys.stdout.write("\n--------------- finished run_validation ---------------\n")

    if stages.create_merged_heatmap_for_trainset_and_testset:
        thoipapy.paper_figures.create_heatmap_from_merge_file.create_merged_heatmap_for_trainset_and_testset(
            s, df_set, logging
        )

    if stages.download_10_homologues_from_ncbi:
        thoipapy.homologues.NCBI_download.download_10_homologues_from_ncbi(
            paths, df_set, s["rerun_existing_blast_results"], logging
        )

    if stages.plot_coev_vs_res_dist:
        thoipapy.paper_figures.retrospective.calc_coev_vs_res_dist(s, dfset, logging)
        thoipapy.paper_figures.retrospective.plot_coev_vs_res_dist(s, logging)

    # close the logger. A new one will be made for the next protein list.
    logging.info(f"FINISHED PROCESSING OF {setname}.")
    logging.shutdown()


if __name__ == "__main__":
    # Hardcoded settings path, no argparse, per the house style for pipeline scripts: strict and
    # reproducible rather than flexible. To run a different settings file, import run() and pass
    # your own settings dict.
    settings_dict = thoipapy.common.create_settingdict(RUN_SETTINGS_CSV)
    run(settings_dict)
