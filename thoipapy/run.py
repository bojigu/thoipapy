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
warnings.filterwarnings("ignore")
import argparse
import sys
import os
import thoipapy
from thoipapy import common
import csv
import platform
import pandas as pd
import re

# read the command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "\Path\to\your\settingsfile.xlsx"')

if __name__ == "__main__":

    sys.stdout.write('\nRun thoipapy as follows:')
    sys.stdout.write(r'python \Path\to\run.py -s \Path\to\your\settingsfile.xlsx -i \Path\to\your\input.fasta -tmd \Path\to\your\input_tmd.txt '
          r'-ts tm_start_position -te tm_end_position -of C:\Path\to\your output_file\ -email_to email_address')
    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    s = common.create_settingdict(args.s)

    ##############################################################################################
    #                                                                                            #
    #                               setname, logging, results folder                             #
    #                                                                                            #
    ##############################################################################################
    sets_folder = s["sets_folder"]

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
        s["set_results_folder"] = os.path.join(s["Results_folder"], setname)
        if not os.path.isdir(s["set_results_folder"]):
            os.makedirs(s["set_results_folder"])

        logging = common.setup_keyboard_interrupt_and_error_logging(s, setname)
        logging.info("STARTING PROCESSING OF {}.".format(setname))

        set_path = thoipapy.common.get_path_of_protein_set(setname, sets_folder)

        ##############################################################################################
        #                                                                                            #
        #                     open and process a set of protein sequences                            #
        #                                                                                            #
        ##############################################################################################
        # load the protein set (e.g. set01.xlsx) as a dataframe
        df_set = pd.read_excel(set_path, sheetname='proteins')

        # create list of uniprot accessions to run
        acc_list = df_set.acc.tolist()
        sys.stdout.write("settings file : {}\nsettings : {}\nprotein set number {}, acc_list : {}\n".format(os.path.basename(args.s), s, s["set_number"], acc_list))
        sys.stdout.flush()

        dfset = thoipapy.common.process_set_protein_seqs(s, setname, df_set, set_path)

        # create a database label. Either crystal, NMR, ETRA or "mixed"
        unique_database_labels = df_set["database"].unique()
        if len(unique_database_labels.shape) == 1:
            database_for_full_set = unique_database_labels[0]
        else:
            database_for_full_set = "mixed"

        ###################################################################################################
        #                                                                                                 #
        #                  calculate closedistance from NMR and crystal structures                        #
        #                                                                                                 #
        ###################################################################################################

        if s["Get_Tmd_Homodimers"] :

            #thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.Download_xml_Get_alphahelix_Get_homo_pair(s, logging)
            thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.Download_trpdb_Calc_inter_rr_pairs(s, logging)
            #thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.create_redundant_interact_homodimer_rm_shorttm(s, logging)
            #thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.extract_crystal_resolv035_interact_pairs_and_create_fasta_file(s, logging)
            thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.create_multiple_bind_closedist_file(s,logging)
            #thoipapy.Get_Tmd_Homodimer.get_tmd_nr_homodimer.create_average_fraction_DI_file(s, logging)

        if s["Atom_Close_Dist"]:
            infor = thoipapy.Atom_Dist.Residu_Closest_Dist.homodimer_residue_closedist_calculate_from_complex(thoipapy, s, logging)
            sys.stdout.write(infor)

        ###################################################################################################
        #                                                                                                 #
        #                   homologues download from NCBI. parse, filter and save                         #
        #                                                                                                 #
        ###################################################################################################

        if s["run_retrieve_NCBI_homologues_with_blastp"]:
            thoipapy.NCBI_BLAST.download.download.download_homologues_from_ncbi_mult_prot(s, df_set, logging)

        if s["run_parse_homologues_xml_into_csv"]:
            thoipapy.NCBI_BLAST.parse.parser.parse_NCBI_xml_to_csv_mult_prot(s, df_set, logging)

        if s["parse_csv_homologues_to_alignment"]:
            thoipapy.NCBI_BLAST.parse.parser.extract_filtered_csv_homologues_to_alignments_mult_prot(s, df_set, logging)


        ###################################################################################################
        #                                                                                                 #
        #                   Random Forest feature calculation                                             #
        #                                                                                                 #
        ###################################################################################################

        if s["pssm_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.create_PSSM_from_MSA_mult_prot(s, df_set, logging)

        if s["entropy_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.entropy_calculation_mult_prot(s, df_set, logging)

        if s["cumulative_coevolution_feature_calculation"]:
            if "Windows" in platform.system():
                sys.stdout.write("\n Freecontact cannot be run in Windows! Skipping coevolution_calculation_with_freecontact_mult_prot.")
                thoipapy.RF_features.feature_calculate.parse_freecontact_coevolution_mult_prot(s, df_set, logging)
            else:
                thoipapy.RF_features.feature_calculate.coevolution_calculation_with_freecontact_mult_prot(s, df_set, logging)
                thoipapy.RF_features.feature_calculate.parse_freecontact_coevolution_mult_prot(s, df_set, logging)

        if s["clac_relative_position"]:
            thoipapy.RF_features.feature_calculate.calc_relative_position_mult_prot(s, df_set, logging)

        if s["calc_lipo_from_pssm"]:
            thoipapy.RF_features.feature_calculate.lipo_from_pssm_mult_prot(s, df_set, logging)

        if s["lips_score_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.LIPS_score_calculation_mult_prot(s, df_set, logging)
            thoipapy.RF_features.feature_calculate.parse_LIPS_score_mult_prot(s, df_set, logging)

        if s["motifs_from_seq"]:
            thoipapy.RF_features.feature_calculate.motifs_from_seq_mult_protein(s, df_set, logging)

        if s["combine_feature_into_train_data"]:
            thoipapy.RF_features.feature_calculate.combine_all_features_mult_prot(s, df_set, logging)
            thoipapy.RF_features.feature_calculate.add_physical_parameters_to_features_mult_prot(s, df_set, logging)
            thoipapy.RF_features.feature_calculate.add_experimental_data_to_combined_features_mult_prot(s, df_set, logging)
            if "add_PREDDIMER_TMDOCK_to_combined_features" in s:
                if s["add_PREDDIMER_TMDOCK_to_combined_features"]:
                    thoipapy.RF_features.feature_calculate.add_PREDDIMER_TMDOCK_to_combined_features_mult_prot(s, df_set, logging)
            thoipapy.RF_features.feature_calculate.combine_all_train_data_for_random_forest(s, df_set, logging)

        ###################################################################################################
        #                                                                                                 #
        #                                    model validation                                             #
        #                                                                                                 #
        ###################################################################################################

        if s["run_10fold_cross_validation"]:
            thoipapy.RF_features.RF_Train_Test.run_10fold_cross_validation(s, logging)
            thoipapy.RF_features.RF_Train_Test.create_10fold_cross_validation_fig(s, logging)


        if s["run_LOO_validation"]:
            thoipapy.RF_features.RF_Train_Test.run_LOO_validation(s, df_set, logging)
            if "create_LOO_validation_figs" in s:
                if s["create_LOO_validation_figs"]:
                    thoipapy.RF_features.RF_Train_Test.create_LOO_validation_fig(s, df_set, logging)

        if s["calculate_variable_importance"]:
            thoipapy.RF_features.RF_Train_Test.calculate_variable_importance(s, logging)
            thoipapy.RF_features.RF_Train_Test.fig_variable_importance(s, logging)

        if s["train_random_forest_model"]:
            thoipapy.RF_features.RF_Train_Test.train_random_forest_model(s, logging)


        if s["run_testset_trainset_validation"] == True:
            thoipapy.figs.Create_Bo_Curve_files.run_testset_trainset_validation(s, logging)

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

        if s["FigZB_07"] == True:
            # barcharts of coevolution values for interface and non-interface
            thoipapy.figs.Average_Fraction_DI.FigZB_07(Fontsize, Width, Size, s)

        if s["FigZB_18"] == True:
            # heatmap of prediction from THOIPA, PREDDIMER, TMDOCK
            thoipapy.figs.Preddimer_TMdock_heatmap.FigZB_18(Fontsize,Width,Size)


        #if s["combine_BO_curve_files_HARDLINKED"] == True:
        #    thoipapy.figs.Combine_Bo_Curve_files.combine_BO_curve_files_HARDLINKED(s)

        # DEPRECATED. USE COMPARE PREDICTORS
        # if s["fig_plot_BO_curve_mult_train_datasets"] == True:
        #     thoipapy.figs.Combine_Bo_Curve_files.fig_plot_BO_curve_mult_train_datasets(s)

        if s["compare_predictors"] == True:
            thoipapy.figs.Combine_Bo_Curve_files.compare_predictors(s)

        if s["run_bocurve_comp_HARDLINKED"] == True:
            thoipapy.figs.BoCurve_ThoipaBest_comp_LIPS_and_Nmr.run_bocurve_comp_HARDLINKED(Fontsize, Width, Size, s, Linewidth)

        if s["calc_PREDDIMER_TMDOCK_closedist"] == True:
            thoipapy.figs.Calc_PREDDIMER_TMDOCK_Closedist.calc_closedist_from_PREDDIMER_TMDOCK_best_model(s)

        if s["add_predictions_to_combined_files"] == True:
            thoipapy.figs.combine_add_3_prediction.combine_file_add_PREDDIMER_TMDOCK_THOIPA_prediction(s, df_set, logging)

        if s["create_AUC_BoAuc_fig_bestset"] == True:
            thoipapy.figs.combine_add_3_prediction.create_AUC_BoAUC_figs_THOIPA_PREDDIMER_TMDOCK(s, df_set, logging)

        if s["create_merged_heatmap"] == True:
            thoipapy.figs.merged_heatmap.create_merged_heatmap(s, df_set, logging)

        if s["create_ROC_4predictors"] == True:
            thoipapy.figs.combine_add_3_prediction.create_ROC_Curve_comp_4predictors(s,df_set,logging)


        if "download_10_homologues_from_ncbi" in s:
            if s["download_10_homologues_from_ncbi"] == True:
                thoipapy.NCBI_BLAST.download.download.download_10_homologues_from_ncbi(s, df_set, logging)

        # close the logger. A new one will be made for the next protein list.
        logging.info("FINISHED PROCESSING OF {}.".format(setname))
        logging.shutdown()

    ###################################################################################################
    #                                                                                                 #
    #            DEPRECATED SINGLE PROTEIN STUFF, REPLACED BY run_THOIPA_prediction.py                #
    #                                                                                                 #
    ###################################################################################################
    # running_single_protein = False
    # if running_single_protein:
    #     """OLD PARSER ARGUMENTS
    #
    #     parser.add_argument("-i",  # "-setting input fasta file location",
    #                         help=r'Full path to your input file.'
    #                              r'E.g. "\Path\to\your\input.fasta"')
    #     parser.add_argument("-tmd",  # "-setting input fasta file location",
    #                         help=r'Full path to your input file contain the tmd sequence.'
    #                              r'E.g. "\Path\to\your\P01908_tmd.txt"')
    #     parser.add_argument("-ts",  # "-setting tm start",
    #                         help=r'integere tm start value'
    #                              r'E.g. "219"')
    #     parser.add_argument("-te",  # "-setting tm end ",
    #                         help=r'integer tm end value.'
    #                              r'E.g. "231"')
    #     parser.add_argument("-of",  # "-setting output file path",
    #                         help=r'Full path to your prediction output file.'
    #                              r'E.g. "\Path\to\your output_file\"')
    #     parser.add_argument("-email_to",  # "-setting output email address",
    #                         help=r'user email given on web server'
    #                              r'E.g. "***REMOVED***"')
    #     """
    #     # create output file, parsed file, and the output figure file names.
    #     if args.of:
    #         output_file_loc = os.path.join(args.of, "output.csv")
    #         output_parse_file = os.path.join(args.of, "output_parse.csv")
    #         output_png_loc = os.path.join(args.of, "output.png")
    #
    #     s["tm_protein_name"] = 'input'
    #     if args.i is not None:
    #         s["input_fasta_file"] = args.i
    #         s["tm_len"] = thoipapy.common.calculate_fasta_file_length(s)
    #     if args.ts is not None:
    #         s["TMD_start"] = args.ts
    #     if args.te is not None:
    #         s["TMD_end"] = args.te
    #     if args.tmd is not None:
    #         s["input_tmd_file"] = args.tmd
    #         s["TMD_start"], s["TMD_end"] = common.tmd_positions_match_fasta(s)
    #     if args.email_to is not None:
    #         s["email_to"] = args.email_to
    #
    #     # when only run one protein each time, s["multiple_tmp_simultaneous"] is false, and create the query protein information file
    #     # according to the arguments inputed by user
    #     if not s["multiple_tmp_simultaneous"]:
    #         query_protein_tmd_file = os.path.join(s["Protein_folder"], "Query_Protein_Tmd.csv")
    #         query_protein_tmd_file_handle = open(query_protein_tmd_file, "w")
    #         writer = csv.writer(query_protein_tmd_file_handle, delimiter=',', quoting=csv.QUOTE_NONE, lineterminator='\n')
    #         writer.writerow(["Protein", "TMD_len", "TMD_Start", "TMD_End"])
    #         writer.writerow([s["tm_protein_name"], s["tm_len"], s["TMD_start"], s["TMD_end"]])
    #         query_protein_tmd_file_handle.close()
    #         s["list_of_tmd_start_end"] = query_protein_tmd_file
    #
    #     # create new fasta file by only keep tmd and surrounded 20 residues for future blastp work
    #     # this function works for both one query protein or multiple proteins simultaneously
    #     # thoipapy.common.create_TMD_surround20_fasta_file(s)
    #
    #     if s["parse_prediction_output"]:
    #         thoipapy.RF_features.Output_Parse.parse_Predicted_Output(thoipapy,s,output_file_loc,output_parse_file,logging)
    #
    #     if s["Send_sine_curve_to_email"]:
    #         sys.stdout.write('begining to run run sine curve fitting')
    #         thoipapy.Sine_Curve.SineCurveFit.Save_Sine_Curve_Result(s,output_file_loc,output_png_loc)
    #         logging.info('the fitting of sine curve is done')
    #
    #     if s["Send_email_finished"]:
    #         thoipapy.Send_Email.Send_Email_Smtp.send_email_when_finished(s, thoipapy, output_parse_file, output_png_loc)

