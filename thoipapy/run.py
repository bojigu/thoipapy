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

import argparse
import sys
import os
import thoipapy
from thoipapy import common
import csv
import platform

# read the command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "\Path\to\your\settingsfile.xlsx"')
parser.add_argument("-i",  # "-setting input fasta file location",
                    help=r'Full path to your input file.'
                         r'E.g. "\Path\to\your\input.fasta"')
parser.add_argument("-tmd",  # "-setting input fasta file location",
                    help=r'Full path to your input file contain the tmd sequence.'
                         r'E.g. "\Path\to\your\P01908_tmd.txt"')
parser.add_argument("-ts",  # "-setting tm start",
                    help=r'integere tm start value'
                         r'E.g. "219"')
parser.add_argument("-te",  # "-setting tm end ",
                    help=r'integer tm end value.'
                         r'E.g. "231"')
parser.add_argument("-of",  # "-setting output file path",
                    help=r'Full path to your prediction output file.'
                         r'E.g. "\Path\to\your output_file\"')
parser.add_argument("-email_to",  # "-setting output email address",
                    help=r'user email given on web server'
                         r'E.g. "***REMOVED***"')

if __name__ == "__main__":

    sys.stdout.write('\nRun thoipapy as follows:')
    sys.stdout.write(r'python \Path\to\run.py -s \Path\to\your\settingsfile.xlsx -i \Path\to\your\input.fasta -tmd \Path\to\your\input_tmd.txt '
          r'-ts tm_start_position -te tm_end_position -of C:\Path\to\your output_file\ -email_to email_address')

    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    set_=common.create_settingdict(args.s)

    # create output file, parsed file, and the output figure file names.
    if args.of:
        output_file_loc=os.path.join(args.of,"output.csv")
        output_parse_file=os.path.join(args.of,"output_parse.csv")
        output_png_loc=os.path.join(args.of,"output.png")

    tm_protein_name=set_["tm_protein_name"]
    Data_type=set_["Datatype"]

    logging=common.setup_keyboard_interrupt_and_error_logging(set_,tm_protein_name)

    tmp_lists=thoipapy.proteins.get_tmp_lists.extract_tmps_from_input_file(set_)
    test_tmp_lists=thoipapy.proteins.get_tmp_lists.extract_test_tmps_from_input_file(set_)

    set_["tm_protein_name"]='input'
    if args.i is not None:
        set_["input_fasta_file"]=args.i
        set_["tm_len"] = thoipapy.common.calculate_fasta_file_length(set_)
    if args.ts is not None:
        set_["tm_start"]=args.ts
    if args.te is not None:
        set_["tm_end"]=args.te
    if args.tmd is not None:
        set_["input_tmd_file"]=args.tmd
        set_["tm_start"], set_["tm_end"] = common.tmd_positions_match_fasta(set_)
    if args.email_to is not None:
        set_["email_to"]=args.email_to


    list_number=int(set_["list_number"])

    if list_number == 42:
        set_["db"] = "Etra"
    elif list_number == 43:
        set_["db"]  = "Nmr"
    elif list_number == 44:
        set_["db"] = "Crystal"
    elif list_number == 45:
        set_["db"] = "All"


    if list_number == 99:
        set_["db"] = "Etra"
        set_["list_of_tmd_start_end"] = os.path.join(os.path.dirname(args.s), "Tmd_Start_End_List_Uniq_New_{}.csv".format("MT1"))
    else:
        set_["list_of_tmd_start_end"] = os.path.join(os.path.dirname(args.s), "Tmd_Start_End_List_Uniq_New_{}.csv".format(set_["db"]))

    # this is important, if user want to run multiple proteins simultaneously, user has to set the tmd start and end list file by themselves
    # example of the tmd input file would look like this:
    # Protein,TMD_Length,TMD_Start,TMD_End
    # O15455,904,705,722
    # P07174,425,253,273


    #set_["list_of_tmd_start_end"] = os.path.join(os.path.dirname(args.s), "Tmd_Start_End_List_Uniq_New_{}.csv".format(set_["db"]))
    # set_["list_of_tmd_start_end"] = os.path.join(set_["data_harddrive"], "Input_data",
    #                                              "Tmd_Start_End_List_Uniq_New_{}.csv".format(set_["db"]))

    # when only run one protein each time, set_["multiple_tmp_simultaneous"] is false, and create the query protein information file
    # according to the arguments inputed by user
    if not set_["multiple_tmp_simultaneous"]:
        query_protein_tmd_file = os.path.join(set_["Protein_folder"], "Query_Protein_Tmd.csv")
        query_protein_tmd_file_handle=open(query_protein_tmd_file,"w")
        writer = csv.writer(query_protein_tmd_file_handle, delimiter=',', quoting = csv.QUOTE_NONE,lineterminator='\n')
        writer.writerow(["Protein","TMD_Length","TMD_Start","TMD_End"])
        writer.writerow([set_["tm_protein_name"],set_["tm_len"],set_["tm_start"],set_["tm_end"]])
        query_protein_tmd_file_handle.close()
        set_["list_of_tmd_start_end"]=query_protein_tmd_file

    #create new fasta file by only keep tmd and surrounded 20 residues for future blastp work
    #this function works for both one query protein or multiple protiens simultaneously
    thoipapy.common.create_TMD_surround20_fasta_file(set_)



                    ###################################################################################################
                    #                                                                                                 #
                    #                   homologous download with hhblits                                              #
                    #                                                                                                 #
                    ###################################################################################################

    # this is not used, since we used NCBI blast instead
    if set_["run_retrieve_homologous_with_hhblits"]:
        thoipapy.hhblits.download.download_homologous_with_hhblits(set_, logging)

        thoipapy.hhblits.download.parse_a3m_alignment( set_, logging)

                    ###################################################################################################
                    #                                                                                                 #
                    #                   homologous download from NCBI                                                 #
                    #                                                                                                 #
                    ###################################################################################################


    if set_["run_retrieve_NCBI_homologous_with_blastp"]:
        thoipapy.NCBI_BLAST.download.download.download_homologous_from_ncbi(set_, logging)


                    ###################################################################################################
                    #                                                                                                 #
                    #                   convert homologous xml file to csv                                            #
                    #                                                                                                 #
                    ###################################################################################################


    if set_["run_parse_homologous_xml_into_csv"]:
        thoipapy.NCBI_BLAST.parse.parser.parse_NCBI_xml_to_csv(set_,logging)

    if set_["parse_csv_homologous_to_alignment"]:
        thoipapy.NCBI_BLAST.parse.parser.extract_filtered_csv_homologous_to_alignments(set_,logging)


                    ###################################################################################################
                    #                                                                                                 #
                    #                   Random Forest feature calculation                                             #
                    #                                                                                                 #
                    ###################################################################################################

    if set_["RandomForest_feature_calculation"]:

        #thoipapy.RF_features.feature_calculate.mem_a3m_homologous_filter(set_, logging)

        if set_["pssm_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.create_PSSM_from_MSA_mult_prot(set_, logging)

        if set_["calc_lipo_from_pssm"]:
            thoipapy.RF_features.feature_calculate.calc_lipo_from_pssm(set_,logging)


        if set_["entropy_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.entropy_calculation(set_, logging)

        if set_["cumulative_coevolution_feature_calculation"]:
            if "Windows" in platform.system():
                sys.stdout.write("\n Freecontact cannot be run in Windows! Skipping coevoluton_calculation_with_freecontact function.")
                thoipapy.RF_features.feature_calculate.cumulative_co_evolutionary_strength_parser(tmp_lists, tm_protein_name, thoipapy, set_, logging)
            else:
                thoipapy.RF_features.feature_calculate.coevoluton_calculation_with_freecontact(set_, logging)
                thoipapy.RF_features.feature_calculate.cumulative_co_evolutionary_strength_parser(tmp_lists, tm_protein_name,thoipapy ,set_, logging)

        if set_["clac_relative_position"]:
            thoipapy.RF_features.feature_calculate.relative_position_calculation(set_,logging)

        if set_["lips_score_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.Lips_score_calculation(tmp_lists,tm_protein_name, thoipapy, set_, logging)
            thoipapy.RF_features.feature_calculate.Lips_score_parsing( set_, logging)

        #thoipapy.RF_features.feature_calculate.convert_bind_data_to_csv(set_, logging)

        if set_["combine_feature_into_train_data"]:
            if set_["db"] == "Crystal" or set_["db"] == "Nmr":
                thoipapy.RF_features.feature_calculate.features_combine_to_traindata( set_, logging)
                thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_train_data( set_, logging)
            if set_["db"] == "Etra":
                thoipapy.RF_features.feature_calculate.features_combine_to_testdata( set_, logging)
                thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_test_data(set_, logging)
            thoipapy.RF_features.feature_calculate.combine_all_train_data_for_random_forest(set_,logging)

    if set_["run_random_forest"]:
        #thoipapy.RF_features.RF_Train_Test.RF_10flod_cross_validation(tmp_lists,thoipapy,set_,logging)
        thoipapy.RF_features.RF_Train_Test.run_Rscipt_random_forest(tmp_lists, thoipapy, set_, output_file_loc,logging)

    if set_["parse_prediciton_output"]:
        thoipapy.RF_features.Output_Parse.parse_Predicted_Output(thoipapy,set_,output_file_loc,output_parse_file,logging)

    if set_["Send_sine_curve_to_email"]:
        print('begining to run run sine curve fitting')
        thoipapy.Sine_Curve.SineCurveFit.Save_Sine_Curve_Result(set_,output_file_loc,output_png_loc)
        logging.info('the fitting of sine curve is done')

                    ###################################################################################################
                    #                                                                                                 #
                    #                  Bind residues calculation for train data                                       #
                    #                                                                                                 #
                    ###################################################################################################


    if set_["Atom_Close_Dist"]:
        infor=thoipapy.Atom_Dist.Residu_Closest_Dist.homodimer_residue_closedist_calculate_from_complex(thoipapy, set_, logging)
        print(infor)

                    ###################################################################################################
                    #                                                                                                 #
                    #                  Bind residues calculation for train data                                       #
                    #                                                                                                 #
                    ###################################################################################################

    if set_["Send_email_finished"]:
        thoipapy.Send_Email.Send_Email_Smtp.send_email_when_finished(set_, thoipapy,output_parse_file,output_png_loc)
