#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         BO ZENG
Created:        Monday December 12 12:33:08 2016
Operation system required: Linux (currently not available for windows)
Dependencies:   Python 3.5
                numpy
                Bio
                freecontact (currently only availble in linux)
                pandas
Purpose:        Self-interacting single-pass membrane protein interface residues prediction
Credits:        All sections by Bo Zeng.
Further Details:

"""

import argparse
import os
import thoipapy
from thoipapy import common
import re

# read the command line arguments
parser = argparse.ArgumentParser()
# add only a single argument, the path to the settings file.
parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "C:\Path\to\your\settingsfile.xlsx"')
parser.add_argument("-i",  # "-setting input fasta file location",
                    help=r'Full path to your input file.'
                         r'E.g. "C:\Path\to\your\P01908.fasta"')
parser.add_argument("-tmd",  # "-setting input fasta file location",
                    help=r'Full path to your input file contain the tmd sequence.'
                         r'E.g. "C:\Path\to\your\P01908_tmd.txt"')
parser.add_argument("-ts",  # "-setting tm start",
                    help=r'integere tm start value'
                         r'E.g. "219"')
parser.add_argument("-te",  # "-setting tm end ",
                    help=r'integer tm end value.'
                         r'E.g. "231"')
parser.add_argument("-of",  # "-setting output file location",
                    help=r'Full path to your output file.'
                         r'E.g. "C:\Path\to\your\output.txt"')
parser.add_argument("-op",  # "-setting output file location",
                    help=r'Full path to your output picture file.'
                         r'E.g. "C:\Path\to\your\output.png"')
parser.add_argument("-email_to",  # "-setting output file location",
                    help=r'user email given on web server'
                         r'E.g. "***REMOVED***"')
if __name__ == "__main__":
#def run_thoipapy(excel_file_with_settings):
    print('\nRun thoipapy as follows:')
    print(r'python "C:\Path\to\run.py" -s "C:\Path\to\your\settingsfile.xlsx"')
    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    # convert the excel settings file to a python dictionary, s
    #s = korbinian.common.create_settingsdict(args.s)
    set_=common.create_settingdict(args.s)
    print(args.s)
    print(args.i)
    print(args.of)
    print(args.op)
    #print(set_)
    output_file_loc=args.of
    output_parse_file=output_file_loc + "1"
    #output_parse_file=output_file_loc
    output_png_loc=args.op
    tm_protein_name=set_["tm_protein_name"]
    Data_type=set_["Datatype"]

    logging=common.setup_keyboard_interrupt_and_error_logging(set_,tm_protein_name)
    #logging.warning("tm_protein_name : {}".format(tm_protein_name))

    base_filename_summaries = os.path.join(set_["homologous_folder"], '%06s' % tm_protein_name, 'List%06s' % tm_protein_name)
    excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')
    pathdict = common.create_pathdict(base_filename_summaries)
    #print(pathdict)
    tmp_lists=thoipapy.proteins.get_tmp_lists.extract_tmps_from_input_file(pathdict,set_,logging)
    test_tmp_lists=thoipapy.proteins.get_tmp_lists.extract_test_tmps_from_input_file(pathdict,set_,logging)
    set_["tm_protein_name"]="QuePro"
    set_["input_fasta_file"]=args.i
    set_["input_tmd_file"]=args.tmd
    set_["tm_start"]=args.ts
    set_["tm_end"]=args.te
    set_["email_to"]=args.email_to
    if args.tmd is not None:
       set_["tm_start"],set_["tm_end"] = common.tmd_positions_match_fasta(pathdict,set_,logging)
    #input_fasta_file_handle=open(args.i,"r")
    # for line in input_fasta_file_handle:
    #     if re.search("^>",line):
    #         if(line[4:10]=="P02724"):
    #             print("tm_protein_name : {}".format(line[4:10]))
    #             if set_["Send_email_finished"]:
    #                 output_file_loc="/scratch2/zeng/thoipapy/OutPut/P02724.csv"
    #                 thoipapy.Send_Email.Send_Email_Smtp.send_email_when_finished(set_, pathdict, output_file_loc)
    #         quit()

    #print(type(args.i[4:6]))
    #if args.i[4:6] == "P02724":
        #


    ###################################################################################################
    #                                                                                                 #
    #                   homologous download with hhblits                                              #
    #                                                                                                 #
    ###################################################################################################


    if set_["run_retrieve_homologous_with_hhblits"]:
        thoipapy.hhblits.download.download_homologous_with_hhblits(pathdict, set_, logging)

        thoipapy.hhblits.download.parse_a3m_alignment(pathdict, set_, logging)

        ###################################################################################################
        #                                                                                                 #
        #                   homologous download from NCBI                                                 #
        #                                                                                                 #
        ###################################################################################################


    if set_["run_retrieve_NCBI_homologous_with_blastp"]:
        thoipapy.NCBI_BLAST.download.download.download_homologous_from_ncbi(pathdict, set_, logging)


        ###################################################################################################
        #                                                                                                 #
        #                   convert homologous xml file to csv                                            #
        #                                                                                                 #
        ###################################################################################################


    if set_["run_parse_homologous_xml_into_csv"]:
        thoipapy.NCBI_BLAST.parse.parser.parse_NCBI_xml_to_csv(pathdict,tm_protein_name, set_,logging)

    if set_["parse_csv_homologous_to_alignment"]:
        thoipapy.NCBI_BLAST.parse.parser.extract_filtered_csv_homologous_to_alignments(tmp_lists,tm_protein_name,pathdict, set_,logging)


        ###################################################################################################
        #                                                                                                 #
        #                   Random Forest feature calculation                                             #
        #                                                                                                 #
        ###################################################################################################


    if set_["RandomForest_feature_calculation"]:

        #thoipapy.RF_features.feature_calculate.mem_a3m_homologous_filter(tmp_lists, tm_protein_name, pathdict, set_, logging)

        if set_["pssm_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.pssm_calculation(tmp_lists,tm_protein_name,pathdict,set_,logging)

        if set_["entropy_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.entropy_calculation(tmp_lists,tm_protein_name, pathdict, set_, logging)

        if set_["cumulative_coevolution_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.coevoluton_calculation_with_freecontact(tmp_lists,tm_protein_name, pathdict, set_, logging)
            thoipapy.RF_features.feature_calculate.cumulative_co_evolutionary_strength_parser(tmp_lists, tm_protein_name,pathdict, set_, logging)


        if set_["lips_score_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.Lips_score_calculation(tmp_lists,tm_protein_name, pathdict, set_, logging)
            thoipapy.RF_features.feature_calculate.Lips_score_parsing(tmp_lists,tm_protein_name, pathdict, set_, logging)

        #thoipapy.RF_features.feature_calculate.convert_bind_data_to_csv(tmp_lists, tm_protein_name, pathdict, set_, logging)

        if set_["combine_feature_into_train_data"]:
            #thoipapy.RF_features.feature_calculate.features_combine_to_traindata(tmp_lists,tm_protein_name, pathdict, set_, logging)
            #thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_train_data(tmp_lists, tm_protein_name, pathdict, set_, logging)
            thoipapy.RF_features.feature_calculate.features_combine_to_testdata(test_tmp_lists, tm_protein_name, pathdict, set_, logging)
            thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_test_data(test_tmp_lists, tm_protein_name, pathdict, set_, logging)
            #thoipapy.RF_features.feature_calculate.combine_all_train_data_for_random_forest(tmp_lists,tm_protein_name,pathdict,set_,logging)

    if set_["run_random_forest"]:
        #thoipapy.RF_features.RF_Train_Test.RF_10flod_cross_validation(tmp_lists,pathdict,set_,logging)
        thoipapy.RF_features.RF_Train_Test.run_Rscipt_random_forest(tmp_lists, pathdict, set_, output_file_loc,logging)

    if set_["parse_prediciton_output"]:
        thoipapy.RF_features.Output_Parse.parse_Predicted_Output(pathdict,set_,output_file_loc,output_parse_file,logging)

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
        infor=thoipapy.Atom_Dist.Residu_Closest_Dist.homodimer_residue_closedist_calculate_from_complex(pathdict, set_, logging)
        print(infor)

    if set_["Send_email_finished"]:
        thoipapy.Send_Email.Send_Email_Smtp.send_email_when_finished(set_, pathdict,output_parse_file,output_png_loc)
