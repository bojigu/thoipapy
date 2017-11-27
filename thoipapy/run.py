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
import glob
import pandas as pd
import re

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
    print("code get to this line")
    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    set_=common.create_settingdict(args.s)

    # create output file, parsed file, and the output figure file names.
    if args.of:
        output_file_loc=os.path.join(args.of,"output.csv")
        output_parse_file=os.path.join(args.of,"output_parse.csv")
        output_png_loc=os.path.join(args.of,"output.png")

    #tm_protein_name=set_["tm_protein_name"]
    #Data_type=set_["Datatype"]

    #tmp_lists=thoipapy.proteins.get_tmp_lists.extract_tmps_from_input_file(set_)
    #test_tmp_lists=thoipapy.proteins.get_tmp_lists.extract_test_tmps_from_input_file(set_)

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

    ##############################################################################################
    #                                                                                            #
    #                          setup routine currently copied from Yao module                    #
    #                                                                                            #
    ##############################################################################################
    set_["data_dir"] = os.path.join(set_["base_dir"], "data_xy")
    set_["figs_dir"] = os.path.join(set_["base_dir"], "figs")
    # get list of all excel files in sets folder
    #D:\\Dropbox\\tm_homodimer_dropbox\\sets
    #sets_folder = os.path.join(set_["base_dir"], "sets")
    sets_folder = set_["sets_folder"]
    xlsx_list = glob.glob(os.path.join(sets_folder, "*.xlsx"))

    # remove temporary open excel files from the list (hidden files that start with ~$)
    xlsx_list = [path for path in xlsx_list if r"~$" not in path]

    # define set name, which should be in the excel file name
    setname = "set{:02d}".format(set_["set_number"])
    # add to the dictionary itself
    set_["setname"] = setname
    # create a results folder for that set
    set_["set_results_folder"] = os.path.join(set_["Results_folder"], setname)
    if not os.path.isdir(set_["set_results_folder"]):
        os.makedirs(set_["set_results_folder"])

    logging = common.setup_keyboard_interrupt_and_error_logging(set_, setname)

    # get subset of excel files that contains e.g. "set01"
    matching_xlsx_file_list = [set_path for set_path in xlsx_list if setname in set_path]
    if len(matching_xlsx_file_list) == 1:
        set_path = matching_xlsx_file_list[0]
    elif len(matching_xlsx_file_list) == 0:
        raise FileNotFoundError("Excel file with this set not found.\nsetname = {}\nexcel files in folder = {}".format(setname, xlsx_list))
    elif len(matching_xlsx_file_list) > 1:
        raise ValueError("More than one excel file in set folder contains '{}' in the filename.\nexcel files in folder = {}".format(setname, xlsx_list))

    # load the protein set (e.g. set01.xlsx) as a dataframe
    df_set = pd.read_excel(set_path, sheetname='proteins', index_col=0)

    # create list of uniprot accessions to run
    acc_list = df_set.index.tolist()
    sys.stdout.write("settings file : {}\nsettings : {}\nlist number {}, acc_list : {}\n".format(os.path.basename(args.s), set_, set_["set_number"], acc_list))
    sys.stdout.flush()

    ##############################################################################################
    #                                                                                            #
    #                          process set of protein sequences                                  #
    #                                                                                            #
    ##############################################################################################

    df_set["seqlen"] = df_set.full_seq.str.len()
    df_set["TMD_len"] = df_set.TMD_seq.str.len()
    # df_set.loc["O75460", "TMD_seq"] = df_set.loc["O75460", "TMD_seq"] + "A"
    for acc in df_set.index:
        TMD_seq = df_set.loc[acc, "TMD_seq"]
        full_seq = df_set.loc[acc, "full_seq"]
        # use regeg to get indices for start and end of TMD in seq
        m = re.search(TMD_seq, full_seq)
        if m:
            # convert from python indexing to unprot indexing
            df_set.loc[acc, "TMD_start"] = m.start() + 1
            df_set.loc[acc, "TMD_end"] = m.end()
        else:
            raise IndexError("TMD seq not found in full_seq.\nacc = {}\nTMD_seq = {}\nfull_seq = {}".format(acc, TMD_seq, full_seq))
    num_of_sur_residues = set_["num_of_sur_residues"]
    df_set["TMD_start_pl_surr"] = df_set.TMD_start - num_of_sur_residues
    df_set.loc[df_set["TMD_start_pl_surr"] < 1, "TMD_start_pl_surr"] = 1
    df_set["TMD_end_pl_surr"] = df_set.TMD_end + num_of_sur_residues
    for acc in df_set.index:
        if df_set.loc[acc, "TMD_end_pl_surr"] > df_set.loc[acc, "seqlen"]:
            df_set.loc[acc, "TMD_end_pl_surr"] = df_set.loc[acc, "seqlen"]

    def slice_TMD_seq_pl_surr(df_set):
        # note that due to uniprot-like indexing, the start index = start-1
        return df_set['full_seq'][int(df_set['TMD_start_pl_surr'] - 1):int(df_set['TMD_end_pl_surr'])]

    df_set["TMD_seq_pl_surr"] = df_set.apply(slice_TMD_seq_pl_surr, axis=1)

    # add the number of included residues in the surrounding seq to the left and right of the TMD
    # e.g. 20 where the TMD is in the centre of the protein, otherwise <20 where TMD is near start or end of full seq
    df_set["tm_surr_left"] = df_set.TMD_start - df_set.TMD_start_pl_surr
    df_set["tm_surr_right"] = df_set.TMD_end_pl_surr - df_set.TMD_end

    """  Rearrange the dataframe columns so that the order is as follows.
    orig Bo file : ['acc', 'TMD_Length', 'TMD_Start', 'TMD_End', 'TMD_Sur_Left', 'TMD_Sur_Right']
    updated file = ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]

    """
    # reorder columns
    df_set = thoipapy.utils.set_column_sequence(df_set, ['seqlen', 'TMD_start', 'TMD_end', "tm_surr_left", "tm_surr_right", "database"])

    # convert the floats to integers
    df_set.iloc[:, 1:5] = df_set.iloc[:, 1:5].astype(int)

    # save to csv, which is opened by other functions
    list_of_tmd_start_end = os.path.join(set_["thoipapy_data_folder"], "Input_data", os.path.basename(set_path)[:-5] + "_processed.csv")
    set_["list_of_tmd_start_end"] = list_of_tmd_start_end
    df_set.to_csv(list_of_tmd_start_end)

    # create a database label. Either crystal, NMR, ETRA or "mixed"
    unique_database_labels = df_set["database"].unique()
    if len(unique_database_labels.shape) == 1:
        database = unique_database_labels[0]
    else:
        database = "mixed"

    # when only run one protein each time, set_["multiple_tmp_simultaneous"] is false, and create the query protein information file
    # according to the arguments inputed by user
    if not set_["multiple_tmp_simultaneous"]:
        query_protein_tmd_file = os.path.join(set_["Protein_folder"], "Query_Protein_Tmd.csv")
        query_protein_tmd_file_handle=open(query_protein_tmd_file,"w")
        writer = csv.writer(query_protein_tmd_file_handle, delimiter=',', quoting = csv.QUOTE_NONE,lineterminator='\n')
        writer.writerow(["Protein","TMD_len","TMD_Start","TMD_End"])
        writer.writerow([set_["tm_protein_name"],set_["tm_len"],set_["tm_start"],set_["tm_end"]])
        query_protein_tmd_file_handle.close()
        set_["list_of_tmd_start_end"]=query_protein_tmd_file

    #create new fasta file by only keep tmd and surrounded 20 residues for future blastp work
    #this function works for both one query protein or multiple protiens simultaneously

    #thoipapy.common.create_TMD_surround20_fasta_file(set_)


                    ###################################################################################################
                    #                                                                                                 #
                    #                   homologues download with hhblits                                              #
                    #                                                                                                 #
                    ###################################################################################################

    # this is not used, since we used NCBI blast instead
    if set_["run_retrieve_homologous_with_hhblits"]:
        thoipapy.hhblits.download.download_homologues_with_hhblits(set_, logging)

        thoipapy.hhblits.download.parse_a3m_alignment(set_, logging)

                    ###################################################################################################
                    #                                                                                                 #
                    #                   homologues download from NCBI                                                 #
                    #                                                                                                 #
                    ###################################################################################################


    if set_["run_retrieve_NCBI_homologues_with_blastp"]:
        thoipapy.NCBI_BLAST.download.download.download_homologues_from_ncbi(set_, df_set, logging)


                    ###################################################################################################
                    #                                                                                                 #
                    #                   convert homologues xml file to csv                                            #
                    #                                                                                                 #
                    ###################################################################################################


    if set_["run_parse_homologues_xml_into_csv"]:
        thoipapy.NCBI_BLAST.parse.parser.parse_NCBI_xml_to_csv_mult_prot(set_, df_set, logging)

    if set_["parse_csv_homologues_to_alignment"]:
        thoipapy.NCBI_BLAST.parse.parser.extract_filtered_csv_homologues_to_alignments_mult_prot(set_, df_set, logging)


                    ###################################################################################################
                    #                                                                                                 #
                    #                   Random Forest feature calculation                                             #
                    #                                                                                                 #
                    ###################################################################################################

    if set_["RandomForest_feature_calculation"]:

        #thoipapy.RF_features.feature_calculate.mem_a3m_homologues_filter(set_, logging)

        if set_["pssm_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.create_PSSM_from_MSA_mult_prot(set_, df_set, logging)

        if set_["calc_lipo_from_pssm"]:
            thoipapy.RF_features.feature_calculate.calc_lipo_from_pssm(set_,df_set, logging)


        if set_["entropy_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.entropy_calculation(set_, logging)

        if set_["cumulative_coevolution_feature_calculation"]:
            if "Windows" in platform.system():
                sys.stdout.write("\n Freecontact cannot be run in Windows! Skipping coevoluton_calculation_with_freecontact function.")
                thoipapy.RF_features.feature_calculate.cumulative_co_evolutionary_strength_parser(thoipapy, set_, logging)
            else:
                thoipapy.RF_features.feature_calculate.coevoluton_calculation_with_freecontact(set_, logging)
                thoipapy.RF_features.feature_calculate.cumulative_co_evolutionary_strength_parser(thoipapy ,set_, logging)

        if set_["clac_relative_position"]:
            thoipapy.RF_features.feature_calculate.relative_position_calculation(set_,logging)

        if set_["lips_score_feature_calculation"]:
            thoipapy.RF_features.feature_calculate.LIPS_score_calculation_mult_prot(thoipapy, set_, logging)
            thoipapy.RF_features.feature_calculate.Lips_score_parsing( set_, logging)

        #thoipapy.RF_features.feature_calculate.convert_bind_data_to_csv(set_, logging)

        if set_["combine_feature_into_train_data"]:
            if database == "Crystal" or database == "Nmr":
                thoipapy.RF_features.feature_calculate.features_combine_to_traindata( set_, logging)
                thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_train_data(set_, logging)
            if database == "Etra":
                thoipapy.RF_features.feature_calculate.features_combine_to_testdata( set_, logging)
                thoipapy.RF_features.feature_calculate.adding_physical_parameters_to_test_data(set_, logging)
            thoipapy.RF_features.feature_calculate.combine_all_train_data_for_random_forest(set_,logging)

    if set_["run_random_forest"]:
        thoipapy.RF_features.RF_Train_Test.thoipa_rfmodel_create(set_,logging)
        #thoipapy.RF_features.RF_Train_Test.RF_10flod_cross_validation(thoipapy,set_,logging)
        #thoipapy.RF_features.RF_Train_Test.run_Rscipt_random_forest(set_, output_file_loc, logging)



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
