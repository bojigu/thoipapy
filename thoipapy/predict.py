import warnings

import thoipapy.features as tf

warnings.filterwarnings("ignore")

import argparse
import glob
import hashlib
import os
import platform
import re
import sys
from io import StringIO
from pathlib import Path
import stat

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from django.utils.text import slugify
import joblib

import thoipapy
from thoipapy.homologues.NCBI_download import download_homologues_from_ncbi
from thoipapy.homologues.NCBI_parser import parse_NCBI_xml_to_csv, extract_filtered_csv_homologues_to_alignments
from thoipapy.utils import normalise_between_2_values, open_csv_as_series

# set matplotlib backend to Agg when run on a server
if os.environ.get('DISPLAY','') == '':
    sys.stdout.write('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

def run_THOIPA_prediction(protein_name, md5, TMD_seq, full_seq, out_dir, create_heatmap=True):
    """Function to run standalone THOIPA prediction for a protein transmembrane domain of interest.

    Parameters
    ----------
    protein_name : str
        Protein name
    TMD_seq : str
        Transmembrane domain sequence
    full_seq : str
        Full protein sequence
    out_dir : str
        Path to folder where the output will be saved

    Saved files and figures
    -----------------------


    Usage
    -----
    import thoipapy
    protein_name = "1c17_A"
    TMD_seq = "AAVMMGLAAIGAAIGIGILG"
    full_seq = "MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA"
    predictions_folder = "/path/to/your/output/folder"
    thoipapy.run_THOIPA_prediction(protein_name, TMD_seq, full_seq, predictions_folder)
    """
    # create settings dict from xlsx file
    thoipapy_module_path = Path(thoipapy.__file__).parent
    settings_path = thoipapy_module_path / "setting/thoipapy_standalone_run_settings.xlsx"
    s = thoipapy.common.create_settingdict(settings_path)

    #folder_with_md5_exists = False
    #if folder_with_md5_exists:
    #    create_indiv_name_and_email = 99
    #    sys.stdout.write("You shuld just email the results here.")
    #acc = md5[0:6]
    #out_dir = os.path.join(predictions_folder, md5)

    ###################################################################################################
    #                                                                                                 #
    #                                     setup all file paths                                        #
    #                                                                                                 #
    ###################################################################################################
    out_dir = Path(out_dir)
    datafiles_dir = out_dir / "datafiles"
    thoipapy.utils.make_sure_path_exists(datafiles_dir)

    # various residue features needs a protein folder as output
    blast_xml_file = datafiles_dir / "BLAST_results.xml"
    xml_txt = datafiles_dir / "BLAST_results_details.txt"
    xml_tar_gz = datafiles_dir / "BLAST_results.xml.tar.gz"
    BLAST_csv_tar = datafiles_dir / "BLAST_results.csv.tar.gz"
    fasta_all_TMD_seqs = datafiles_dir /"homologues.redundant.fas"
    path_uniq_TMD_seqs_for_PSSM_FREECONTACT = datafiles_dir / "homologues.uniq.for_PSSM_FREECONTACT.txt"
    path_uniq_TMD_seqs_no_gaps_for_LIPS = datafiles_dir / "homologues.uniq.for_LIPS.txt"
    path_uniq_TMD_seqs_surr5_for_LIPO = datafiles_dir / "homologues.uniq.for_LIPO.txt"
    fasta_uniq_TMD_seqs_surr5_for_LIPO = datafiles_dir / "homologues.uniq.for_LIPO.fas"
    pssm_csv = datafiles_dir / "pssm.csv"
    pssm_surr5_csv = datafiles_dir / "pssm_surr5.csv"
    lipo_csv = datafiles_dir / "lipo.csv"
    entropy_file = datafiles_dir / "entropy.csv"
    rate4site_csv = datafiles_dir / "rate4site.csv"
    freecontact_file = datafiles_dir / "freecontact_out.csv"
    freecontact_parsed_csv = datafiles_dir / "freecontact_parsed.csv"
    relative_position_file = datafiles_dir / "relative_position.csv"
    LIPS_output_file = datafiles_dir / "LIPS_output.csv"
    LIPS_parsed_csv = datafiles_dir / "LIPS_output_parsed.csv"
    motifs_file = datafiles_dir / "motifs.csv"
    full_seq_fasta_file = datafiles_dir / "protein.fasta"
    full_seq_phobius_output_file = datafiles_dir / "protein.phobius"
    feature_combined_file = datafiles_dir / "features_combined.csv"
    alignment_summary_csv = datafiles_dir / "homologues.alignment_summary.csv"
    THOIPA_full_out_csv = datafiles_dir / "THOIPA_full_out.csv"
    THOIPA_pretty_out_xlsx = out_dir / "THOIPA_out.xlsx"
    THOIPA_pretty_out_txt = out_dir / "THOIPA_out.csv"
    heatmap_path = out_dir / "heatmap.png"

    model_features_txt = thoipapy_module_path / "ML_model/model_features.txt"
    machine_learning_model = thoipapy_module_path / "ML_model/THOIPA_trained_ML_model.lpkl"

    logfile = out_dir / "logfile.txt"
    logging = thoipapy.common.setup_error_logging(logfile, "INFO", "INFO", print_system_info=False)

    logging.info("Starting THOIPA standalone prediction for {}.".format(protein_name))

    if os.path.isfile(THOIPA_full_out_csv):
        logging.info("{} already analysed. Previous results will be overwritten.".format(protein_name))

    ###################################################################################################
    #                                                                                                 #
    #                       check validity of the TMD and full sequence                               #
    #                                                                                                 #
    ###################################################################################################

    logging.info("md5 checksum of TMD and full sequence = {}".format(md5))

    seqlen = len(full_seq)
    hitlist = re.findall(TMD_seq, full_seq)
    if len(hitlist) > 1:
        logging.warning("TMD sequence is found multiple times in the full sequence.")
        return
    elif len(hitlist) == 0:
        logging.warning("TMD sequence was not found in full protein sequence. Please check input sequences.")
        return
        # the following code will only continue if the TMD seq is found ONCE in the full protein seq

    m = re.search(TMD_seq, full_seq)
    TMD_start = m.start()
    TMD_end = m.end()
    logging.info("TMD found in full sequence. start = {}, end = {}".format(TMD_start, TMD_end))

    ###################################################################################################
    #                                                                                                 #
    #                   get residues surrounding the TMD from the full sequence                       #
    #                                                                                                 #
    ###################################################################################################

    surr = s["num_of_sur_residues"]
    TMD_start_pl_surr, TMD_end_pl_surr = get_start_end_pl_surr(TMD_start, TMD_end, seqlen, surr=surr)
    TMD_start_pl_5, TMD_end_pl_5 = get_start_end_pl_surr(TMD_start, TMD_end, seqlen, surr=5)

    logging.info("TMD_start_pl_surr = {}, TMD_end_pl_surr = {}".format(TMD_start_pl_surr, TMD_end_pl_surr))

    # due to uniprot indexing, the TMD seq starts 1 residue earlier in the python index
    TMD_seq_pl_surr = full_seq[TMD_start_pl_surr: TMD_end_pl_surr]
    # currently I don't know why this should not be TMD_start_pl_5 - 1
    query_TMD_seq_surr5 = full_seq[TMD_start_pl_5: TMD_end_pl_5]
    logging.info("TMD_seq : {}".format(TMD_seq))
    logging.info("query_TMD_seq_surr5 : {}".format(query_TMD_seq_surr5))
    logging.info("TMD_seq_pl_surr : {}".format(TMD_seq_pl_surr))

    tm_surr_left = TMD_start - TMD_start_pl_surr
    tm_surr_right = TMD_end_pl_surr - TMD_end

    # reset to 5. AM NOT SURE WHY.
    if tm_surr_left >= 5:
        tm_surr_left_lipo = 5
    else:
        tm_surr_left_lipo = tm_surr_left
    if tm_surr_right >= 5:
        tm_surr_right_lipo = 5
    else:
        tm_surr_right_lipo = tm_surr_right

    ###################################################################################################
    #                                                                                                 #
    #                                    get BLAST homologues                                         #
    #                                                                                                 #
    ###################################################################################################
    # most scripts use uniprot accession as the protein name
    acc = protein_name
    # n_TMDs = rp.return_num_tmd(s, acc, full_seq, full_seq_fasta_file, full_seq_phobius_output_file, logging)

    expect_value = s["expect_value"]
    hit_list_size = s["hit_list_size"]
    if not os.path.isfile(xml_tar_gz):
        download_homologues_from_ncbi(acc, TMD_seq_pl_surr, blast_xml_file, xml_txt, xml_tar_gz, expect_value, hit_list_size, logging)

    if not os.path.isfile(BLAST_csv_tar):
        parse_NCBI_xml_to_csv(s, acc, xml_tar_gz, BLAST_csv_tar, TMD_start, TMD_end, logging)

    ###################################################################################################
    #                                                                                                 #
    #                              parse BLAST homologues to create MSA                               #
    #                                                                                                 #
    ###################################################################################################

    #if not os.path.isfile(path_uniq_TMD_seqs_for_PSSM_FREECONTACT):
    extract_filtered_csv_homologues_to_alignments(s, acc, len(TMD_seq), fasta_all_TMD_seqs, path_uniq_TMD_seqs_for_PSSM_FREECONTACT,
                                                      path_uniq_TMD_seqs_no_gaps_for_LIPS, path_uniq_TMD_seqs_surr5_for_LIPO, BLAST_csv_tar,
                                                      TMD_seq, query_TMD_seq_surr5, logging)

    tf.pssm.create_PSSM_from_MSA(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, pssm_csv, acc, TMD_seq, logging)
    tf.pssm.create_PSSM_from_MSA(path_uniq_TMD_seqs_surr5_for_LIPO, pssm_surr5_csv, acc, query_TMD_seq_surr5, logging)

    ###################################################################################################
    #                                                                                                 #
    #              get all necessary residue properties (machine learning features)                   #
    #                                                                                                 #
    ###################################################################################################

    tf.lipophilicity.lipo_from_pssm(acc, pssm_surr5_csv, lipo_csv, tm_surr_left_lipo, tm_surr_right_lipo, s["lipophilicity_scale"], logging, plot_linechart=True)

    tf.entropy.entropy_calculation(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, TMD_seq, entropy_file, logging)

    if "Windows" in platform.system():
        logging.warning("\n Freecontact, CD-HIT and rate4site cannot be run in Windows! Skipping coevolution_calculation_with_freecontact and rate4site_calculation.\n"
                        f"For testing, copy output files from Linux and rename to {freecontact_file} and {rate4site_csv}")
    else:
        tf.freecontact.coevolution_calculation_with_freecontact(path_uniq_TMD_seqs_for_PSSM_FREECONTACT, freecontact_file, s["freecontact_dir"], logging)
        tf.rate4site.rate4site_calculation(TMD_seq, acc, fasta_uniq_TMD_seqs_surr5_for_LIPO, rate4site_csv, logging, rerun_rate4site=False)

    tf.freecontact.parse_freecontact_coevolution(acc, freecontact_file, freecontact_parsed_csv, TMD_start, TMD_end, logging)

    tf.relative_position.calc_relative_position(acc, path_uniq_TMD_seqs_for_PSSM_FREECONTACT, relative_position_file, TMD_start, seqlen, logging)

    tf.lips.LIPS_score_calculation(path_uniq_TMD_seqs_no_gaps_for_LIPS, LIPS_output_file)

    tf.lips.parse_LIPS_score(acc, LIPS_output_file, LIPS_parsed_csv, logging)
    tf.motifs.motifs_from_seq(TMD_seq, TMD_seq_pl_surr, tm_surr_left, tm_surr_right, motifs_file, logging)

    database = "standalone_prediction"
    tf.combine_features.combine_all_features(s, full_seq, acc, database, TMD_seq, TMD_start, feature_combined_file, entropy_file, rate4site_csv, pssm_csv,
                                                            lipo_csv, freecontact_parsed_csv, relative_position_file, LIPS_parsed_csv, motifs_file,
                                                            alignment_summary_csv, full_seq_fasta_file, full_seq_phobius_output_file, logging)
    tf.physical_parameters.add_physical_parameters_to_features(acc, feature_combined_file, logging)

    ###################################################################################################
    #                                                                                                 #
    #             Apply trained machine learning model to get prediction values                       #
    #                                                                                                 #
    ###################################################################################################

    df_data = pd.read_csv(feature_combined_file, index_col=0)
    with open(model_features_txt, "r") as f:
        model_features = [x.strip() for x in f.readlines()]

    df_ML_input = df_data[model_features]

    cols_to_keep = ["residue_num", "residue_name", "res_num_full_seq"]

    df_out = pd.concat([df_data[cols_to_keep], df_ML_input], axis=1, join="outer")

    # load the trained machine learning model, saved in the thoipapy directory
    fit = joblib.load(machine_learning_model)

    # get prediction values for each residue
    df_out["THOIPA"] = fit.predict_proba(df_ML_input)[:, 1]

    ###################################################################################################
    #                                                                                                 #
    #         Prepare some "pretty" output files with THOIPA prediction values                        #
    #                                                                                                 #
    ###################################################################################################

    df_pretty_out = df_out[["residue_num", "residue_name", "THOIPA"]].copy()
    df_pretty_out.columns = ["residue number", "residue name", "THOIPA"]

    df_pretty_out["THOIPA"] = df_pretty_out["THOIPA"].round(3)
    df_pretty_out.to_excel(THOIPA_pretty_out_xlsx, index=False)

    # pad all content with spaces so it lines up with the column name
    df_pretty_out["residue number"] = df_pretty_out["residue number"].apply(lambda x: "{: >14}".format(x))
    df_pretty_out["residue name"] = df_pretty_out["residue name"].apply(lambda x: "{: >12}".format(x))
    df_pretty_out["THOIPA"] = df_pretty_out["THOIPA"].apply(lambda x : "{:>6.03f}".format(x))

    df_pretty_out.set_index("residue number", inplace=True)

    df_out.to_csv(THOIPA_full_out_csv)
    df_pretty_out.to_csv(THOIPA_pretty_out_txt, sep="\t")

    # print exactly what the CSV looks like
    out = StringIO()
    df_pretty_out.to_csv(out, sep="\t")
    logging.info("\n\nTHOIPA homotypic TMD interface prediction:\n\n{}".format(out.getvalue()))

    if create_heatmap:

        ###################################################################################################
        #                                                                                                 #
        #         Save a heatmap with THOIPA output, conservation, polarity, and coevolution              #
        #                                                                                                 #
        ###################################################################################################
        fontsize = 12
        tum_blue4_as_python_color = np.array([0, 82, 147]) / 255
        cmap = sns.light_palette(tum_blue4_as_python_color, as_cmap=True)

        cols_to_plot = ['THOIPA', 'conservation', 'relative_polarity', 'DImax']
        cols_to_plot_renamed = ['THOIPA', 'conservation', 'relative polarity', 'coevolution']

        # transpose dataframe so that "interface" etc is on the left
        df = df_out[cols_to_plot].copy()
        df.columns = cols_to_plot_renamed
        # start residue indexing from 1
        df.index = df.index + 1

        # normalise values so that they can all be plotted with same colour scheme and range of shading
        df["THOIPA"] = normalise_between_2_values(df["THOIPA"], 0.15, 0.5, invert=False)
        df["conservation"] = normalise_between_2_values(df["conservation"], 1.25, 3, invert=False)
        df["relative polarity"] = normalise_between_2_values(df["relative polarity"], 0.5, 2.5, invert=False)

        # transpose
        df = df.T

        """
        IMPORTANT!!
        The default fontsize controls the spacing between the subplots, EVEN IF THERE ARE NO TITLES or XLABELS!      
        """
        plt.rcParams['font.size'] = fontsize / 2

        # create plot
        fig, ax = plt.subplots(figsize=(16, 2))
        # duplicate plot so it's possible to add label at top
        ax2 = ax.twiny()
        ax.set_xticklabels(df.columns, rotation=0, fontsize=fontsize)

        # create heatmap
        sns.heatmap(df, ax=ax, cmap=cmap)

        # format residue numbering at bottom

        ax.tick_params(axis="x", direction='out', pad=1.5, tick2On=False)
        ax.set_xlabel("position in TMD", fontsize=fontsize)
        ax.set_yticklabels(df.index, rotation=0, fontsize=fontsize)

        # plot the same heatmap again, and format residue single-letter AA code at top
        sns.heatmap(df, ax=ax2, cmap=cmap)
        ax2.set_xlabel("{}, residue in TMD".format(protein_name), fontsize=fontsize)
        ax2.set_xticks(ax.get_xticks())
        ax2.xaxis.tick_top()
        ax2.set_xticklabels(df_out.residue_name, fontsize=fontsize)
        ax2.tick_params(axis="x", direction='out', pad=-0.1, tick2On=False)

        plt.tight_layout()
        fig.savefig(heatmap_path, dpi=240)
        fig.savefig(str(heatmap_path)[:-4] + ".pdf")
        logging.info("Output heatmap : {}".format(heatmap_path))

    logging.info("Output file : {}\n"
                 "THOIPA standalone completed successfully for {}.\n\n".format(THOIPA_pretty_out_xlsx, protein_name))

def get_start_end_pl_surr(TMD_start, TMD_end, seqlen, surr):
    """Get start and end of TMD plus surrounding sequence

    Parameters
    ----------
    TMD_start
    TMD_end
    seqlen
    surr

    Returns
    -------

    """
    TMD_start_pl_surr = TMD_start - surr
    if TMD_start_pl_surr < 1:
        TMD_start_pl_surr = 1
    TMD_end_pl_surr = TMD_end + surr
    if TMD_end_pl_surr > seqlen:
        TMD_end_pl_surr = seqlen
    return TMD_start_pl_surr, TMD_end_pl_surr


def get_md5_checksum(TMD_seq, full_seq):
    """Create md5 checksum from a concatemer of the TMD_seq and full_seq

    Parameters
    ----------
    TMD_seq : str
        TMD sequence
    full_seq : str
        Full protein sequence

    Returns
    -------
    md5 : str
        md5 checksum
    """
    TMD_plus_full_seq = TMD_seq + "_" + full_seq
    # adjust encoding for md5 creation
    #TMD_plus_full_seq = unicodedata.normalize('NFKD', TMD_plus_full_seq).encode('ascii', 'ignore')
    TMD_plus_full_seq = TMD_plus_full_seq.encode('ascii', 'ignore')
    hash_object = hashlib.md5(TMD_plus_full_seq)
    md5 = hash_object.hexdigest()
    return md5


def print_help():
    sys.stdout.write("\nUsage example:\n")
    sys.stdout.write(r"python thoipa.py -i D:\data\Q12983.txt -f D:\data\predictions")
    sys.stdout.write("\n\nOR process every input file in the -d input folder. "
                     "You can specify the output folder. Otherwise, a default 'output' folder will be "
                     "created in the same directory as the input folder. \n")
    sys.stdout.write(r"python thoipa.py -d D:\your\directory\with\input_text_files")
    sys.stdout.write("\n\n")
    sys.stdout.flush()

usage = "\nthoipa.py [-h] [-d D] [-i I] [-f F]\n\nexample: \npython thoipa.py -d /Path/to/your/input/csv/text/files"

# read the command line arguments
parser = argparse.ArgumentParser(usage=usage)

#parser.add_argument("-h", help = "Prints help")

parser.add_argument("-d",  # "-directory",
                    help=r'Path to directory containing input csv text files.')
parser.add_argument("-i",  # "-input_file",
                    help=r'Path to a single input csv text file, E.g. "/Path/to/your/file/Q12983.txt"')
parser.add_argument("-f",  # "-folder",
                    help='Optional path to an output folder.\nIf no output folder is specified, '
                         'a new output folder will be created\nin same directory as the input folder.')

if __name__ == "__main__":
    """
    Example input csv file:
    -----------------------
    
    name,1c17_A
    TMD_seq,AAVMMGLAAIGAAIGIGILG
    full_seq,MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA
    """

    # get the command-line arguments
    args = parser.parse_args()

    sys.stdout.write("args ", args)
    sys.stdout.flush()

    if args.d is not None:
        # process every input file in the args.d input folder
        input_dir = Path(args.d)
        infile_names = glob.glob(os.path.join(input_dir, "*.txt"))
        infile_list = [file for file in infile_names]
    elif args.i is not None:
        # process only a single input file
        infile_list = [Path(args.i)]
        input_dir = Path(args.i).parent
    elif args.d is not None and args.i is not None:
        raise ValueError("Please include either an input directory of files to process (-d D:\data),"
                         "or an input file (-i D:\data\Q12983.txt), but not both.")
    else:
        sys.stdout.write("usage: " + usage + "\n\n")
        raise ValueError("Input error. The command should include an input directory of files to process (-d directory), "
                         "or an input file (-i D:\data\Q12983.txt).")

    # use the output directory given as -f
    # or Q
    if args.f is not None:
        output_dir = Path(args.f)
    else:
        output_dir = Path(os.path.split(input_dir)[0]).joinpath("output")
        if not output_dir.is_dir():
            os.makedirs(output_dir)


    for input_csv in infile_list:
        input_ser = open_csv_as_series(input_csv)

        # convert protein_name to file-format-friendly text, without symbols etc, max 20 characters
        protein_name = slugify(input_ser["name"])[0:20]
        if protein_name != input_ser["name"]:
            sys.stdout.write("\nprotein name modified from {} to directory-folder-friendly {}\n".format(input_ser["name"], protein_name))

        input_ser["slugified_name"] = protein_name

        TMD_seq = input_ser["TMD_seq"]
        full_seq = input_ser["full_seq"]

        #predictions_folder = os.path.normpath(args.f)

        # get checksum
        md5 = get_md5_checksum(TMD_seq, full_seq)
        input_ser["md5"] = md5

        # create output directory based on protein name
        # save the original csv
        out_dir = output_dir.joinpath(md5)
        #out_dir = os.path.join(output_dir, protein_name)
        thoipapy.utils.make_sure_path_exists(out_dir)
        input_ser.to_csv(out_dir / "input.csv")

        run_THOIPA_prediction(protein_name, md5, TMD_seq, full_seq, out_dir)
        os.chmod(out_dir, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP| stat.S_IXGRP)