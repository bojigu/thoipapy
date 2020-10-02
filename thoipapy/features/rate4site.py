import os
import platform
import sys
from pathlib import Path
from typing import Union

import pandas as pd

from thoipapy import utils as utils


def rate4site_calculation_mult_prot(s, df_set, logging):
    """Calculates conservation of positions using rate4site.

    install rate4site for linux: sudo apt-get install rate4site
    install cd-hit for linux: sudo apt-get install cd-hit

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if not "Linux" in platform.system():
        logging.warning("Aborting rate4site calculation, ao is only implemented on linux.")
        return False

    max_n_gaps_in_TMD_subject_seq = s["max_n_gaps_in_TMD_subject_seq"]

    for i in df_set.index:
        acc = df_set.at[i, "acc"]
        database = df_set.at[i, "database"]
        TMD_seq = df_set.at[i, "TMD_seq"]
        alignments_dir = Path(s["thoipapy_data_folder"]) / f"homologues/alignments/{database}"

        # input
        fasta_uniq_TMD_seqs_surr5_for_LIPO = alignments_dir / f"{acc}.surr5.gaps{max_n_gaps_in_TMD_subject_seq}.uniq.for_LIPO.fas"

        # output
        rate4site_orig_output: Path = Path(s["thoipapy_data_folder"]).joinpath("features", "rate4site", database, f"{acc}.rate4site_orig_output.txt")
        rate4site_csv: Path = Path(s["thoipapy_data_folder"]).joinpath("features", "rate4site", database, f"{acc}_rate4site.csv")

        rate4site_calculation(TMD_seq, acc, fasta_uniq_TMD_seqs_surr5_for_LIPO, rate4site_csv, logging)

    logging.info("rate4site_calculation finished")


def rate4site_calculation(TMD_seq: str, acc: str, fasta_uniq_TMD_seqs_surr5_for_LIPO: Union[Path, str], rate4site_csv: Path, logging, rerun_rate4site: bool = False):
    output_dir: Union[Path, str] = rate4site_csv.parent
    assert output_dir.is_dir()
    # temp output files
    rate4site_orig_output: Union[Path, str] = output_dir/ f"{acc}.rate4site_orig_output.txt"
    cons_cdhit_input_fasta: Union[Path, str] = output_dir/ f"{acc}.lipo_seqs_cdhit_input.fas"
    cons_cdhit_output_fasta: Union[Path, str] = output_dir/  f"{acc}.lipo_seqs_cdhit_output.fas"
    rate4site_input: Union[Path, str] = output_dir/ f"{acc}.rate4site_input.fas"
    with open(cons_cdhit_input_fasta, "w") as f_out:
        with open(fasta_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
            for line in f_in:
                f_out.write(line.replace("-", ""))
    # delete output file if it exists
    if cons_cdhit_output_fasta.is_file():
        cons_cdhit_output_fasta.unlink()
    len_cdhit_cluster_reps = 1000
    max_n_sequences_for_rate4site = 200
    cutoff = 1.0
    cutoff_decrease_per_round = 0.01
    rerun = False
    if not rate4site_orig_output.is_file() or (rerun_rate4site in [True, 1]):

        sys.stdout.write(f"\ndecreasing cdhit cutoff for {acc}: ")
        sys.stdout.flush()

        while len_cdhit_cluster_reps > max_n_sequences_for_rate4site:
            if rerun:
                temp = str(cons_cdhit_output_fasta)[:-4] + "temp.fas"
                if Path(temp).is_file():
                    Path(temp).unlink()
                os.rename(cons_cdhit_output_fasta, temp)
                cdhit_cluster_reps = run_cdhit(temp, cons_cdhit_output_fasta, cutoff)
            else:
                cdhit_cluster_reps = run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)

            len_cdhit_cluster_reps = len(cdhit_cluster_reps)
            sys.stdout.write(f"{cutoff:0.2f}({len_cdhit_cluster_reps}), ")
            sys.stdout.flush()
            cutoff -= cutoff_decrease_per_round
            if cutoff <= 0.20:
                to_be_truncated_fasta: str = str(cons_cdhit_output_fasta)[:-4] + "to_be_truncated.fas"
                os.rename(cons_cdhit_output_fasta, to_be_truncated_fasta)
                with open(to_be_truncated_fasta, "r") as f_in:
                    with open(cons_cdhit_output_fasta, "w") as f_out:
                        for n, line in enumerate(f_in):
                            f_out.write(line)
                            if n >= 400:
                                break
                break
            rerun = True

        final_cutoff_used = cutoff + cutoff_decrease_per_round
        sys.stdout.write("\n")
        logging.info(f"cd-hit for rate4site finished. Final cutoff = {final_cutoff_used:0.2f}. Clusters = {len_cdhit_cluster_reps}. Output = {cons_cdhit_output_fasta}")

        # if len(cdhit_cluster_reps) > 200:
        #    cutoff = 0.7
        #    run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)
        #    logging.info(f"repeating cd-hit with a stricter cutoff\nafter cd-hit analysis, {cons_cdhit_output_fasta} has {len(cdhit_cluster_reps)} clusters")

        copy_sequence = False
        with open(fasta_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
            with open(str(rate4site_input), "w") as f_out:
                for line in f_in:
                    if line[0] == ">":
                        if line[1:] in cdhit_cluster_reps:
                            f_out.write(line)
                            copy_sequence = True
                    else:
                        if copy_sequence:
                            f_out.write(line)
                            copy_sequence = False

        if os.path.isfile(rate4site_input):

            if not rate4site_orig_output.parent.is_dir():
                rate4site_orig_output.parent.mkdir(parents=True)

            exect_str = f"rate4site -s {rate4site_input} -o {rate4site_orig_output}"
            sys.stdout.write(exect_str)
            sys.stdout.flush()
            command = utils.Command(exect_str)
            command.run(timeout=1200, log_stderr=False)

            if not Path(rate4site_orig_output).is_file():
                raise FileNotFoundError("rate4site output file is not found")

            # cleanup temp files
            temp_output_files = ["r4s.res", "r4sOrig.res", "TheTree.txt"]
            for temp_output_file in temp_output_files:
                if Path(temp_output_file).is_file():
                    Path(temp_output_file).unlink()

            if rate4site_orig_output.stat().st_size == 0:
                rate4site_orig_output.unlink()
                raise Exception("rate4site output is empty, file has been deleted, please check input file")

            logging.info('{} rate4site finished ({})'.format(acc, rate4site_orig_output))

        else:
            logging.warning("{} rate4site failed. {} input file not found".format(acc, fasta_uniq_TMD_seqs_surr5_for_LIPO))
    else:
        logging.info(f"skipping rate4site algo for existing file {rate4site_orig_output}. Set 'rerun_rate4site' to True to rerun calculation.")
    # convert text output to standard csv
    df = pd.read_csv(rate4site_orig_output, skiprows=range(13), index_col=0, header=None, delim_whitespace=True, error_bad_lines=False, comment="#")
    df.columns = ["seq", "score", "qq-interval", "std", "msa-data"]
    df.to_csv(str(rate4site_orig_output)[:-4] + ".orig.csv")
    # convert standard csv to csv for thoipa features
    df_rate4site = df.reindex(["seq", "score"], axis=1)
    df_rate4site.columns = ["residue_name", "rate4site"]
    df_rate4site.index.name = "residue_num"
    # since the lipophilicity alignment was padded by 5 residues on each side, the index doesn't match the TMD residue number.
    # first get the offset
    alignment_TMD_plus_offset = "".join(df_rate4site["residue_name"].to_list())
    TMD_seq_len = len(TMD_seq)
    for offset in range(0, (len(alignment_TMD_plus_offset) - TMD_seq_len)):
        TMD_seq_slice = alignment_TMD_plus_offset[offset: offset + TMD_seq_len]
        if TMD_seq_slice == TMD_seq:
            break
    logging.info(f"alignment offset found ({offset})")
    # minus all indices by the number of residues surrounding the TMD on the left
    df_rate4site.index = pd.Series(df_rate4site.index) - offset
    residue_num_to_keep = range(1, TMD_seq_len + 1)
    df_rate4site = df_rate4site.reindex(residue_num_to_keep)
    saved_TMD_seq = "".join(df_rate4site["residue_name"].to_list())
    assert saved_TMD_seq == TMD_seq
    df_rate4site.to_csv(rate4site_csv)
    logging.info(f"{rate4site_csv} saved")


def get_word_size(cutoff):
    word_size_dict = {
        0.70: 5,
        0.60: 4,
        0.50: 3,
        0.40: 2,
        0.00: 1
    }
    for threshold, word_size in word_size_dict.items():
        if cutoff >= threshold:
            return word_size

    return "error"


def run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff):
    word_size = get_word_size(cutoff)
    word_size_command = "" if cutoff == 1.0 else f"-n {word_size}"
    exect_str = f"cdhit -i {cons_cdhit_input_fasta} -o {cons_cdhit_output_fasta} -c {cutoff:0.2f} {word_size_command}"
    command = utils.Command(exect_str)
    command.run(timeout=120, log_stderr=False)
    assert cons_cdhit_output_fasta.is_file()
    cdhit_cluster_reps = []
    with open(str(cons_cdhit_output_fasta), "r") as f:
        for line in f:
            if line[0] == ">":
                cdhit_cluster_reps.append(line[1:])
    return cdhit_cluster_reps