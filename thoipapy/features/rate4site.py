import os
import platform
import sys
from pathlib import Path

import pandas as pd

from thoipapy import utils as utils


def rate4site_calculation(s, df_set, logging):
    """Calculates conservation of positions using rate4site.

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

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        alignments_dir = os.path.join(s["thoipapy_data_folder"], "homologues", "alignments", database)
        path_uniq_TMD_seqs_surr5_for_LIPO = os.path.join(alignments_dir, "{}.surr5.gaps{}.uniq.for_LIPO.fas".format(acc, s["max_n_gaps_in_TMD_subject_seq"]))
        cons_cdhit_input_fasta: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.lipo_seqs_cdhit_input.fas")
        cons_cdhit_output_fasta: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.lipo_seqs_cdhit_output.fas")

        rate4site_input: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.rate4site_input.fas")
        rate4site_orig_output: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}.rate4site_orig_output.txt")

        rate4site_csv: Path = Path(s["thoipapy_data_folder"]).joinpath("Features", "rate4site", database, f"{acc}_rate4site.txt")

        with open(str(cons_cdhit_input_fasta), "w") as f_out:
            with open(path_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
                for line in f_in:
                    f_out.write(line.replace("-",""))

        # delete output file if it exists
        if cons_cdhit_output_fasta.is_file():
            cons_cdhit_output_fasta.unlink()

        len_cdhit_cluster_reps = 1000

        cutoff = 1.0
        cutoff_decrease_per_round = 0.01
        rerun = False

        print(f"decreasing cdhit cutoff for {acc}: ")

        while len_cdhit_cluster_reps > 200:
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

        #if len(cdhit_cluster_reps) > 200:
        #    cutoff = 0.7
        #    run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)
        #    logging.info(f"repeating cd-hit with a stricter cutoff\nafter cd-hit analysis, {cons_cdhit_output_fasta} has {len(cdhit_cluster_reps)} clusters")


        copy_sequence = False
        with open(path_uniq_TMD_seqs_surr5_for_LIPO, "r") as f_in:
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


        if not rate4site_orig_output.is_file() or (s["rerun_rate4site"] in [True, 1]):

            if os.path.isfile(rate4site_input):

                if not rate4site_orig_output.parent.is_dir():
                    rate4site_orig_output.parent.mkdir(parents=True)

                exect_str = f"rate4site -s {rate4site_input} -o {rate4site_orig_output}"
                print(exect_str)
                command = utils.Command(exect_str)
                command.run(timeout=1200, log_stderr=False)

                if not Path(rate4site_orig_output).is_file():
                    raise FileNotFoundError("rate4site output file is not found")

                # cleanup temp files
                temp_output_files = ["r4s.res", "r4sOrig.res", "TheTree.txt"]
                for temp_output_file in temp_output_files:
                    if Path(temp_output_file).is_file():
                        Path(temp_output_file).unlink()

                logging.info('{} rate4site finished ({})'.format(acc, rate4site_orig_output))

            else:
                logging.warning("{} rate4site failed. {} input file not found".format(acc, path_uniq_TMD_seqs_surr5_for_LIPO))
        else:
            logging.info(f"skipping rate4site algo for existing file {rate4site_orig_output}. Set 'rerun_rate4site' to True to rerun calculation.")

        # convert text output to standard csv
        df = pd.read_csv(rate4site_orig_output, skiprows=range(13), index_col=0, header=None, delim_whitespace=True, error_bad_lines=False, comment="#")
        df.columns = ["seq", "score", "qq-interval", "std", "msa-data"]
        df.to_csv(str(rate4site_orig_output)[:-4] + ".orig.csv")

        # convert standard csv to csv for thoipa features
        df_cons = df.reindex(["seq", "score"], axis=1)
        df_cons.columns = ["residue_name", "rate4site"]
        df_cons.index.name = "residue_num"

        df_cons.to_csv(rate4site_csv)

        logging.info(f"{rate4site_csv} saved")

    logging.info("rate4site_calculation finished")


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