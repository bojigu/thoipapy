import os
import platform
import re
import sys
import time
from pathlib import Path

import pandas as pd

from thoipapy import utils as utils
from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import SurroundingSequence

# cd-hit rejects any sequence identity threshold below 0.40 with a fatal error, whatever word size
# it is given, so this is the floor of the redundancy-reduction loop and not a tuning choice.
MIN_CDHIT_CUTOFF = 0.40


def rate4site_calculation_mult_prot(paths: ArtefactPaths, df_set, logging):
    """Calculates conservation of positions using rate4site.

    install rate4site for linux: sudo apt-get install rate4site
    install cd-hit for linux: sudo apt-get install cd-hit

    Parameters
    ----------
    paths : ArtefactPaths
        Locations of the pipeline's input and output files.
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    if "Linux" not in platform.system():
        logging.warning("Aborting rate4site calculation, ao is only implemented on linux.")
        return False

    max_n_gaps_in_TMD_subject_seq = paths.max_n_gaps_in_TMD_subject_seq

    for i in df_set.index:
        acc = df_set.at[i, "acc"]
        database = df_set.at[i, "database"]
        TMD_seq = df_set.at[i, "TMD_seq"]
        alignments_dir = paths.alignment_dir(database)

        # input
        fasta_uniq_TMD_seqs_surr5_for_LIPO = (
            alignments_dir / f"{acc}.surr5.gaps{max_n_gaps_in_TMD_subject_seq}.uniq.for_LIPO.fas"
        )

        # output
        rate4site_csv: Path = paths.rate4site_csv(database, acc)

        full_seq_len = len(df_set.at[i, "full_seq"])
        # number of residues surrounding the TMD on the left in alignment
        # (in alignments for lipo and rate4site, this is hard-coded to 5 residues)
        n_residues_surrounding_tmd_in_alignment = 5
        surrounding_sequence_5 = SurroundingSequence(
            df_set.at[i, "TMD_start"], df_set.at[i, "TMD_end"], full_seq_len, n_residues_surrounding_tmd_in_alignment
        )
        surrounding_seq_len_n_term_offset = surrounding_sequence_5.n_term_offset
        rate4site_calculation(
            TMD_seq, acc, fasta_uniq_TMD_seqs_surr5_for_LIPO, rate4site_csv, surrounding_seq_len_n_term_offset, logging
        )

    logging.info("rate4site_calculation finished")


def read_fasta_records(fasta_path: Path | str) -> list[tuple[str, str]]:
    """Read a FASTA file into an ordered list of (header_without_gt, sequence) pairs.

    The alignments thoipapy writes are one header line and one sequence line per record, but a
    record is accumulated across lines anyway so that a wrapped file does not silently lose
    residues.
    """
    records: list[tuple[str, list[str]]] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append((line[1:], []))
            elif records:
                records[-1][1].append(line)
    return [(header, "".join(parts)) for header, parts in records]


def write_cdhit_input(alignment_records: list[tuple[str, str]], cons_cdhit_input_fasta: Path) -> None:
    """Write the alignment for cd-hit, with the gaps stripped from the sequences only.

    This used to run line.replace("-", "") over every line of the file, deflines included. The
    query record is the only one whose defline is not a bare integer, so an accession containing a
    hyphen came back from cd-hit under a different name, failed the membership test that rebuilds
    the alignment, and was dropped. rate4site then scored whichever homologue was left at the top,
    and every residue where that homologue differed from the query was reported under the wrong
    letter -- which combine_all_features refused to merge, two stages later.
    """
    with open(cons_cdhit_input_fasta, "w") as f_out:
        for header, seq in alignment_records:
            f_out.write(f">{header}\n{seq.replace('-', '')}\n")


def truncate_alignment(cons_cdhit_output_fasta: Path, max_n_sequences: int) -> list[str]:
    """Keep the first max_n_sequences records of an alignment, returning the headers kept.

    Reached when cd-hit cannot get below max_n_sequences at the lowest identity threshold it
    accepts. Truncation is by record: the previous version cut at a fixed line count, which left a
    dangling defline with no sequence under it whenever the count was even.
    """
    records = read_fasta_records(cons_cdhit_output_fasta)[:max_n_sequences]
    with open(cons_cdhit_output_fasta, "w") as f_out:
        for header, seq in records:
            f_out.write(f">{header}\n{seq}\n")
    return [header for header, _seq in records]


def rate4site_calculation(
    TMD_seq: str,
    acc: str,
    fasta_uniq_TMD_seqs_surr5_for_LIPO: Path | str,
    rate4site_csv: Path,
    surrounding_seq_len_n_term_offset: int,
    logging,
    rerun_rate4site: bool = False,
):
    """Calculate a conservation score for each residue of the TMD with rate4site.

    rate4site scores positions of a *reference* sequence, which is the query here. Getting the
    query into the alignment as record zero, and keeping it there, is the whole game: everything
    downstream assumes the score at row N belongs to residue N of the query.
    """
    # resolve(): rate4site runs with cwd=output_dir (it writes r4s.res into the working directory
    # whatever -o says), so a relative output path would be resolved a second time against that
    # directory and the run would fail with "unable to open output file for writing".
    output_dir: Path = rate4site_csv.parent.resolve()
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)
    # temp output files
    rate4site_orig_output: Path = output_dir / f"{acc}.rate4site_orig_output.txt"
    cons_cdhit_input_fasta: Path = output_dir / f"{acc}.lipo_seqs_cdhit_input.fas"
    cons_cdhit_output_fasta: Path = output_dir / f"{acc}.lipo_seqs_cdhit_output.fas"
    rate4site_input: Path = output_dir / f"{acc}.rate4site_input.fas"

    alignment_records = read_fasta_records(fasta_uniq_TMD_seqs_surr5_for_LIPO)
    query_header, query_seq = query_record_from_alignment(
        alignment_records, acc, TMD_seq, surrounding_seq_len_n_term_offset, fasta_uniq_TMD_seqs_surr5_for_LIPO
    )

    write_cdhit_input(alignment_records, cons_cdhit_input_fasta)

    # delete output file if it exists
    if cons_cdhit_output_fasta.is_file():
        cons_cdhit_output_fasta.unlink()
    cdhit_cluster_reps: list[str] = []
    if not rate4site_orig_output.is_file() or (rerun_rate4site in [True, 1]):

        sys.stdout.write(f"decreasing cdhit cutoff for {acc}: ")
        sys.stdout.flush()

        # Whole percentage points, not repeated subtraction of 0.01 from 1.0. Sixty subtractions
        # of 0.01 leave 0.39999999999999947, so the loop stopped one usable round early -- at an
        # identity threshold that formatted as 0.41 -- and get_word_size picked a word size one
        # step small at each of its boundaries.
        cutoff_percent = 100
        final_cutoff_used = 1.0
        rerun = False
        len_cdhit_cluster_reps = 1000
        max_n_sequences_for_rate4site = 200

        while len_cdhit_cluster_reps > max_n_sequences_for_rate4site:
            cutoff = cutoff_percent / 100
            if rerun:
                temp: Path = cons_cdhit_output_fasta.with_name(cons_cdhit_output_fasta.stem + ".temp.fas")
                temp.unlink(missing_ok=True)
                os.rename(cons_cdhit_output_fasta, temp)
                cdhit_cluster_reps = run_cdhit(temp, cons_cdhit_output_fasta, cutoff)
                temp.unlink(missing_ok=True)
            else:
                cdhit_cluster_reps = run_cdhit(cons_cdhit_input_fasta, cons_cdhit_output_fasta, cutoff)

            len_cdhit_cluster_reps = len(cdhit_cluster_reps)
            final_cutoff_used = cutoff
            sys.stdout.write(f"cutoff={cutoff:0.2f}(n_reps={len_cdhit_cluster_reps}), ")
            sys.stdout.flush()
            cutoff_percent -= 1
            # cd-hit refuses any threshold below 0.40 outright ("Fatal Error", exit 1), so this
            # loop cannot go lower however many representatives are left. It used to walk down to
            # 0.20, which meant a protein with a deep enough alignment aborted the whole
            # prediction with "cd-hit failed" instead of reaching the truncation below. Deep
            # alignments became normal when the webserver moved from NCBI nr to a local UniRef90.
            if cutoff_percent < round(MIN_CDHIT_CUTOFF * 100):
                cdhit_cluster_reps = truncate_alignment(cons_cdhit_output_fasta, max_n_sequences_for_rate4site)
                logging.warning(
                    f"{acc} still had {len_cdhit_cluster_reps} cd-hit representatives at the "
                    f"lowest cutoff cd-hit accepts ({MIN_CDHIT_CUTOFF:0.2f}). The alignment was "
                    f"truncated to {len(cdhit_cluster_reps)} sequences for rate4site."
                )
                len_cdhit_cluster_reps = len(cdhit_cluster_reps)
                break
            rerun = True

        sys.stdout.write("\n")
        logging.info(
            f"cd-hit for rate4site finished. Final cutoff = {final_cutoff_used:0.2f}. Clusters = {len_cdhit_cluster_reps}. Output = {cons_cdhit_output_fasta}"
        )

        write_rate4site_input(alignment_records, query_header, cdhit_cluster_reps, rate4site_input)

        if not rate4site_orig_output.parent.is_dir():
            rate4site_orig_output.parent.mkdir(parents=True)

        # -a names the reference sequence rate4site scores. Without it rate4site takes the first
        # record, which is what this function guarantees anyway -- but if the guarantee is ever
        # broken, rate4site exits nonzero on an unknown name rather than silently scoring some
        # other sequence. The two together are why the query cannot quietly stop being the
        # reference again.
        exect_str = [
            "rate4site",
            "-s",
            str(rate4site_input),
            "-o",
            str(rate4site_orig_output),
            "-a",
            query_header,
        ]
        logging.info("starting rate4site")
        start = time.time()
        # cwd: rate4site writes r4s.res, r4sOrig.res and TheTree.txt into the working directory
        # whatever -o says. Run it in its own output directory so two predictions running at once
        # cannot overwrite each other's, and so the cleanup below removes the right files.
        command = utils.Command(exect_str, cwd=output_dir)
        command.run(timeout=1200, log_stderr=False)
        if not command.succeeded():
            raise RuntimeError(
                f"rate4site failed: returncode={command.returncode}, "
                f"timed_out={command.timed_out}. Command: {command.command_string}\n{command.stderr}"
            )
        duration = time.time() - start
        logging.info(f"rate4site finished after {duration} seconds")

        if not Path(rate4site_orig_output).is_file():
            raise FileNotFoundError("rate4site output file is not found")

        # cleanup temp files
        temp_output_files = ["r4s.res", "r4sOrig.res", "TheTree.txt"]
        for temp_output_file in temp_output_files:
            if (output_dir / temp_output_file).is_file():
                (output_dir / temp_output_file).unlink()

        if rate4site_orig_output.stat().st_size == 0:
            rate4site_orig_output.unlink()
            raise Exception("rate4site output is empty, file has been deleted, please check input file")

        logging.info(f"{acc} rate4site finished ({rate4site_orig_output})")
    else:
        logging.info(
            f"skipping rate4site algo for existing file {rate4site_orig_output}. Set 'rerun_rate4site' to True to rerun calculation."
        )
    df_rate4site = parse_rate4site_output(rate4site_orig_output, TMD_seq, surrounding_seq_len_n_term_offset, acc)
    df_rate4site.to_csv(rate4site_csv)
    logging.info(f"{rate4site_csv} saved")


def query_record_from_alignment(
    alignment_records: list[tuple[str, str]],
    acc: str,
    TMD_seq: str,
    surrounding_seq_len_n_term_offset: int,
    fasta_path: Path | str,
) -> tuple[str, str]:
    """Return the query record of a homologue alignment, checking it is what rate4site needs.

    NCBI_parser.save_fasta_from_array always writes the query first, as "{acc}_orig_seq". Three
    things about it are load-bearing further down and none of them were ever checked.
    """
    if not alignment_records:
        raise ValueError(f"{acc} rate4site_calculation failed: {fasta_path} contains no sequences.")

    query_header, query_seq = alignment_records[0]
    expected_header = f"{acc}_orig_seq"
    if query_header != expected_header:
        raise ValueError(
            f"{acc} rate4site_calculation failed: the first record of {fasta_path} is "
            f"'{query_header}', not the query '{expected_header}'. rate4site scores positions of "
            f"the reference sequence, so scoring anything but the query gives the wrong residues."
        )
    if query_header.startswith("-"):
        raise ValueError(
            f"{acc} rate4site_calculation failed: the query record is named '{query_header}', "
            f"which rate4site would read as a command-line option rather than a sequence name."
        )
    if "-" in query_seq:
        raise ValueError(
            f"{acc} rate4site_calculation failed: the query record in {fasta_path} contains gaps "
            f"({query_seq}). rate4site numbers its output by ungapped position in the reference, "
            f"so a gapped query would shift every score relative to the TMD."
        )
    offset = surrounding_seq_len_n_term_offset
    tmd_in_query = query_seq[offset : offset + len(TMD_seq)]
    if tmd_in_query != TMD_seq:
        raise ValueError(
            f"{acc} rate4site_calculation failed: the TMD is not at offset {offset} of the query "
            f"record in {fasta_path}.\n  expected: {TMD_seq}\n  found:    {tmd_in_query}\n"
            f"  query record: {query_seq}"
        )
    return query_header, query_seq


def write_rate4site_input(
    alignment_records: list[tuple[str, str]],
    query_header: str,
    cdhit_cluster_reps: list[str],
    rate4site_input: Path,
) -> None:
    """Write the gapped alignment for rate4site: the query first, then the cd-hit representatives.

    The query is written whether or not cd-hit returned it. cd-hit reduces redundancy, which is a
    reason to drop a sequence that has near-identical neighbours -- and the query is exactly such
    a sequence. Losing it changes which sequence rate4site scores.
    """
    reps = set(cdhit_cluster_reps)
    with open(rate4site_input, "w") as f_out:
        for header, seq in alignment_records:
            if header == query_header or header in reps:
                f_out.write(f">{header}\n{seq}\n")


def parse_rate4site_output(
    rate4site_orig_output: Path, TMD_seq: str, surrounding_seq_len_n_term_offset: int, acc: str
) -> pd.DataFrame:
    """Convert rate4site's text output to one conservation score per TMD residue."""
    df = pd.read_csv(
        rate4site_orig_output,
        skiprows=count_rate4site_header_lines(rate4site_orig_output),
        index_col=0,
        header=None,
        sep=r"\s+",
        on_bad_lines="skip",
        comment="#",
    )
    df.columns = ["seq", "score", "qq-interval", "std", "msa-data"]
    df.to_csv(str(rate4site_orig_output)[:-4] + ".orig.csv")
    # convert standard csv to csv for thoipa features
    df_rate4site = df.reindex(["seq", "score"], axis=1)
    df_rate4site.columns = ["residue_name", "rate4site"]
    df_rate4site.index.name = "residue_num"
    # since the lipophilicity alignment was padded by 5 residues on each side,
    # the index doesn't match the TMD residue number.
    TMD_seq_len = len(TMD_seq)
    # minus all indices by the number of residues surrounding the TMD on the left (for lipo/rate4site, this is usually 5)
    df_rate4site.index = df_rate4site.index.to_series() - surrounding_seq_len_n_term_offset
    residue_num_to_keep = range(1, TMD_seq_len + 1)
    df_rate4site = df_rate4site.reindex(residue_num_to_keep)

    # reindex() fills missing positions with NaN. That used to travel two stages downstream and
    # surface as combine_all_features complaining that the merged sequence was the wrong length.
    missing = df_rate4site.index[df_rate4site["rate4site"].isna()].tolist()
    if missing:
        raise ValueError(
            f"{acc} rate4site produced no score for TMD residue(s) {missing} of "
            f"{TMD_seq_len}. Output file: {rate4site_orig_output}"
        )

    reported_seq = "".join(df_rate4site["residue_name"].tolist())
    if reported_seq != TMD_seq:
        differences = [
            f"pos {n}: query {query_residue}, rate4site {reported_residue}"
            for n, (query_residue, reported_residue) in enumerate(zip(TMD_seq, reported_seq), start=1)
            if query_residue != reported_residue
        ]
        raise ValueError(
            f"{acc} rate4site scored a sequence that is not the query TMD.\n"
            f"  query TMD: {TMD_seq}\n  rate4site: {reported_seq}\n  " + "\n  ".join(differences)
        )

    # Take the residue identity from the query rather than from rate4site's echo of its reference.
    # rate4site contributes the score; the residue is already known.
    df_rate4site["residue_name"] = list(TMD_seq)
    return df_rate4site


def count_rate4site_header_lines(rate4site_orig_output: Path) -> int:
    """Number of banner lines before rate4site's first data row.

    Was hardcoded to 13. A banner one line longer or shorter would eat the first data row, or
    admit a comment line, and the resulting error would point at the wrong residue.
    """
    with open(rate4site_orig_output) as f:
        for n, line in enumerate(f):
            if re.match(r"\s*\d+\s", line):
                return n
    raise ValueError(f"rate4site output {rate4site_orig_output} contains no data rows.")


def get_word_size(cutoff: float) -> int:
    """cd-hit's word size for a given identity threshold, per its own documentation."""
    word_size_dict = {0.70: 5, 0.60: 4, 0.50: 3, 0.40: 2}
    for threshold, word_size in word_size_dict.items():
        if cutoff >= threshold:
            return word_size
    # Was `return "error"`, which cd-hit would then have been handed as `-n error`.
    raise ValueError(
        f"cd-hit has no word size for an identity threshold of {cutoff:0.2f}; it rejects "
        f"anything below {MIN_CDHIT_CUTOFF:0.2f}."
    )


def run_cdhit(cons_cdhit_input_fasta: Path, cons_cdhit_output_fasta: Path, cutoff: float) -> list:
    """Cluster an alignment at the given identity threshold, returning the representatives' headers."""
    cdhit_command = [
        "cd-hit",
        "-i",
        str(cons_cdhit_input_fasta),
        "-o",
        str(cons_cdhit_output_fasta),
        "-c",
        f"{cutoff:0.2f}",
    ]
    if cutoff != 1.0:
        cdhit_command += ["-n", str(get_word_size(cutoff))]
    command = utils.Command(cdhit_command)
    command.run(timeout=120, log_stderr=False)
    if not command.succeeded():
        raise RuntimeError(
            f"cd-hit failed: returncode={command.returncode}, timed_out={command.timed_out}. "
            f"Command: {command.command_string}\n{command.stderr}"
        )
    if not cons_cdhit_output_fasta.is_file():
        raise Exception(
            f"cd-hit output not found: input={cons_cdhit_input_fasta}, output={cons_cdhit_output_fasta}, "
            f"commandstring={command.command_string}"
        )
    return [header for header, _seq in read_fasta_records(cons_cdhit_output_fasta)]
