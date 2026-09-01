#!/usr/bin/env python
"""
Utilities file containing useful functions.
More recent functions are at the top.
Many of these are taken from he korbinian python package by Mark Teese. This is allowed under the permissive MIT license.
"""

import ctypes
import logging
import os
import platform
import re as re
import shlex
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

import matplotlib.colors as colors
import numpy as np
import pandas as pd

# Above this fraction of residues lost to NaN, a dropna() is treated as a pipeline failure rather
# than as data cleaning. See dropna_with_report.
MAX_FRACTION_ROWS_DROPPED = 0.05


class Command:
    """Run an external command-line program and record how it exited.

    The command is a list of arguments, run without a shell. It used to be a single string run
    with ``shell=True``, which meant every character of an interpolated path reached ``/bin/sh``.
    Since rate4site's and cd-hit's filenames are built from the protein name, and the protein name
    is whatever the user put in their FASTA header, an ordinary UniProt header ("sp|Q8WWF3|...")
    was split into shell commands, and anything else could be injected the same way.

    Every caller runs an external binary (freecontact, rate4site, cd-hit). The exit status used to
    be discarded, and a shell redirect creates the output file before the binary is resolved -- so
    a missing tool left a 0-byte file behind and looked exactly like success. Callers must check
    succeeded() rather than the mere existence of an output file.

    Parameters
    ----------
    cmd : sequence of str
        Program and arguments. A str is rejected: without a shell it would be taken as one
        enormous program name, and the failure would be reported as a missing file.
    stdin_bytes : bytes, optional
        Fed to the program on stdin. Replaces the shell pipelines that used to feed it.
    stdout_path : Path, optional
        Where the program's stdout is written. Replaces the shell "> file" redirects.
    cwd : Path, optional
        Working directory for the program. rate4site writes r4s.res, r4sOrig.res and TheTree.txt
        into the current directory whatever the -o argument says, so two predictions running in
        one directory overwrite each other's temporary files.
    """

    def __init__(
        self,
        cmd: Sequence[str],
        stdin_bytes: bytes | None = None,
        stdout_path: Path | str | None = None,
        cwd: Path | str | None = None,
    ):
        if isinstance(cmd, str):
            raise TypeError(
                "Command takes a list of arguments and does not use a shell. "
                f"Pass ['prog', '-i', str(path)] rather than the string {cmd!r}."
            )
        self.cmd: list[str] = [str(argument) for argument in cmd]
        self.stdin_bytes: bytes | None = stdin_bytes
        self.stdout_path: Path | None = Path(stdout_path) if stdout_path is not None else None
        self.cwd: Path | None = Path(cwd) if cwd is not None else None
        self.returncode: int | None = None
        self.timed_out: bool = False
        self.stderr: str = ""

    @property
    def command_string(self) -> str:
        """The command as a copy-pasteable shell line, for log and exception messages only."""
        return shlex.join(self.cmd)

    def run(self, timeout: float, log_stderr: bool = True) -> int | None:
        """Run the command, recording its exit status in self.returncode and self.timed_out."""
        stdout_file = open(self.stdout_path, "wb") if self.stdout_path is not None else None
        try:
            completed = subprocess.run(
                self.cmd,
                input=self.stdin_bytes,
                stdout=stdout_file if stdout_file is not None else subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                cwd=self.cwd,
                check=False,
            )
            self.returncode = completed.returncode
            self.stderr = completed.stderr.decode("utf-8", errors="replace")
        except subprocess.TimeoutExpired as e:
            # subprocess.run kills the child before re-raising, which the previous implementation
            # did not: it called terminate() and then joined the reader thread with no timeout, so
            # a child that ignored SIGTERM hung the pipeline forever.
            self.timed_out = True
            self.returncode = None
            self.stderr = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
            logging.warning(f"Terminating process after {timeout}s timeout: {self.command_string}")
        except OSError as e:
            # Without a shell there is nothing to report "command not found" as an exit status:
            # the missing binary raises here instead. 127 is what the shell used to return, and a
            # missing external tool is the most common way this fails on a fresh container, so it
            # is reported through the same channel as any other failure rather than as a traceback
            # from inside a utility class.
            self.returncode = 127
            self.stderr = f"{type(e).__name__}: {e}"
        finally:
            if stdout_file is not None:
                stdout_file.close()

        # if the console prints anything longer than 5 characters, log it
        if len(self.stderr) > 5 and log_stderr:
            logging.warning(f"FAULTS: {self.stderr}")

        # A nonzero exit is reported even when log_stderr is False. log_stderr suppresses the
        # noisy stderr *content* of tools that chatter on success; it must not suppress failure.
        if not self.succeeded():
            logging.warning(
                f"Command failed (returncode={self.returncode}, timed_out={self.timed_out}): "
                f"{self.command_string}\n{self.stderr}"
            )

        return self.returncode

    def succeeded(self) -> bool:
        """True only if the command ran to completion with a zero exit status."""
        return self.returncode == 0 and not self.timed_out


def residuals(constants, function, x, y):
    """Cost function (Kostenfunktion) used to optimise the fit of a curve to data.

    Returns the distance between the y-values of the real data and the y-values produced by the
    fitted function (sigmoid, sine, etc.). This is the quantity that scipy.optimize.leastsq
    minimises.

    Copied from pytoxr.mathfunctions (https://github.com/teese/pytoxr) by Mark Teese, so that
    thoipapy does not take a dependency on pytoxr for two short functions. This is allowed under
    the permissive MIT license.

    Parameters
    ----------
    constants : arraylike
        The parameters of `function` that are being optimised.
    function : callable
        A function of the form f(constants, x), whose fit to the data is being optimised.
    x : np.ndarray
        x-values of the real data.
    y : arraylike
        y-values of the real data.

    Returns
    -------
    residuals : np.ndarray
        Element-wise difference between the real and the fitted y-values.
    """
    return y - function(constants, x)


def sine_perfect_helix(sine_constants_cd, x):
    """Sine equation constrained to alpha-helical periodicity, 3.6 residues per turn.

    f(x) = a * sin(b * x + c) + d, where a and b are fixed and only c and d are fitted.

    Why is a fixed to 0.2?
        This is arbitrary, resulting in a curve that is 0.4 in height.
    Why is b fixed to 1.745?
        Since periodicity = 2 * np.pi / b, for a perfect alpha helix b = 2 * np.pi / 3.6 = 1.745.

    Copied from pytoxr.mathfunctions (https://github.com/teese/pytoxr) by Mark Teese, so that
    thoipapy does not take a dependency on pytoxr for two short functions. This is allowed under
    the permissive MIT license.

    Parameters
    ----------
    sine_constants_cd : arraylike
        The two fitted constants: c, the phase shift, and d, the vertical offset.
    x : np.ndarray
        x-values, in residue numbers.

    Returns
    -------
    y : np.ndarray
        y-values of the constrained sine curve at each position in `x`.
    """
    a = 0.2
    b = 1.745
    c, d = sine_constants_cd
    y = a * np.sin(b * x + c) + d
    return y


def make_sure_path_exists(input_path: Path | str, isfile: bool = False):
    """If path to directory or folder doesn't exist, creates the necessary directory.
    Set isfile=True to indicate a filepath, where the parent directory needs to be created.
    """
    input_path = Path(input_path)

    if isfile:
        directory = input_path.parent
    else:
        directory = input_path

    if directory.exists():
        return
    else:
        directory.mkdir(parents=True)


def reorder_dataframe_columns(dataframe, cols, front=True):
    """Takes a dataframe and a subsequence of its columns,
       returns dataframe with seq as first columns if "front" is True,
       and seq as last columns if "front" is False.
       taken from https://stackoverflow.com/questions/12329853/how-to-rearrange-pandas-column-sequence

    Parameters
    ----------
    df : pd.DataFrame
        original dataframe
    cols : list
        List of cols to put first or last.
    front : bool
        Whether to put the columns at the front or the back.

    Usage
    -----
    df = reorder_dataframe_columns(df, ["TMD_start", "database", "whatever other col I want first"])
    """
    for col in cols:
        if col not in dataframe.columns:
            cols.remove(col)
    cols_to_place_first = cols[:]  # copy so we don't mutate seq
    for existing_col in dataframe.columns:
        if existing_col not in cols_to_place_first:
            if front:  # we want "seq" to be in the front
                # so append current column to the end of the list
                cols_to_place_first.append(existing_col)
            else:
                # we want "seq" to be last, so insert this
                # column in the front of the new column list
                # "cols" we are building:
                cols_to_place_first.insert(0, existing_col)
    return dataframe[cols_to_place_first]


def delete_BLAST_xml(blast_xml_file):
    """Small function to remove files that are already compressed in a tarball.

    Also deletes the associated text file with the date.

    Parameters
    ----------
    blast_xml_file : str
        Path to BLAST xml file
    """
    xml_txt = str(blast_xml_file)[:-4] + "_details.txt"

    # delete the original files
    try:
        os.remove(blast_xml_file)
    except OSError:
        sys.stdout.write(f"{blast_xml_file} could not be deleted")
    try:
        os.remove(xml_txt)
    except OSError:
        sys.stdout.write(f"{xml_txt} could not be deleted")


def get_n_of_gaps_at_start_and_end_of_seq(seq):
    start = 0
    end = 0
    for aa in seq:
        if aa == "-":
            start += 1
        else:
            break
    for aa in seq[::-1]:
        if aa == "-":
            end += 1
        else:
            break
    return start, end


def get_list_residues_in_motif(seq, motif_ss, motif_len):
    # list of positions that matched the start of a motif
    match_start_list = []
    match_end_list = [0] * motif_len
    # counter for the matches

    for start in range(len(seq) - motif_len):
        # for a SmallxxxSmall motif, the end is 4 residues later
        end = start + motif_len
        # get the matched segment
        segment = seq[start : end + 1]
        # check if the segment contains a motif
        match = re.match(motif_ss, segment)
        if match:
            # classify position as start of a motif
            match_start_list.append(1)
            match_end_list.append(1)
        else:
            match_start_list.append(0)
            match_end_list.append(0)
    # add the final zeros on the end of the list, so it's length matches the original sequence
    match_start_list = match_start_list + [0] * motif_len

    match_start_arr = np.array(match_start_list)
    match_end_arr = np.array(match_end_list)

    list_residues_in_motif = match_start_arr + match_end_arr
    list_residues_in_motif[list_residues_in_motif > 1] = 1

    # sys.stdout.write(seq, "seq")
    # sys.stdout.write("".join([str(x) for x in match_start_list]), "match_start_list")
    # sys.stdout.write("".join([str(x) for x in match_end_list]), "match_end_list")
    # sys.stdout.write("".join([str(x) for x in list_residues_in_motif]), "list_residues_in_motif")

    return list_residues_in_motif


def slice_TMD_seq_pl_surr(df_set):
    # note that due to uniprot-like indexing, the start index = start-1
    return df_set["full_seq"][int(df_set["TMD_start_pl_surr"] - 1) : int(df_set["TMD_end_pl_surr"])]


def create_column_with_TMD_plus_surround_seq(df_set, num_of_sur_residues):
    df_set["TMD_start_pl_surr"] = df_set.TMD_start - num_of_sur_residues
    df_set.loc[df_set["TMD_start_pl_surr"] < 1, "TMD_start_pl_surr"] = 1
    df_set["TMD_end_pl_surr"] = df_set.TMD_end + num_of_sur_residues
    for i in df_set.index:
        # acc = df_set.loc[i, "acc"]
        if df_set.loc[i, "TMD_end_pl_surr"] > df_set.loc[i, "seqlen"]:
            df_set.loc[i, "TMD_end_pl_surr"] = df_set.loc[i, "seqlen"]
    TMD_seq_pl_surr_series = df_set.apply(slice_TMD_seq_pl_surr, axis=1)
    return df_set, TMD_seq_pl_surr_series


def create_namedict(names_csv_path, style="shortname [acc-db]"):
    """Create protein name dictionary from the CSV of detailed protein info.

    e.g. namedict[P02724-NMR

    Parameters
    ----------
    names_csv_path : str
        Path to the CSV of manually edited protein names
    style : str
        Style of protein name in output dictionary. Current options are "shortname [acc-db]" or "shortname [acc]"
        "shortname [acc-db]" = 'Q9Y286-ETRA': 'Siglec7 [Q9Y286-ETRA]'
        "shortname [acc]" = 'Q9Y286-ETRA': 'Siglec7 [Q9Y286]'

    Returns
    -------
    namedict : dict
        Dictionary in format namedict[acc_db] = "formatted protein name"
    """
    #################################################################
    #                EXTRACT NAMES FROM THE NAMES CSV               #
    #################################################################
    df_names = pd.read_csv(names_csv_path, index_col=0)
    # restrict names dict to only that database
    df_names["acc"] = df_names.index
    df_names["acc_db"] = df_names.acc + "-" + df_names.database
    df_names.set_index("acc_db", inplace=True, drop=False)
    df_names.index.name = "acc_db_index"

    # df_names.acc_db_for_figs = df_names.acc_db.replace("crystal", "X-ray")

    # add old names in index "e.g. Q13563-crystal", so that they are replaced with new "X-ray" names in figs
    xray_row_bool_ser = df_names.acc_db.str.contains("X-ray")
    df_xray = df_names.loc[xray_row_bool_ser == True].copy()  # noqa: E712  bool mask, not truthiness
    df_xray.index = df_xray["PDB acc"] + "-crystal"
    df_xray["acc_db"] = df_xray["PDB acc"] + "-" + df_xray.database
    df_names = pd.concat([df_names.loc[xray_row_bool_ser == False], df_xray])  # noqa: E712  bool mask

    # df_names = df_names.loc[df_names.database == database]
    if style == "shortname [acc-db]":
        df_names["label"] = df_names.shortname + " [" + df_names.acc_db + "]"
    elif style == "shortname [acc]":
        df_names["label"] = df_names.shortname + " [" + df_names.acc + "]"
    else:
        raise ValueError("other styles not implemented")

    namedict = df_names["label"].to_dict()
    return namedict


def pdf_subpath(png_path):
    """Create a subfolder "pdf" where the png is saved.

    Also checks that the path exists.

    Parameters
    ----------
    png_path : str
        Path to png file

    Returns
    -------

    """
    pdf_path = os.path.join(os.path.dirname(png_path), "pdf", os.path.basename(png_path)[:-4] + ".pdf")
    make_sure_path_exists(pdf_path, isfile=True)
    return pdf_path


def drop_redundant_proteins_from_list(df_set, logging):
    """Simply drops the proteins labeled "False" as CD-HIT cluster representatives.

    Relies on the dataframe containing the cdhit_cluster_rep column.

    Parameters
    ----------
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Returns
    -------
    df_set_nonred : pd.DataFrame
        df_set with only non-redundant proteins, based on redundancy settings
    """
    # `"x" in series` tests the INDEX, not the values, so this warning never fired: the index is
    # a RangeIndex. Use .eq().any() to actually inspect the column.
    no_cdhit_results = df_set.cdhit_cluster_rep.eq("no_cdhit_results").any()
    if no_cdhit_results:
        logging.warning(
            "No CD-HIT results were found, so no redundancy reduction was applied. All proteins "
            "are being used to train the model. To apply redundancy reduction, run CD-HIT and "
            "ensure the .clstr.sorted.txt output is present in results/<setname>/clusters/."
        )

    n_prot_initial = df_set.shape[0]
    # create model only using CD-HIT cluster representatives.
    # `!= False` was comparing against the string "no_cdhit_results", which is never equal to
    # False, so nothing was ever dropped even when CD-HIT results were absent. Compare against
    # the boolean explicitly, and only when there are real results to compare.
    if no_cdhit_results:
        df_set_nonred = df_set
    else:
        df_set_nonred = df_set.loc[df_set.cdhit_cluster_rep.astype(bool)]
    n_prot_final = df_set_nonred.shape[0]
    logging.info(
        f"CDHIT redundancy reduction : n_prot_initial = {n_prot_initial}, n_prot_final = {n_prot_final}, n_prot_dropped = {n_prot_initial - n_prot_final}"
    )
    return df_set_nonred


def add_res_num_full_seq_to_df(acc, df, TMD_seq, full_seq, prediction_name, file):
    """

    Parameters
    ----------
    df : pd.DataFrame
    TMD_seq : str
        TMD sequence according to df_set
    full_seq : str
        Full protein sequence according to df_set
    Returns
    -------

    """
    # re.escape: TMD_seq is sequence data, not a pattern. Same defect as in predict.py and
    # common.py -- an unescaped "*" or "[" raises re.error instead of reporting a bad input.
    m = re.search(re.escape(TMD_seq), full_seq)
    if m:
        # convert from python indexing to unprot indexing
        TMD_start = m.start() + 1
        m.end()
        df["res_num_full_seq"] = np.array(range(df.shape[0])) + TMD_start
    else:
        raise IndexError(
            f"TMD seq not found in full_seq.\nacc = {acc}\nTMD_seq = {TMD_seq}\nfull_seq = {full_seq}\n"
            f"prediction_name={prediction_name},file={file}"
        )
    return df


def shorten(x):
    """
    convert 3-letter amino acid name to 1-letter form
    Parameters
    ----------
    x : str , thrre letter aa sequence

    Returns
    y : str, one-letter aa sequence
    -------

    """
    d = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        "UNK": "U",
    }
    if len(x) % 3 != 0:
        raise ValueError("Input length should be a multiple of three")

    y = ""
    for i in range(int(len(x) / 3)):
        y += d[x[3 * i : 3 * i + 3]]
    return y


def add_mutation_missed_residues_with_na(combined_features_csv, acc, database, df):
    acc_combind_feature_file = combined_features_csv
    df_feature = pd.read_csv(acc_combind_feature_file, engine="python", index_col=0)
    not_df_index = [element for element in df_feature["residue_num"].values if element not in df.index.values]
    for element in not_df_index:
        df.loc[element] = [df_feature.loc[element - 1, "residue_name"], np.nan]
    df = df.sort_index()
    return df


def normalise_0_1(arraylike):
    """Normalise an array to values between 0 and 1.

    The following linear formula is used.
    norm_array = (orig_array - array_min)/(array_max - array_min)

    The use of this simple linear formula allows the normalised data to be "denormalised" later, so long as
    the min and max values of the original array are known.

    Parameters
    ----------
    arraylike : array
        Numpy array (or other arraylike) dataset of floats or ints to be normalised.

    Returns
    -------
    normalised : array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # or, if denormalisation is not necessary
    normalised_array = normalise_0_1(original_array)[0]
    # for further usage examples, see the docstring for denormalise_0_1
    """
    array_min = np.min(arraylike)
    array_max = np.max(arraylike)
    normalised = (arraylike - array_min) / (array_max - array_min)
    # convert to float
    normalised = np.array(normalised).astype(float)
    return normalised, array_min, array_max


def denormalise_0_1(value_or_array, array_min, array_max):
    """Denormalise a value or array to orig values.

    For use after normalisation between 0 and 1 with the normalise_0_1 function.

    The normalisation formula (normalise_0_1):
        norm_array = (orig_array - array_min)/(array_max - array_min)

    The denormalisation formula (denormalise_0_1):
        denormalised_array = norm_array*(array_max - array_min) + array_min

    Parameters
    ----------
    value_or_array : int, float or arraylike
        Int or float to be denormalised.
        Numpy array (or other arraylike) of data (float, int, etc) to be denormalised.

    Returns
    -------
    normalised : float, or numpy array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    from thoipapy.utils import normalise_0_1, denormalise_0_1
    import numpy as np
    original_array = np.linspace(10,130,10)
    original_array[2], original_array[4] = 3, 140
    # normalise original array
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # do stuff to normalised array (e.g., multiply by 0.5)
    normalised_array_halved = normalised_array * 0.5
    # denormalise values to match the equivalents in the original array.
    # Note that the min value (3) was normalised to zero, and was therefore not affected by multiplication.
    normalised_array_halved_denorm = denormalise_0_1(normalised_array_halved, min_, max_)
    # now calculate average values, and check that they match
    norm_array_mean = np.mean(normalised_array)
    norm_array_mean_denormalised = denormalise_0_1(norm_array_mean, min_, max_)
    orig_array_mean = np.mean(original_array)
    # print the two mean values. They should be equal.
    norm_array_mean_denormalised == orig_array_mean
    """
    if isinstance(value_or_array, list):
        raise ValueError(
            "this function accepts arraylike data, not a list. " "Please check data or convert list to numpy array"
        )
    elif isinstance(value_or_array, float):
        denormalised = value_or_array * (array_max - array_min) + array_min
    elif isinstance(value_or_array, np.ndarray):
        denormalised = value_or_array * (array_max - array_min) + array_min
    elif isinstance(value_or_array, pd.Series):
        denormalised = value_or_array * (array_max - array_min) + array_min
    else:
        sys.stdout.write(
            "Unknown datatype. denormalise_0_1 has been given an input that does not appear to be "
            "an int, float, np.ndarray or pandas Series\n"
            "Attempting to process as if it is arraylike....."
        )
    return denormalised


def normalise_between_2_values(arraylike, min_value, max_value, invert=False):
    """Normalises an array of data between two desired values.

    Any values below min_value will be converted to 0.
    Any values above max_value will be converted to 1.
    Optionally, the normalised array can be inverted, so that the original highest
    values are 0, and the original lowest values are now 1.

    Parameters
    ----------
    arraylike : np.ndarray
        Arraylike original data (numpy array or pandas Series)
    min_value : float
        Desired minimum value for normalisation
    max_value : float
        Desired max value for normalisation
    invert : bool
        If True, normalised data will be inverted (former highest value = 0)

    Returns
    -------
    normalised : np.ndarray
        Normalised array of data

    Usage
    -----
    from thoipapy.utils import normalise_between_2_values
    # for array
    orig_array = np.array(range(0, 15))
    norm_array = normalise_between_2_values(orig_array, 3, 10)
    # for pandas Dataframe
    df["norm_data"] = normalise_between_2_values(df["orig_data"], 3, 10)
    """
    # normalise array between min and max values
    normalised = (arraylike - min_value) / (max_value - min_value)
    # replace anything above 1 with 1
    normalised[normalised > 1] = 1
    # replace anything below 0 with 0
    normalised[normalised < 0] = 0
    # if desired, invert the normalised values
    if invert:
        normalised = abs(normalised - 1)
    return normalised


def create_colour_lists():
    """
    Converts several lists of rgb colours to the python format (normalized to between 0 and 1)
    Returns a dictionary that contains dictionaries of palettes with named colours (eg. TUM blues)
    and also lists of unnamed colours (e.g. tableau20)
    (copied from tlabtools 2016.08.08)
    """
    output_dict = {}

    matplotlib_150 = list(colors.cnames.values())
    output_dict["matplotlib_150"] = matplotlib_150

    # define colour dictionaries. TUM colours are based on the style guide.
    colour_dicts = {
        "TUM_colours": {
            "TUMBlue": (34, 99, 169),
            "TUM1": (100, 160, 200),
            "TUM2": (1, 51, 89),
            "TUM3": (42, 110, 177),
            "TUM4": (153, 198, 231),
            "TUM5": (0, 82, 147),
        },
        "TUM_oranges": {
            "TUM0": (202, 101, 10),
            "TUM1": (213, 148, 96),
            "TUM2": (102, 49, 5),
            "TUM3": (220, 108, 11),
            "TUM4": (247, 194, 148),
            "TUM5": (160, 78, 8),
        },
        "TUM_accents": {
            "green": (162, 183, 0),
            "orange": (227, 114, 34),
            "ivory": (218, 215, 203),
        },
    }

    # convert the nested dicts to python 0 to 1 format
    for c_dict in colour_dicts:
        for c in colour_dicts[c_dict]:
            # define r, g, b as ints
            r, g, b = colour_dicts[c_dict][c]
            # normalise r, g, b and add to dict
            colour_dicts[c_dict][c] = (r / 255.0, g / 255.0, b / 255.0)
        # add normalised colours to output dictionary
        output_dict[c_dict] = colour_dicts[c_dict]

    # define colour lists
    colour_lists = {
        "tableau20": [
            (31, 119, 180),
            (174, 199, 232),
            (255, 127, 14),
            (255, 187, 120),
            (44, 160, 44),
            (152, 223, 138),
            (214, 39, 40),
            (255, 152, 150),
            (148, 103, 189),
            (197, 176, 213),
            (140, 86, 75),
            (196, 156, 148),
            (227, 119, 194),
            (247, 182, 210),
            (127, 127, 127),
            (199, 199, 199),
            (188, 189, 34),
            (219, 219, 141),
            (23, 190, 207),
            (158, 218, 229),
        ],
        "tableau20blind": [
            (0, 107, 164),
            (255, 128, 14),
            (171, 171, 171),
            (89, 89, 89),
            (95, 158, 209),
            (200, 82, 0),
            (137, 137, 137),
            (163, 200, 236),
            (255, 188, 121),
            (207, 207, 207),
        ],
    }
    # normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list]) / 255.0
        colour_array_tup = tuple(map(tuple, colour_array))
        colour_lists[rgb_list] = colour_array_tup
        # add normalised colours to output dictionary
        output_dict[rgb_list] = colour_lists[rgb_list]
    # create a mixed blue/grey colour list, with greys in decreasing darkness
    TUM_colours_list_with_greys = []
    grey = 0.7
    for c in colour_dicts["TUM_colours"].values():
        TUM_colours_list_with_greys.append(f"{grey:0.2f}")
        TUM_colours_list_with_greys.append(c)
        grey -= 0.1
    output_dict["TUM_colours_list_with_greys"] = TUM_colours_list_with_greys

    output_dict["HTML_list01"] = [
        "#808080",
        "#D59460",
        "#005293",
        "#A1B11A",
        "#9ECEEC",
        "#0076B8",
        "#454545",
        "#7b3294",
        "#c2a5cf",
        "#008837",
        "#a6dba0",
    ]
    return output_dict


def get_free_space(folder, format="MB"):
    """
    Return folder/drive free space
    """
    fConstants = {"GB": 1073741824, "MB": 1048576, "KB": 1024, "B": 1}
    if platform.system() == "Windows":
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return (int(free_bytes.value / fConstants[format.upper()]), format)
    else:
        return (int(os.statvfs(folder).f_bfree * os.statvfs(folder).f_bsize / fConstants[format.upper()]), format)


def create_regex_string(inputseq):
    """adds '-*' between each aa or nt/aa in a DNA or protein sequence, so that a particular
    aligned sequence can be identified via a regex search, even if it contains gaps
    inputseq : 'LQQLWNA'
    output   : 'L-*Q-*Q-*L-*W-*N-*A'
    """
    # re.escape each letter. inputseq is a user-supplied sequence, and it is being used here as
    # a regex pattern: an unvalidated "*" or "[" raises re.error from the middle of the pipeline.
    search_string = ""
    for letter in inputseq:
        letter_with_underscore = re.escape(letter) + "-*"
        search_string += letter_with_underscore
    return search_string[:-2]


def convert_truelike_to_bool(input_item, convert_int=False, convert_float=False, convert_nontrue=False):
    """Converts true-like values ("true", 1, True", "WAHR", etc) to python boolean True.

    Taken from korbinian python package by Mark Teese. This is allowed under the permissive MIT license.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "true", 1, "WAHR" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "1.0" will be converted to True
    convert_nontrue : bool
        If True, the output for input_item not recognised as "True" will be False.
        If True, the output for input_item not recognised as "True" will be the original input_item.

    Returns
    -------
    return_value : True, or input_item
        If input_item is True-like, returns python bool True. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_truelike_to_bool("true")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_truelike_to_bool)
    """
    list_True_items = [
        True,
        "True",
        "true",
        "TRUE",
        "T",
        "t",
        "wahr",
        "WAHR",
        "prawdziwy",
        "verdadeiro",
        "sann",
        "istinit",
        "veritable",
        "Pravda",
        "sandt",
        "vrai",
        "igaz",
        "veru",
        "verdadero",
        "sant",
        "gwir",
        "PRAWDZIWY",
        "VERDADEIRO",
        "SANN",
        "ISTINIT",
        "VERITABLE",
        "PRAVDA",
        "SANDT",
        "VRAI",
        "IGAZ",
        "VERU",
        "VERDADERO",
        "SANT",
        "GWIR",
        "bloody oath",
        "BLOODY OATH",
        "nu",
        "NU",
        "damn right",
        "DAMN RIGHT",
    ]

    # if you want to accept 1 or 1.0 as a true value, add it to the list
    if convert_int:
        list_True_items += ["1"]
    if convert_float:
        list_True_items += [1.0, "1.0"]
    # check if the user input string is in the list_True_items
    input_item_is_true = input_item in list_True_items
    # if you want to convert non-True values to "False", then nontrue_return_value = False
    if convert_nontrue:
        nontrue_return_value = False
    else:
        # otherwise, for strings not in the True list, the original string will be returned
        nontrue_return_value = input_item
    # return True if the input item is in the list. If not, return either False, or the original input_item
    return_value = input_item_is_true if input_item_is_true else nontrue_return_value
    # special case: decide if 1 as an integer is True or 1
    if input_item == 1:
        if convert_int:
            return_value = True
        else:
            return_value = 1
    return return_value


def convert_falselike_to_bool(input_item, convert_int=False, convert_float=False):
    """Converts false-like values ("false", 0, FALSE", "FALSCH", etc) to python boolean False.

    Taken from korbinian python package by Mark Teese. This is allowed under the permissive MIT license.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "FALSE", 0, "FALSCH" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "0.0" will be converted to True

    Returns
    -------
    return_value : False, or input_item
        If input_item is False-like, returns python bool False. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_falselike_to_bool("false")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_falselike_to_bool)
    """
    list_False_items = [
        False,
        "False",
        "false",
        "FALSE",
        "F",
        "f",
        "falsch",
        "FALSCH",
        "valse",
        "lažna",
        "fals",
        "NEPRAVDA",
        "falsk",
        "vals",
        "faux",
        "pa vre",
        "tsis tseeb",
        "hamis",
        "palsu",
        "uongo",
        "ngeb",
        "viltus",
        "klaidinga",
        "falz",
        "falso",
        "USANN",
        "wartosc false",
        "falošné",
        "falskt",
        "yanlis",
        "sai",
        "ffug",
        "VALSE",
        "LAŽNA",
        "FALS",
        "FALSK",
        "VALS",
        "FAUX",
        "PA VRE",
        "TSIS TSEEB",
        "HAMIS",
        "PALSU",
        "UONGO",
        "NGEB",
        "VILTUS",
        "KLAIDINGA",
        "FALZ",
        "FALSO",
        "WARTOSC FALSE",
        "FALOŠNÉ",
        "FALSKT",
        "YANLIS",
        "SAI",
        "FFUG",
    ]

    # if you want to accept 0 or 0.0 as a false value, add it to the list
    if convert_int:
        list_False_items += [0, "0"]
    if convert_float:
        list_False_items += [0.0, "0.0"]
    # return boolean False if the input item is in the list. If not, return the original input_item
    return_value = False if input_item in list_False_items else input_item

    return return_value


class HardDriveSpaceException(Exception):
    def __init__(self, value):
        self.parameter = value

    def __str__(self):
        # name changed to allow p-r( to be unique to the print function
        canonical_string_representation = repr
        return canonical_string_representation(self.parameter)


class LogOnlyToConsole:
    def __init__(self):
        pass

    def info(self, message):
        sys.stdout.write(f"\n{message}")

    def warning(self, message):
        sys.stdout.write(f"\n{message}")

    def critical(self, message):
        sys.stdout.write(f"\n{message}")


def open_csv_as_series(input_csv):
    # extract name and sequences from input csv
    input_df = pd.read_csv(input_csv, header=None, index_col=0)
    input_df.columns = ["data"]
    input_ser = input_df["data"]
    return input_ser


def get_testsetname_trainsetname_from_run_settings(s):
    test_set_list, train_set_list = get_test_and_train_set_lists(s)
    assert len(test_set_list) == 1
    assert len(train_set_list) == 1
    testsetname = f"set{int(test_set_list[0]):02d}"
    trainsetname = f"set{int(train_set_list[0]):02d}"
    return testsetname, trainsetname


def get_test_and_train_set_lists(s):
    if s["test_datasets"] is True:
        test_set_list = ["1"]
    elif s["test_datasets"] is False:
        test_set_list = ["0"]
    elif isinstance(s["test_datasets"], int):
        test_set_list = [str(s["test_datasets"])]
    elif isinstance(s["test_datasets"], str):
        test_set_list = s["test_datasets"].split(",")
    else:
        raise ValueError(
            "test_datasets type is not correct {} ({})".format(s["test_datasets"], type(s["test_datasets"]))
        )

    if s["train_datasets"] is True:
        train_set_list = ["1"]
    elif s["train_datasets"] is False:
        train_set_list = ["0"]
    elif isinstance(s["train_datasets"], int):
        train_set_list = [str(s["train_datasets"])]
    elif isinstance(s["train_datasets"], str):
        train_set_list = s["train_datasets"].split(",")
    else:
        raise ValueError(
            "train_datasets type is not correct {} ({})".format(s["train_datasets"], type(s["train_datasets"]))
        )

    return test_set_list, train_set_list


class SurroundingSequence:
    def __init__(self, tmd_start: int, tmd_end: int, full_seq_len: int, num_of_sur_residues: int):
        # int(): the training pipeline reads these out of a DataFrame column, where they are
        # numpy floats, and every offset computed from them was a float too. The annotations here
        # have always said int, and the offsets are used to slice sequences, so a float offset
        # raises TypeError as soon as anything indexes with it rather than with .iloc.
        self.tmd_start: int = int(tmd_start)
        self.tmd_end: int = int(tmd_end)
        self.full_seq_len: int = int(full_seq_len)
        self.num_of_sur_residues: int = int(num_of_sur_residues)
        self.tmd_start_pl_surr: int = self.get_tmd_start_pl_surr()
        self.tmd_end_pl_surr: int = self.get_tmd_end_pl_surr()
        self.n_term_offset: int = self.tmd_start - self.tmd_start_pl_surr
        self.c_term_offset: int = self.tmd_end_pl_surr - self.tmd_end

    def get_tmd_end_pl_surr(self):
        tmd_end_pl_surr = self.tmd_end + self.num_of_sur_residues
        if tmd_end_pl_surr > self.full_seq_len:
            tmd_end_pl_surr = self.full_seq_len
        return tmd_end_pl_surr

    def get_tmd_start_pl_surr(self):
        tmd_start_pl_surr = self.tmd_start - self.num_of_sur_residues
        if tmd_start_pl_surr < 1:
            tmd_start_pl_surr = 1
        return tmd_start_pl_surr


# Characters kept in a protein name that becomes part of a filename. Deliberately narrow:
# anything else is replaced, so no name can contain a path separator, a shell metacharacter, or
# anything else that changes the meaning of a path or a command line.
SAFE_FILENAME_CHARACTERS = re.compile(r"[^A-Za-z0-9._-]+")

# A protein name is a label, not a path component. thoipapy appends suffixes of up to about 40
# characters to it (".lipo_seqs_cdhit_output.fas" and friends), and NAME_MAX is 255 on Linux and
# on every filesystem thoipapy is run on, so a name is capped well below that.
MAX_SAFE_NAME_LEN = 60


def safe_filename_component(value: str, max_length: int = MAX_SAFE_NAME_LEN) -> str:
    """Reduce a protein name to something safe to put in a filename.

    Replaces the previous slugify(), which lowercased and dropped dots and underscores, so
    "P21860_ERBB3" became "p21860_erbb3" and "sp|Q8WWF3|SSMM1_HUMAN" collapsed into one run of
    letters. Case and the separators users actually rely on are worth keeping; what has to go is
    anything that could change the meaning of a path.

    The rule is the same one the THOIPA webserver applies at its own boundary, so a name that has
    already passed through that is unchanged here.

    Parameters
    ----------
    value : str
        Text to convert, e.g. a user-supplied protein name or a full FASTA header.
    max_length : int
        Maximum length of the result.

    Returns
    -------
    str
        The name reduced to [A-Za-z0-9._-], truncated, and never empty.
    """
    safe = SAFE_FILENAME_CHARACTERS.sub("_", str(value)).strip("._-")
    # Strip again after truncating: the cut can land mid-separator and leave a trailing "." or
    # "-", which is legal on Linux but awkward on Windows and SMB shares.
    safe = safe[:max_length].strip("._-")
    # "protein" rather than "": an empty component would silently produce filenames that begin
    # with the suffix, and two differently-named proteins would collide.
    return safe or "protein"


def dropna_with_report(df_data: pd.DataFrame, context: str, logging) -> pd.DataFrame:
    """Drop rows containing NaN, reporting exactly what was lost.

    An untargeted dropna() before training silently removed residues whenever a single feature
    was missing, with no count and no indication of which feature was responsible. House rule for
    ML pipelines is to fail on bad data rather than quietly shrink the dataset, so anything beyond
    a small fraction is an error rather than a warning.

    Parameters
    ----------
    df_data : pd.DataFrame
        Training data, one row per residue.
    context : str
        Where this is being called from, used in the log and error messages.
    logging : logging.Logger
        Python object with settings for logging to console and file.

    Returns
    -------
    pd.DataFrame
        df_data with NaN-containing rows removed.

    Raises
    ------
    ValueError
        If more than MAX_FRACTION_ROWS_DROPPED of the rows would be discarded, or if every row
        would be discarded.
    """
    n_before = len(df_data)
    if n_before == 0:
        raise ValueError(f"{context}: training data is empty before dropping NaN rows.")

    n_nan_per_col = df_data.isna().sum()
    cleaned = df_data.dropna()
    n_dropped = n_before - len(cleaned)

    if n_dropped == 0:
        return cleaned

    worst = n_nan_per_col[n_nan_per_col > 0].sort_values(ascending=False)
    detail = ", ".join(f"{col}={int(n)}" for col, n in worst.items())
    fraction = n_dropped / n_before

    if len(cleaned) == 0:
        raise ValueError(
            f"{context}: every one of the {n_before} rows contains a NaN, so there is no "
            f"training data left. NaN counts per feature: {detail}"
        )
    if fraction > MAX_FRACTION_ROWS_DROPPED:
        raise ValueError(
            f"{context}: dropping NaN rows would discard {n_dropped} of {n_before} residues "
            f"({fraction:.1%}), above the {MAX_FRACTION_ROWS_DROPPED:.0%} limit. This usually "
            f"means a feature failed to be calculated for some proteins rather than that the data "
            f"is genuinely incomplete. NaN counts per feature: {detail}"
        )

    logging.warning(
        f"{context}: dropped {n_dropped} of {n_before} residues ({fraction:.1%}) containing NaN. "
        f"NaN counts per feature: {detail}"
    )
    return cleaned
