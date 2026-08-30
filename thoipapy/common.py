import glob
import json
import logging
import logging.config
import os
import platform
import re
import signal
import sys
from pathlib import Path
from time import strftime
import numpy as np
import pandas as pd
import psutil
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import thoipapy
from thoipapy.paths import BASE_DIR, DATA_DIR, SETS_DIR
from thoipapy.utils import convert_truelike_to_bool, convert_falselike_to_bool


# from korbinian.utils import convert_truelike_to_bool, convert_falselike_to_bool
# from thoipapy.utils import convert_truelike_to_bool, convert_falselike_to_bool


    # return s["tmp_surr_left"],s["tmp_surr_right"]


def calc_lipophilicity(seq, method="mean"):
    """ Calculates the average hydrophobicity of a sequence according to the Hessa biological scale.

    Function taken from korbinian.utils (allowed under MIT license).

    Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81

    The Hessa scale has been calculated empirically, using the glycosylation assay of TMD insertion.
    Negative values indicate hydrophobic amino acids with favourable membrane insertion.

    Other hydrophobicity scales are in the settings folder. They can be generated as follows.
    hydrophob_scale_path = "thoipapy/setting/hydrophobicity_scales.csv"
    df_hs = pd.read_csv(hydrophob_scale_path, skiprows=2)
    df_hs.set_index("1aa", inplace=True)
    dict_hs = df_hs.Hessa.to_dict()
    hessa_scale = np.array([value for (key, value) in sorted(dict_hs.items())])
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
     'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
     'W', 'Y']

    Parameters:
    -----------
    seq : string
        Sequence to be analysed. Gaps (-) and unknown amino acids (x) should be ignored.
    method : string
        Method to be used to average the hydrophobicity values over the whole sequence.
        The hydrophobicity score is positive for polar/charged aa, negative for hydrophobic aa.
            "sum" will return the sum of the hydrophobicity scores over the sequence
            "mean" will return the mean of the hydrophobicity scores over the sequence

    Returns:
    --------
    mean hydrophobicity value for the sequence entered

    Usage:
    ------
    from korbinian.utils import calc_lipophilicity
    # for a single sequence
    s = "SAESVGEVYIKSTETGQYLAG"
    calc_lipophilicity(s)
    # for a series of sequences
    TMD_ser = df2.TM01_SW_match_seq.dropna()
    hydro = TMD_ser.apply(lambda x : calc_lipophilicity(x))

    Notes:
    ------
    %timeit results:
    for a 20aa seq: 136 µs per loop
    for a pandas series with 852 tmds: 118 ms per loop
    """
    # hydrophobicity scale
    hessa_scale = np.array([0.11, -0.13, 3.49, 2.68, -0.32, 0.74, 2.06, -0.6, 2.71,
                            -0.55, -0.1, 2.05, 2.23, 2.36, 2.58, 0.84, 0.52, -0.31,
                            0.3, 0.68])
    # convert to biopython analysis object
    analysed_seq = ProteinAnalysis(seq)
    # biopython count_amino_acids returns a dictionary.
    aa_counts_dict = analysed_seq.count_amino_acids()
    # get the number of AA residues used to calculated the hydrophobicity
    # this is not simply the sequence length, as the sequence could include gaps or non-natural AA
    aa_counts_excluding_gaps = np.array(list(aa_counts_dict.values()))
    number_of_residues = aa_counts_excluding_gaps.sum()
    # if there are no residues, don't attempt to calculate a mean. Return np.nan.
    if number_of_residues == 0:
        return np.nan
    # convert dictionary to array, sorted by aa
    aa_counts_arr = np.array([value for (key, value) in sorted(aa_counts_dict.items())])
    multiplied = aa_counts_arr * hessa_scale
    sum_of_multiplied = multiplied.sum()
    if method == "mean":
        return sum_of_multiplied / number_of_residues
    if method == "sum":
        return sum_of_multiplied



def _coerce_setting(value, declared_type, parameter: str):
    """Convert a CSV settings value to the type the pipeline expects.

    Everything in a CSV is a string, and the old Excel loader relied on Excel having already
    typed the cells. Without this, "0" would arrive as a non-empty string and therefore be
    truthy, so a setting of 0 would behave as if it were 1.

    The `type` column is used when present, but is not trusted to be complete: it only ever
    existed on one of the three original sheets and was blank on several rows there. Values
    without a usable declared type are inferred.

    Parameters
    ----------
    value : Any
        Raw value from the CSV.
    declared_type : str or None
        Contents of the `type` column, e.g. "bool", "int", "float".
    parameter : str
        Parameter name, used in the error message.

    Returns
    -------
    bool, int, float or str
    """
    if not isinstance(value, str):
        return value
    text = value.strip()
    declared = (declared_type or "").strip().lower() if isinstance(declared_type, str) else ""

    if declared == "bool":
        # convert_int so that 0/1, which is how most of these flags are actually written, is read
        # as a boolean rather than inferred as an integer
        as_bool = convert_falselike_to_bool(
            convert_truelike_to_bool(text, convert_int=True, convert_float=True),
            convert_int=True, convert_float=True)
        if not isinstance(as_bool, bool):
            raise ValueError(f"setting {parameter!r} is declared bool but has the value {value!r}")
        return as_bool
    if declared == "int":
        return int(float(text))
    if declared == "float":
        return float(text)

    # no usable declared type: infer, most specific first
    as_bool = convert_falselike_to_bool(convert_truelike_to_bool(text))
    if isinstance(as_bool, bool):
        return as_bool
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def create_settingdict(settings_path):
    """Load the settings CSV into a dictionary.

    The settings used to live in a multi-sheet Excel workbook. CSV replaces it because the
    settings are the pipeline's most important input and were the hardest thing to inspect: a
    workbook cannot be diffed in git, cannot be edited without Excel, and programmatic row
    deletion silently shifted values against the wrong parameters. A flat CSV with a `section`
    column carries exactly the same information and can be read by anything.

    Parameters
    ----------
    settings_path : str or Path
        Path to the settings CSV.

    Returns
    -------
    dict
        Settings, keyed by parameter name.

    Raises
    ------
    FileNotFoundError
        If the settings file does not exist.
    ValueError
        If a parameter appears more than once, which would silently make one value win.
    """
    settings_path = Path(settings_path)
    if not settings_path.is_file():
        raise FileNotFoundError(f"settings file not found: {settings_path}")

    dfset = pd.read_csv(settings_path)
    dfset = dfset.dropna(subset=["parameter", "value"])

    duplicated = dfset["parameter"][dfset["parameter"].duplicated()].tolist()
    if duplicated:
        raise ValueError(f"{settings_path.name} defines these parameters more than once: {duplicated}")

    s = {"settings_path": settings_path}
    types = dfset.set_index("parameter")["type"].to_dict() if "type" in dfset.columns else {}
    for parameter, value in dfset.set_index("parameter")["value"].to_dict().items():
        s[parameter] = _coerce_setting(value, types.get(parameter), parameter)

    list_paths_to_normalise = ['MiRMAK_data_folder', 'base_dir', 'sets_dir', 'data_dir',
                               'Rcode', 'hhblits_dir', 'uniprot_database_dir', 'Rscript_dir']
    # normalise the paths for selected columns, so that they are appropriate for the operating system
    for path in list_paths_to_normalise:
        if path in s:
            s[path] = os.path.normpath(s[path])
            # if not os.path.exists(s[path]):
            #     os.makedirs(s[path])

    # Convert every remaining true/false-like string, regardless of the sheet's "type" column.
    # That column only exists on the run_settings sheet, and eight rows there leave it blank, so
    # values like the string "TRUE" reached the pipeline unconverted. They happened to work
    # because a non-empty string is truthy -- but so is "FALSE", which meant a flag switched off
    # in the spreadsheet could silently behave as if it were on.
    for key, value in s.items():
        if isinstance(value, str):
            as_bool = convert_falselike_to_bool(convert_truelike_to_bool(value))
            if isinstance(as_bool, bool):
                s[key] = as_bool

    # The data, sets and base directories are properties of this repository, not user settings.
    # Neither shipped settings file has ever defined them, so before this they were simply absent
    # unless a caller injected them, and the pipeline could not be run from a clean checkout.
    # Set after the loop above so they stay Path objects; os.path.normpath would return str.
    s["data_dir"] = DATA_DIR
    s["sets_dir"] = SETS_DIR
    s["base_dir"] = BASE_DIR

    return s


def setup_keyboard_interrupt_and_error_logging(s, setname):
    ''' -------Setup keyboard interrupt----------
        Taken from korbinian python package by Mark Teese. This is allowed under the permissive MIT license.
    '''

    # import arcgisscripting

    def ctrlc(sig, frame):
        raise KeyboardInterrupt("CTRL-C!")

    signal.signal(signal.SIGINT, ctrlc)
    '''+++++++++++++++LOGGING++++++++++++++++++'''
    date_string = strftime("%Y%m%d_%H_%M_%S")

    # designate the output logfile
    logfile = os.path.join(s["data_dir"], "Logging", '%s_%s_logfile.log' % (setname, date_string))

    # # if multiprocessing is used, disable logging except for critical messages.
    # if s["use_multiprocessing"]:
    #     level_console = "CRITICAL"
    #     level_logfile = "CRITICAL"
    # else:
    #     level_console = s["logging_level_console"]
    #     level_logfile = s["logging_level_logfile"]

    level_console = s["logging_level_console"]
    level_logfile = s["logging_level_logfile"]

    logging = thoipapy.common.setup_error_logging(logfile, level_console, level_logfile)
    return logging


def setup_error_logging(logfile, level_console="DEBUG", level_logfile="DEBUG", print_system_info=True):
    """ Sets up error logging, and logs a number of system settings.

    Taken from korbinian python package by Mark Teese. This is allowed under the permissive MIT license.

    Parameters:
    -----------
    logfile : str
        Path to output logfile. If size exceeds limit set below in JSON settings, path.1, path.2 etc will be created.
    level_console : str
        Logging level for printing to console. DEBUG, WARNING or CRITICAL
    level_logfile : str
        Logging level for printing to logfile. DEBUG, WARNING or CRITICAL
    """
    # load the log settings in json format
    logsettings = json.dumps({
        "handlers": {
            "console": {
                "formatter": "brief",
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
                "level": "DEBUG"
            },
            "file": {
                "maxBytes": 10000000,
                "formatter": "precise",
                "backupCount": 3,
                "class": "logging.handlers.RotatingFileHandler",
                "level": "DEBUG",
                "filename": "logfile.txt"
            }
        },
        "version": 1,
        "root": {
            "handlers": [
                "console",
                "file"
            ],
            "propagate": "no",
            "level": "DEBUG"
        },
        "formatters": {
            "simple": {
                "format": "format=%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            },
            "precise": {
                "format": "%(asctime)s %(name)-15s %(levelname)-8s %(message)s"
            },
            "brief": {
                "format": "%(levelname)-8s: %(name)-15s: %(message)s"
            }
        }
    }, skipkeys=True, sort_keys=True, indent=4, separators=(',', ': '))

    config = json.loads(logsettings)
    # add user parameters to the logging settings (logfile, and logging levels)
    config['handlers']['file']['filename'] = logfile
    config['handlers']['console']['level'] = level_console
    config['handlers']['file']['level'] = level_logfile

    # create a blank logging file
    thoipapy.utils.make_sure_path_exists(logfile, isfile=True)
    with open(logfile, 'w') as f:
        pass

    # clear any previous logging handlers that might have been previously run in the console
    logging.getLogger('').handlers = []
    # load the logging settings from the modified json string
    logging.config.dictConfig(config)
    # collect a number of system settings that could be useful for troubleshooting
    system_settings_dict = {}
    system_settings_dict["system description"] = platform.uname()
    system_settings_dict["system"] = platform.system()
    system_settings_dict["architecture"] = platform.architecture()
    system_settings_dict["network_name"] = platform.node()
    system_settings_dict["release"] = platform.release()
    system_settings_dict["version"] = platform.version()
    system_settings_dict["machine"] = platform.machine()
    system_settings_dict["processor"] = platform.processor()
    system_settings_dict["python_version"] = platform.python_version()
    system_settings_dict["python_build"] = platform.python_build()
    system_settings_dict["python_compiler"] = platform.python_compiler()
    system_settings_dict["argv"] = sys.argv
    system_settings_dict["dirname(argv[0])"] = os.path.abspath(os.path.expanduser(os.path.dirname(sys.argv[0])))
    system_settings_dict["pwd"] = os.path.abspath(os.path.expanduser(os.path.curdir))
    system_settings_dict["total_ram"] = "{:0.2f} GB".format(psutil.virtual_memory()[0] / 1000000000)
    system_settings_dict["available_ram"] = "{:0.2f} GB ({}% used)".format(psutil.virtual_memory()[1] / 1000000000, psutil.virtual_memory()[2])
    # log the system settings
    if print_system_info:
        logging.warning(system_settings_dict)
    # test error message reporting
    # logging.warning('LOGGING TEST:')
    # try:
    #    open('/path/to/does/not/exist', 'rb')
    # except (SystemExit, KeyboardInterrupt):
    #    raise
    # except Exception:
    #    logging.error('Failed to open file', exc_info=True)
    logging.warning('LOGGING SETUP IS SUCCESSFUL (logging levels: console={}, logfile={}). \n'.format(level_console, level_logfile))
    return logging


def get_path_of_protein_set(setname, sets_dir):
    """Get the path of a protein set by searching the sets folder for its name.

    Parameters
    ----------
    setname : str
        Name of the protein set. E.g. set08
    sets_dir : str
        Path to protein set folder

    Returns
    -------
    set_path : str
        Path to particular protein set.
    """
    set_file_list = glob.glob(os.path.join(sets_dir, "*.csv"))

    matching = [set_path for set_path in set_file_list if setname in os.path.basename(set_path)]
    if len(matching) == 1:
        return matching[0]
    if len(matching) == 0:
        raise FileNotFoundError(f"No protein set file found for setname '{setname}'.\nFiles in {sets_dir}: {[os.path.basename(p) for p in set_file_list]}")
    raise ValueError(f"More than one file in the set folder contains '{setname}' in the filename.\nmatching = {matching}")


def process_set_protein_seqs(s, setname, df_set, set_path):
    df_set["seqlen"] = df_set.full_seq.str.len()
    df_set["TMD_len"] = df_set.TMD_seq.str.len()
    df_set["acc_db"] = df_set.acc + "-" + df_set.database

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        TMD_seq = df_set.loc[i, "TMD_seq"]
        full_seq = df_set.loc[i, "full_seq"]

        # use regex to get indices for start and end of TMD in seq
        m = re.search(TMD_seq, full_seq)
        if m:
            # convert from python indexing to unprot indexing
            df_set.loc[i, "TMD_start"] = m.start() + 1
            df_set.loc[i, "TMD_end"] = m.end()
        else:
            raise IndexError("TMD seq not found in full_seq.\nacc = {}\nTMD_seq = {}\nfull_seq = {}".format(acc, TMD_seq, full_seq))

    # first get TMD plus 5 surrounding residues (for TMD_lipo script)
    num_of_sur_residues = 5
    # TODO to improve consistency, replace create_column_with_TMD_plus_surround_seq with new
    # thoipapy.utils.SurroundingSequence developed for the standalone predictor
    df_set, TMD_seq_pl_surr_series = thoipapy.utils.create_column_with_TMD_plus_surround_seq(df_set, num_of_sur_residues)
    df_set["TMD_seq_pl_surr5"] = TMD_seq_pl_surr_series

    # Repeat for the actual surrounding number of residues chosen in the settings file
    # this overwrites the indexing columns created for the surr5 above, except for the final sequence
    num_of_sur_residues = s["num_of_sur_residues"]
    df_set, TMD_seq_pl_surr_series = thoipapy.utils.create_column_with_TMD_plus_surround_seq(df_set, num_of_sur_residues)
    df_set["TMD_seq_pl_surr"] = TMD_seq_pl_surr_series

    # add the number of included residues in the surrounding seq to the left and right of the TMD
    # e.g. 20 where the TMD is in the centre of the protein, otherwise <20 where TMD is near start or end of full seq
    df_set["tm_surr_left"] = df_set.TMD_start - df_set.TMD_start_pl_surr
    df_set["tm_surr_right"] = df_set.TMD_end_pl_surr - df_set.TMD_end

    # save the full sequences in fasta format for CD-HIT, etc.
    protein_set_full_seq_fasta = Path(s["data_dir"]) / f"results/{s['setname']}/clusters/{setname}_full_seqs.fas"
    thoipapy.utils.make_sure_path_exists(protein_set_full_seq_fasta, isfile=True)
    with open(protein_set_full_seq_fasta, "w") as f:
        for n, acc in enumerate(df_set.index):
            f.write(">{}-{}\n{}\n".format(n, df_set.loc[acc, "acc_db"], df_set.loc[acc, "full_seq"]))

    protein_set_tmd_seq_fasta = Path(s["data_dir"]) / f"results/{s['setname']}/clusters/{setname}_tmd_seqs.fas"
    thoipapy.utils.make_sure_path_exists(protein_set_tmd_seq_fasta, isfile=True)
    with open(protein_set_tmd_seq_fasta, "w") as f:
        for n, acc in enumerate(df_set.index):
            f.write(">{}-{}\n{}\n".format(n, df_set.loc[acc, "acc_db"], df_set.loc[acc, "TMD_seq"]))

    # open previously saved CD-hit results
    cdhit_cluster_txt = Path(s["data_dir"]) / f"results/{s['setname']}/clusters/{setname}.fas.1.clstr.sorted.txt"
    if os.path.isfile(cdhit_cluster_txt):
        lines_with_ref_seq = []
        with open(cdhit_cluster_txt, "r") as f:
            for line in f:
                if "*" in line:
                    lines_with_ref_seq.append(line)

        cluster_rep_lines_ser = pd.Series(lines_with_ref_seq)
        # extracts the number between > and - (i.e., the index), sorts and returns as a list
        # redundant sequences will be excluded
        cluster_rep_list = cluster_rep_lines_ser.str.extract(">(\d*)-", expand=False).astype(int).sort_values().tolist()
        df_set.loc[cluster_rep_list, "cdhit_cluster_rep"] = True
        df_set["cdhit_cluster_rep"] = df_set["cdhit_cluster_rep"].fillna(False)
    else:
        logging.warning("No CD-HIT results found for automatic redundancy reduction. It is assumed that dataset is non-redundant. Further CD-HIT clustering may be used for predictor validation.")
        df_set["cdhit_cluster_rep"] = "no_cdhit_results"

    """  Rearrange the dataframe columns so that the order is as follows.
    orig Bo file : ['acc', 'TMD_Length', 'TMD_Start', 'TMD_End', 'TMD_Sur_Left', 'TMD_Sur_Right']
    updated file = ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]

    """
    # reorder columns
    df_set = thoipapy.utils.reorder_dataframe_columns(df_set, ['acc', 'seqlen', 'TMD_start', 'TMD_end', "tm_surr_left", "tm_surr_right", "database"])

    # Residue positions are whole numbers and are written to the processed CSV, so they must not
    # be serialised as "192.0". Assigning back through .iloc no longer changes the dtype under
    # pandas copy-on-write, which silently made this a no-op; name the columns explicitly instead
    # of relying on their position.
    for column in ["seqlen", "TMD_start", "TMD_end", "tm_surr_left", "tm_surr_right"]:
        df_set[column] = df_set[column].astype(int)

    # save to csv, which is opened by other functions
    list_of_tmd_start_end = os.path.join(s["data_dir"], "Input_data", os.path.basename(set_path)[:-5] + "_processed.csv")
    s["list_of_tmd_start_end"] = list_of_tmd_start_end
    thoipapy.utils.make_sure_path_exists(list_of_tmd_start_end, isfile=True)
    df_set.set_index("acc").to_csv(list_of_tmd_start_end)

    return df_set
