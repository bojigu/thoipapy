from pathlib import Path
from typing import Union

from Bio.Blast import NCBIWWW
from thoipapy.utils import get_free_space, HardDriveSpaceException
import os
import subprocess
import tarfile
import platform
import time
from time import strftime
from thoipapy.utils import delete_BLAST_xml, make_sure_path_exists
from thoipapy.artefacts import ArtefactPaths


def download_homologues_from_ncbi_mult_prot(paths: ArtefactPaths, df_set, expect_value, hit_list_size: int, rerun_existing_blast_results: bool, logging):
    """Runs download_homologues_from_ncbi for a set of proteins

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
    hit_list_size = int(hit_list_size)

    OS_description = platform.system()

    # try to detect if there's not enough HD space for the download.
    # currently not working for
    if "Linux" in OS_description or "Windows" in OS_description:
        try:
            byteformat = "GB"
            data_dir = paths.data_dir

            size = get_free_space(data_dir, byteformat)
            # logging.info('Hard disk remaining space = {}'.format(size))

            if size[0] < 5:
                raise HardDriveSpaceException(
                    "Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))

        # log the exception, so you actually can see what goes on in the logfile
        except HardDriveSpaceException as e:
            logging.critical(e)
            # now stop all processes and raise an error
            raise HardDriveSpaceException("process stopped")
    else:
        logging.warning("Your system does not seem to be Linux or Windows. Harddrive space check not conducted.")

    ##############################################################################################
    #                                                                                            #
    #      download homologues from protein lists file                                           #
    #                                                                                            #
    ##############################################################################################

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        TMD_seq_pl_surr = df_set.loc[i, "TMD_seq_pl_surr"]
        # TMD_seq_pl_surr = df_set.loc[i, "full_seq"]
        database = df_set.loc[i, "database"]

        # run online server NCBI blastp with biopython module
        blast_xml_file = str(paths.blast_xml(database, acc))
        xml_tar_gz = blast_xml_file[:-4] + ".xml.tar.gz"
        xml_txt = blast_xml_file[:-4] + "_details.txt"

        if not os.path.isfile(xml_tar_gz):
            run_download = True
        else:
            if rerun_existing_blast_results:
                run_download = True
                logging.info('{} starting download_homologues_from_ncbi_mult_prot (EXISTING xml.tar.gz FILE WILL BE OVERWRITTEN)'.format(acc))
            elif rerun_existing_blast_results in [False, 0]:
                run_download = False
                logging.info('{} download_homologues_from_ncbi_mult_prot skipped (EXISTING xml.tar.gz FILE)'.format(acc))
                # skip protein
                continue
            else:
                raise ValueError('rerun_existing_blast_results does not seem to be True or False')

        if run_download:
            download_homologues_from_ncbi(acc, TMD_seq_pl_surr, blast_xml_file, xml_txt, xml_tar_gz, expect_value, hit_list_size, logging)

    logging.info('homologues download was finished')


NCBI_CONTACT_EMAIL_ENV_VAR = "THOIPA_NCBI_CONTACT_EMAIL"

# Point this at a local BLAST database (the -db argument, i.e. the path without the file
# extensions) to search it with a local blastp instead of querying NCBI over the network.
LOCAL_BLAST_DB_ENV_VAR = "THOIPA_LOCAL_BLAST_DB"
LOCAL_BLAST_TIMEOUT_S = 1800


def get_local_blast_db() -> str:
    """Return the configured local BLAST database, or "" if predictions should query NCBI."""
    return os.environ.get(LOCAL_BLAST_DB_ENV_VAR, "").strip()


def run_local_blastp(query_fasta_string: str, blast_xml_file: Union[Path, str], expect_value: float,
                     hit_list_size: int, db: str, logging) -> None:
    """Search a local BLAST database with blastp, writing NCBI-format XML.

    Output format 5 is not optional: NCBI_parser reads BLAST XML, so a local search has to produce
    exactly what the web service would have. -evalue and -max_target_seqs are the local spellings
    of the expect_value and hit_list_size settings, so behaviour matches the remote path.

    A local search of Swiss-Prot returns in well under a second, against a wait at NCBI that has
    been measured in hours when their queue is deep. The trade is alignment depth: Swiss-Prot is
    roughly a quarter of what nr yields, which is fine for tests and not enough for production.
    """
    query_file = Path(blast_xml_file).with_suffix(".query.fasta")
    query_file.write_text(query_fasta_string + "\n")

    cmd = [
        "blastp",
        "-query", str(query_file),
        "-db", db,
        "-outfmt", "5",
        "-evalue", str(expect_value),
        "-max_target_seqs", str(hit_list_size),
        "-num_threads", str(min(8, os.cpu_count() or 1)),
        "-out", str(blast_xml_file),
    ]
    logging.info(f"running local blastp against {db}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=LOCAL_BLAST_TIMEOUT_S)
        if result.returncode != 0:
            raise RuntimeError(f"local blastp failed (returncode={result.returncode}): {result.stderr.strip()}")
    finally:
        # Ran only on success before, so every failed search left a stray .query.fasta behind.
        query_file.unlink(missing_ok=True)


def identify_to_ncbi() -> None:
    """Tell NCBI who is making the BLAST request, as their usage policy requires.

    NCBI asks that automated clients send a contact address, and throttles or blocks those that do
    not: repeated identical queries from an unidentified client get refused outright or pushed
    into a slow queue. Biopython's qblast takes no email argument, but it reads the module-level
    `NCBIWWW.email` and passes it through as a URL parameter alongside `tool`, so setting it here
    is the supported route.

    Raises if unset rather than quietly sending an anonymous query, because an anonymous query is
    what gets the whole server throttled.
    """
    email = os.environ.get(NCBI_CONTACT_EMAIL_ENV_VAR, "").strip()
    if not email:
        raise ValueError(
            f"{NCBI_CONTACT_EMAIL_ENV_VAR} is not set. NCBI requires automated BLAST clients to "
            f"identify themselves with a contact address, and throttles those that do not. "
            f"Set it in .env, or export it, before running a prediction."
        )
    NCBIWWW.email = email
    NCBIWWW.tool = "thoipapy"


# The remote database stays "nr". UniRef is a UniProt resource and NCBI does not host it: the
# BLAST URL API serves nr, refseq_protein, swissprot, pdb, env_nr, tsa_nr, landmark and pataa.
# Of those only nr has depth comparable to UniRef90, and NCBI's swissprot is as shallow as a local
# Swiss-Prot build. To search UniRef90, build it locally and set THOIPA_LOCAL_BLAST_DB.
def download_homologues_from_ncbi(acc: str, TMD_seq_pl_surr: str, blast_xml_file: Union[Path, str], xml_txt: Union[Path, str], xml_tar_gz: Union[Path, str], expect_value: float, hit_list_size: int, logging, db: str = "nr"):
    """Downloads homologue xml file using NCBI BLAST via the biopython qBLAST wrapper.

    Parameters
    ----------
    acc : str
        Protein accession number. (e.g. UniProt acc, PDB_acc + chain letter)
    TMD_seq_pl_surr : str
        TMD sequence plus surrouding residues (usually 20) for BLAST
    blast_xml_file : str
        Path to xml file with BLAST results from NCBI.
    xml_txt : str
        Path to a text file that saves the download date, etc.
    xml_tar_gz : str
        Path to compressed tar file containing blast_xml_file, with BLAST results from NCBI.
    expect_value : float or int
        BLAST parameter "expect", which gives the expected number of random matches after BLAST with these parameters.
        Typically 10.
    hit_list_size : int
        Maximum number of BLAST hits.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    db : database
        NCBI database to search.
        Default is "nr"
    """
    logging.info('{} starting download_homologues_from_ncbi'.format(acc))

    # create an empty text file with the download date
    date = strftime("%Y%m%d")
    database_searched = get_local_blast_db() or f"ncbi_{db}"
    with open(xml_txt, "w") as f:
        f.write("acc\t{}\ndownload_date\t{}\ndatabase\t{}\nexpect_value\t{}\n".format(acc, date, database_searched, expect_value))

    query_fasta_string = ">{} TMD add surround 20 residues\n{}".format(acc, TMD_seq_pl_surr)
    make_sure_path_exists(blast_xml_file, isfile=True)

    local_db = get_local_blast_db()

    start = time.time()

    if not local_db:
        # Outside the try: a missing contact email is a configuration error, and reporting it as
        # "blast_xml_file not found" further down sends the reader after the wrong problem.
        identify_to_ncbi()

    try:
        if local_db:
            run_local_blastp(query_fasta_string, blast_xml_file, expect_value, hit_list_size, local_db, logging)
        else:
            tmp_protein_homologues_xml_handle = NCBIWWW.qblast("blastp", db, query_fasta_string,
                                                               expect=expect_value,
                                                               hitlist_size=hit_list_size)
            with open(blast_xml_file, "w") as save_tmp_xml_file:
                save_tmp_xml_file.write(tmp_protein_homologues_xml_handle.read())

            tmp_protein_homologues_xml_handle.close()

    except Exception:
        # Was a bare `except`. A hang is not an exception, so this never caught the failure mode
        # that matters; what it did catch, it caught silently.
        logging.exception(f"{acc} BLAST search failed. "
                          f"A typical error message is 'query string not found in the CGI context in qblast'. "
                          f"A refusal shortly after an identical query usually means NCBI is throttling "
                          f"this client.")

    duration = time.time() - start

    if os.path.isfile(blast_xml_file):
        with tarfile.open(xml_tar_gz, mode='w:gz') as tar:
            # add the files to the compressed tarfile
            tar.add(blast_xml_file, arcname=os.path.basename(blast_xml_file))
            tar.add(xml_txt, arcname=os.path.basename(xml_txt))

        delete_BLAST_xml(blast_xml_file)
    else:
        raise Exception(f"{acc} download_homologues_from_ncbi failed, blast_xml_file not found ({blast_xml_file})")

    logging.info("Output file: {}. (time taken = {:0.3f} min)".format(xml_tar_gz, duration / 60))


def download_10_homologues_from_ncbi(paths: ArtefactPaths, df_set, rerun_existing_blast_results: bool, logging):
    """Runs download_homologues_from_ncbi for a set of proteins

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
    expect_value = 1
    hit_list_size = 10

    ##############################################################################################
    #                                                                                            #
    #      download homologues from protein lists file                                           #
    #                                                                                            #
    ##############################################################################################

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        full_seq = df_set.loc[i, "full_seq"]
        database = df_set.loc[i, "database"]

        # run online server NCBI blastp with biopython module
        blast_xml_file = str(paths.data_dir / "homologues" / "xml" / "10_hits" / database / f"{acc}.{paths.surr_suffix}.BLAST.xml")
        xml_tar_gz = blast_xml_file[:-4] + ".xml.tar.gz"
        xml_txt = blast_xml_file[:-4] + "_details.txt"

        if not os.path.isfile(xml_tar_gz):
            run_download = True
        else:
            if rerun_existing_blast_results:
                run_download = True
                logging.info('{} starting download_homologues_from_ncbi_mult_prot (EXISTING xml.tar.gz FILE WILL BE OVERWRITTEN)'.format(acc))
            elif rerun_existing_blast_results in [False, 0]:
                run_download = False
                logging.info('{} download_homologues_from_ncbi_mult_prot skipped (EXISTING xml.tar.gz FILE)'.format(acc))
                # skip protein
                continue
            else:
                raise ValueError('rerun_existing_blast_results does not seem to be True or False')

        if run_download:
            download_homologues_from_ncbi(acc, full_seq, blast_xml_file, xml_txt, xml_tar_gz, expect_value, hit_list_size, logging)

    logging.info('homologues download was finished')
