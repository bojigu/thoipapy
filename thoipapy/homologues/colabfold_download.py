"""Fetch homologues from the ColabFold MSA server instead of NCBI BLAST.

Why a third homologue source, after ``nr`` and a local UniRef90 BLAST database:

The tracked alignments date from May 2020 and came from a remote ``nr`` query, which now queues
for hours (a probe of a three-residue query returned RTOE=11165, their own estimate of 3.1 hours).
A local UniRef90 blastp fixes the wait but is still a single-pass search, and single-pass BLAST is
the weakest part of the pipeline: conservation and coevolution features are computed from the
alignment, so alignment depth is the ceiling on everything downstream.

The ColabFold server runs MMseqs2 with ``--num-iterations 3`` against UniRef30 and, optionally, an
environmental database assembled from BFD, MGnify, MetaEuk and SMAG. Three profile iterations
against 30%-identity clusters find homologues a single blastp pass cannot, and the clusters are
then expanded back to their members, so the returned alignment is deeper than the cluster count
suggests. For the 59-residue query of 1xioA4 it returns 387 UniRef hits plus 176 environmental
hits, against 114 HSPs from the 2020 ``nr`` search.

It is also fast: that search completed in 20 seconds, against hours of NCBI queue.

Fair use
--------
api.colabfold.com is a free academic service, explicitly "a limited shared resource only capable
of processing a few thousand MSAs per day", and it asks that queries be serial and come from a
single IP. THOIPA's sets are small -- 40 proteins in set08, 10 in set07 -- so a whole set is three
orders of magnitude inside that budget, and the webserver has never carried enough traffic to come
near it either. The public server is therefore a legitimate default for both, and not only for
one-off local runs.

What makes that acceptable is the access pattern rather than the total: this module submits one
protein at a time and never searches in parallel, which is what the server asks for and what keeps
a low-volume client from looking like a batch job. The server was overwhelmed in August 2025 when
a single user submitted a large batch, so a pipeline that depends on it must tolerate the queue
rather than assume it, hence the retries and the RATELIMIT handling below.

A deployment that did start seeing steady traffic should host its own MMseqs2 server and point
THOIPA_COLABFOLD_HOST at it. Nothing else changes: the protocol is the same. Note that this trades
one cost for another rather than removing it, since the ColabFold databases need roughly as much
disk as the local UniRef90 BLAST database they would replace.
"""

import os
import random
import tarfile
import time
from importlib import metadata
from pathlib import Path
from time import strftime

import requests

from thoipapy.artefacts import ArtefactPaths
from thoipapy.utils import make_sure_path_exists

COLABFOLD_HOST_ENV_VAR = "THOIPA_COLABFOLD_HOST"
COLABFOLD_CONTACT_EMAIL_ENV_VAR = "THOIPA_COLABFOLD_CONTACT_EMAIL"

DEFAULT_HOST = "https://api.colabfold.com"

# The search modes ColabFold's own client sends. "env" adds the environmental database to the
# UniRef30 search; the "-nofilter" variants switch off the MMseqs2 diversity filter.
#
# The filter is not a detail. It applies --qsc 0.8 and --max-seq-id 0.95 to thin near-identical
# hits, which is what an alignment for structure prediction wants and not always what a
# conservation or coevolution feature wants. On set08 it is the reason six proteins came back
# with a shallower alignment than the tracked 2020 nr search; switching it off recovered most of
# that (4r0cA7 went from 2740 unique TMD homologues to 6268, against 3144 from nr).
VALID_MODES = ("env", "all", "env-nofilter", "nofilter")


def mode_searches_env_db(mode: str) -> bool:
    """Whether a mode searches the environmental database as well as UniRef30.

    Both "env" and "env-nofilter" do, so this cannot be a test for equality with "env" -- getting
    it wrong silently discards every environmental hit at the parsing step, long after the search
    that paid for them.
    """
    validate_mode(mode)
    return mode.startswith("env")


def validate_mode(mode: str) -> None:
    """Reject an unknown search mode before it reaches the server.

    The server answers an unrecognised mode by quietly falling back to a default, so an unnoticed
    typo in the settings spreadsheet would produce a whole set of alignments from a search nobody
    chose.
    """
    if mode not in VALID_MODES:
        raise ValueError(f"colabfold_mode {mode!r} is not one of {VALID_MODES}")


# Connect timeouts are set "slightly larger than a multiple of 3", per the requests documentation
# and ColabFold's own client.
REQUEST_TIMEOUT_S = 30.02
DOWNLOAD_TIMEOUT_S = 300.02

# The server returns a ticket immediately and the search runs asynchronously, so this bounds the
# wait for the queue, not the search. A 59-residue query completed in 20 s on an idle server; the
# ceiling is for the days when it is not idle.
MAX_WAIT_S = 3600

MAX_HTTP_RETRIES = 5


def get_colabfold_host() -> str:
    """The MSA server to query. Override to point at a self-hosted MMseqs2 server."""
    return os.environ.get(COLABFOLD_HOST_ENV_VAR, "").strip() or DEFAULT_HOST


def get_user_agent() -> str:
    """Identify THOIPA to the MSA server, as ColabFold asks automated clients to do.

    Their client warns when no user agent is set and says the warning will become an error. An
    unidentified client is also the one that gets blocked first when the server is overloaded,
    because the operators have no way to ask it to stop. Raising here rather than sending an
    anonymous query keeps THOIPA a good citizen of a free service it depends on.
    """
    email = os.environ.get(COLABFOLD_CONTACT_EMAIL_ENV_VAR, "").strip()
    if not email:
        raise ValueError(
            f"{COLABFOLD_CONTACT_EMAIL_ENV_VAR} is not set. The ColabFold MSA server asks "
            f"automated clients to identify themselves with a contact address, and blocks those "
            f"that do not when the queue is deep. Set it in .env, or export it, before running a "
            f"prediction."
        )
    try:
        version = metadata.version("thoipapy")
    except metadata.PackageNotFoundError:
        version = "unknown"
    return f"thoipapy/{version} {email}"


def _post_with_retries(url: str, data: dict, headers: dict, logging) -> dict:
    """POST and return the decoded json, retrying transport errors with a short backoff."""
    error_count = 0
    while True:
        try:
            response = requests.post(url, data=data, timeout=REQUEST_TIMEOUT_S, headers=headers)
        except requests.exceptions.Timeout:
            logging.warning("timeout submitting to the MSA server, retrying")
            continue
        except requests.exceptions.RequestException as e:
            error_count += 1
            logging.warning(f"error talking to the MSA server ({error_count}/{MAX_HTTP_RETRIES}): {e}")
            if error_count >= MAX_HTTP_RETRIES:
                raise
            time.sleep(5)
            continue
        break

    try:
        return response.json()
    except ValueError as e:
        raise RuntimeError(f"MSA server did not reply with json: {response.text[:500]}") from e


def _get_with_retries(url: str, headers: dict, logging) -> requests.Response:
    """GET, retrying transport errors with a short backoff."""
    error_count = 0
    while True:
        try:
            return requests.get(url, timeout=DOWNLOAD_TIMEOUT_S, headers=headers)
        except requests.exceptions.Timeout:
            logging.warning("timeout fetching from the MSA server, retrying")
            continue
        except requests.exceptions.RequestException as e:
            error_count += 1
            logging.warning(f"error talking to the MSA server ({error_count}/{MAX_HTTP_RETRIES}): {e}")
            if error_count >= MAX_HTTP_RETRIES:
                raise
            time.sleep(5)


def submit_msa_job(query_seq: str, mode: str, host: str, headers: dict, logging) -> str:
    """Submit one sequence and return the ticket id, resubmitting while the server rate-limits.

    The query is named "101" because that is what ColabFold's own client sends, and the server's
    a3m output carries the query name through. Nothing downstream reads it.
    """
    query = f">101\n{query_seq}\n"
    out = _post_with_retries(f"{host}/ticket/msa", {"q": query, "mode": mode}, headers, logging)

    waited = 0
    while out.get("status") in ("UNKNOWN", "RATELIMIT"):
        pause = 5 + random.randint(0, 5)
        logging.info(f"MSA server returned {out['status']}, resubmitting in {pause}s")
        time.sleep(pause)
        waited += pause
        if waited > MAX_WAIT_S:
            raise TimeoutError(f"MSA server kept returning {out['status']} for {waited}s")
        out = _post_with_retries(f"{host}/ticket/msa", {"q": query, "mode": mode}, headers, logging)

    status = out.get("status")
    if status == "MAINTENANCE":
        raise RuntimeError("the ColabFold MSA server is undergoing maintenance, try again shortly")
    if status == "ERROR":
        raise RuntimeError(f"the ColabFold MSA server rejected the query: {out}")
    if "id" not in out:
        raise RuntimeError(f"unexpected reply from the MSA server: {out}")

    return out["id"]


def wait_for_msa_job(ticket_id: str, host: str, headers: dict, logging) -> None:
    """Poll a ticket until the search completes, raising if it fails or outlasts MAX_WAIT_S."""
    status = "PENDING"
    waited = 0
    while status in ("UNKNOWN", "RUNNING", "PENDING"):
        pause = 5 + random.randint(0, 5)
        time.sleep(pause)
        waited += pause
        response = _get_with_retries(f"{host}/ticket/{ticket_id}", headers, logging)
        try:
            out = response.json()
        except ValueError as e:
            raise RuntimeError(f"MSA server did not reply with json: {response.text[:500]}") from e
        status = out.get("status", "ERROR")
        if waited > MAX_WAIT_S:
            raise TimeoutError(f"MSA server job {ticket_id} still {status} after {waited}s")

    if status != "COMPLETE":
        raise RuntimeError(f"MSA server job {ticket_id} finished with status {status}")


def download_homologues_from_colabfold(
    acc: str,
    query_seq: str,
    a3m_tar_gz: Path,
    details_txt: Path,
    mode: str,
    logging,
) -> None:
    """Search the ColabFold MSA server for one protein and save the returned archive.

    The archive is stored as the server returned it. It holds ``uniref.a3m``, the environmental
    ``bfd.mgnify30.metaeuk30.smag30.a3m`` when the mode asks for it, and ``msa.sh`` -- the exact
    MMseqs2 command line the server ran. Keeping msa.sh is the point: it is the only record of the
    database versions and search parameters behind an alignment, and it is what a methods section
    has to cite. The BLAST path had no equivalent.

    Parameters
    ----------
    acc : str
        Protein accession (e.g. UniProt acc, or PDB acc plus chain letter).
    query_seq : str
        Sequence to search with. The pipeline sends the TMD plus surrounding residues, exactly as
        the BLAST path does, so the two sources stay comparable.
    a3m_tar_gz : Path
        Where to save the archive returned by the server.
    details_txt : Path
        Text file recording what was searched, when, and against which server.
    mode : str
        ColabFold search mode. "env" searches UniRef30 plus the environmental database with the
        diversity filter on; "all" searches UniRef30 only.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    logging.info(f"{acc} starting download_homologues_from_colabfold")

    validate_mode(mode)
    make_sure_path_exists(a3m_tar_gz, isfile=True)
    make_sure_path_exists(details_txt, isfile=True)

    host = get_colabfold_host()
    headers = {"User-Agent": get_user_agent()}

    start = time.time()
    ticket_id = submit_msa_job(query_seq, mode, host, headers, logging)
    logging.info(f"{acc} MSA server ticket {ticket_id}, waiting for the search to finish")
    wait_for_msa_job(ticket_id, host, headers, logging)

    response = _get_with_retries(f"{host}/result/download/{ticket_id}", headers, logging)
    if response.status_code != 200:
        raise RuntimeError(f"{acc} MSA download failed with HTTP {response.status_code}")
    Path(a3m_tar_gz).write_bytes(response.content)

    duration = time.time() - start

    # Written after the download, not before: a details file for an archive that does not exist
    # is how the BLAST path ended up with orphaned metadata.
    date = strftime("%Y%m%d")
    Path(details_txt).write_text(
        f"acc\t{acc}\n"
        f"download_date\t{date}\n"
        f"database\tcolabfold_{mode}\n"
        f"host\t{host}\n"
        f"ticket\t{ticket_id}\n"
        f"query_len\t{len(query_seq)}\n"
        f"query_seq\t{query_seq}\n"
        f"duration_s\t{duration:0.1f}\n"
    )

    # Fail loudly on an archive that is not a readable tarball, rather than letting the parser
    # discover it later against a protein whose name is no longer on screen.
    try:
        with tarfile.open(a3m_tar_gz) as tar:
            members = tar.getnames()
    except tarfile.TarError as e:
        raise RuntimeError(f"{acc} MSA server returned something that is not a tar archive") from e

    if "uniref.a3m" not in members:
        raise RuntimeError(f"{acc} MSA archive has no uniref.a3m, only {members}")

    logging.info(f"Output file: {a3m_tar_gz} ({', '.join(members)}). (time taken = {duration / 60:0.3f} min)")


def download_homologues_from_colabfold_mult_prot(
    paths: ArtefactPaths,
    df_set,
    mode: str,
    rerun_existing_blast_results: bool,
    logging,
) -> None:
    """Run download_homologues_from_colabfold for a set of proteins, one at a time.

    Serial by design. The server asks for serial queries from a single IP, and a set of 40
    proteins at ~20 s each is a 15-minute job, so there is nothing to gain by pushing.

    Parameters
    ----------
    paths : ArtefactPaths
        Locations of the pipeline's input and output files.
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and
        full-length sequences.
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right',
        'database', 'TMD_seq_pl_surr', ....]
    mode : str
        ColabFold search mode, "env" or "all".
    rerun_existing_blast_results : bool
        Re-download proteins whose archive is already present.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    logging.info("~~~~~~~~~~~~   starting download_homologues_from_colabfold_mult_prot   ~~~~~~~~~~~~")

    # Raise before the first request, not on protein 27 of 40, if the client is unidentified.
    get_user_agent()

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        query_seq = df_set.loc[i, "TMD_seq_pl_surr"]
        database = df_set.loc[i, "database"]

        a3m_tar_gz = paths.colabfold_a3m_tar(database, acc)
        details_txt = Path(str(a3m_tar_gz)[: -len(".a3m.tar.gz")] + "_details.txt")

        if a3m_tar_gz.is_file() and not rerun_existing_blast_results:
            logging.info(f"{acc} download_homologues_from_colabfold skipped (EXISTING a3m.tar.gz FILE)")
            continue

        if a3m_tar_gz.is_file():
            logging.info(f"{acc} starting download_homologues_from_colabfold (EXISTING FILE WILL BE OVERWRITTEN)")

        download_homologues_from_colabfold(acc, query_seq, a3m_tar_gz, details_txt, mode, logging)

    logging.info("~~~~~~~~~~~~   finished download_homologues_from_colabfold_mult_prot   ~~~~~~~~~~~~")
