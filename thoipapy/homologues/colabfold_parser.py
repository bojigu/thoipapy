"""Turn the ColabFold MSA server's a3m output into the homologue csv the pipeline already reads.

The point of this module is that nothing downstream changes. ``extract_filtered_csv_homologues_to
_alignments`` reads four columns out of the homologue csv -- ``query_align_seq``,
``subject_align_seq``, ``match_markup_seq`` and ``description`` -- slices the TMD out of them with
a regex, and filters on gaps and identity. Produce those four columns from an a3m and every
feature downstream is computed by the same code from a deeper alignment.

a3m as a pairwise alignment
---------------------------
An a3m is query-anchored, which is exactly the shape BLAST HSPs are already reduced to. Each
record is a string over three kinds of character:

    uppercase   aligned to the query residue in that column
    '-'         the subject has nothing aligned to that query residue
    lowercase   an insertion in the subject, with no query residue opposite it

So the match-state columns (uppercase and '-') number exactly the query length, and the pairwise
alignment comes back by walking the string once: a lowercase character puts a gap in the query and
an uppercase one opposite it; anything else consumes the next query residue. This was verified
against all 561 records returned for 1xioA4, including the 211 with insertions and the 387 with
partial coverage: the reconstructed identity matches the server's own ``fident`` field and the
reconstructed span matches its qstart/qend, for every record.

The header carries the rest. ``--msa-format-mode 6`` writes a tab-separated line:

    target  score  fident  evalue  qstart  qend  qlen  tstart  tend  tlen

with coordinates 0-based and inclusive, so they are converted to BLAST's 1-based coordinates here.
That is more than the XML path had -- BLAST XML gave no per-hit identity fraction -- so no
information is lost by the switch.
"""

import csv
import tarfile
from pathlib import Path

from thoipapy.artefacts import ArtefactPaths

# The a3m files inside the archive returned by the server. uniref.a3m is always present; the
# environmental one only when the search ran in a mode that includes it.
UNIREF_A3M = "uniref.a3m"
ENV_A3M = "bfd.mgnify30.metaeuk30.smag30.a3m"

# The columns the NCBI path writes, so the two sources produce interchangeable files. Sorted,
# because the NCBI parser sorts its header and the downstream reader is a DataFrame either way.
CSV_COLUMNS = sorted(
    [
        "hit_num",
        "description",
        "organism",
        "FASTA_expectation",
        "FASTA_identity",
        "query_align_seq",
        "subject_align_seq",
        "match_markup_seq",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
    ]
)

N_A3M_HEADER_FIELDS = 10


def read_a3m(a3m_text: str) -> list[tuple[str, str]]:
    """Parse a3m text into (header, sequence) pairs, the query first.

    MMseqs2 terminates its database entries with a NUL byte, and the last record of an a3m written
    by ``result2msa`` carries it through into the file. Left in place it makes that one record a
    column longer than the query and the reconstruction below walks off the end of the query, so
    it is stripped here rather than guarded against in every caller.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    for raw_line in a3m_text.splitlines():
        line = raw_line.replace("\x00", "").strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header, seq_parts = line[1:], []
        else:
            seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def pairwise_alignment_from_a3m_record(query_seq: str, subject_a3m_seq: str) -> tuple[str, str, str]:
    """Expand one a3m record into (query_align_seq, subject_align_seq, match_markup_seq).

    The markup uses the residue letter where query and subject agree and a space where they do
    not. BLAST also writes '+' for a conservative substitution, but the only consumer of this
    column replaces '+' with a space before counting, so letter-or-space carries exactly the
    information that is read and none that is not.
    """
    query_chars: list[str] = []
    subject_chars: list[str] = []
    query_index = 0

    for char in subject_a3m_seq:
        if char.islower():
            # An insertion in the subject: no query residue is opposite it.
            query_chars.append("-")
            subject_chars.append(char.upper())
        else:
            if query_index >= len(query_seq):
                raise ValueError(
                    f"a3m record has more match-state columns than the {len(query_seq)}-residue "
                    f"query, which means it is not anchored to this query"
                )
            query_chars.append(query_seq[query_index])
            subject_chars.append(char)
            query_index += 1

    if query_index != len(query_seq):
        raise ValueError(
            f"a3m record covers {query_index} of the query's {len(query_seq)} residues; an a3m "
            f"record must have one match-state column per query residue"
        )

    query_align_seq = "".join(query_chars)
    subject_align_seq = "".join(subject_chars)
    markup = "".join(q if q == s and q != "-" else " " for q, s in zip(query_align_seq, subject_align_seq, strict=True))

    return query_align_seq, subject_align_seq, markup


def _rows_from_a3m(a3m_text: str, source: str, e_value_cutoff: float, logging) -> tuple[list[dict], str, int]:
    """Build csv rows from one a3m file. Returns (rows, query_seq, n_excluded_by_e_value)."""
    records = read_a3m(a3m_text)
    if not records:
        raise ValueError(f"{source} a3m is empty")

    query_seq = records[0][1]
    if any(c.islower() or c == "-" for c in query_seq):
        raise ValueError(f"{source} a3m: the first record should be the ungapped query, got {query_seq!r}")

    rows: list[dict] = []
    n_excluded = 0

    for header, subject_a3m_seq in records[1:]:
        fields = header.split("\t")
        if len(fields) != N_A3M_HEADER_FIELDS:
            # A header without the alignment statistics cannot be filtered on e-value, and
            # silently keeping it would put an unfiltered hit into a filtered alignment.
            logging.warning(f"{source}: skipping a3m record with {len(fields)} header fields: {header[:80]}")
            continue

        target, _score, fident, evalue, qstart, qend, _qlen, tstart, tend, _tlen = fields
        expectation = float(evalue)
        if expectation > e_value_cutoff:
            n_excluded += 1
            continue

        query_align_seq, subject_align_seq, markup = pairwise_alignment_from_a3m_record(query_seq, subject_a3m_seq)

        rows.append(
            {
                "description": f"{target} {source}",
                # a3m carries no organism, and the NCBI path wrote "no_organism" unconditionally
                # anyway -- it assigned the parsed taxonomy and then overwrote it on the next line.
                "organism": "no_organism",
                "FASTA_expectation": expectation,
                # The real identity fraction from MMseqs2. The NCBI path wrote
                # `hsp.identities / 100`, an identity *count* divided by 100, which is not a
                # fraction of anything; nothing reads the column, so the bug never surfaced.
                "FASTA_identity": float(fident),
                "query_align_seq": query_align_seq,
                "subject_align_seq": subject_align_seq,
                "match_markup_seq": markup,
                # a3m coordinates are 0-based and inclusive; BLAST's are 1-based.
                "query_start": int(qstart) + 1,
                "query_end": int(qend) + 1,
                "subject_start": int(tstart) + 1,
                "subject_end": int(tend) + 1,
            }
        )

    return rows, query_seq, n_excluded


def parse_a3m_to_csv(
    acc: str,
    a3m_tar_gz: Path,
    blast_csv_tar: Path,
    e_value_cutoff: float,
    use_env: bool,
    logging,
) -> tuple[str, bool, str]:
    """Convert one downloaded MSA archive into the homologue csv the pipeline reads.

    The query itself is written as hit 0 with a perfect self-alignment, because the NCBI path
    did the same (its first hit was the query matching itself in nr) and downstream code depends
    on it: ``extract_filtered_csv_homologues_to_alignments`` takes the length of row 0's TMD-plus-5
    slice as the ungapped reference length and drops every row that differs from it.
    """
    if not Path(a3m_tar_gz).is_file():
        warning = f"{acc} parse_a3m_to_csv failed, a3m archive not found = {a3m_tar_gz}"
        logging.warning(warning)
        return acc, False, warning

    with tarfile.open(a3m_tar_gz, "r:gz") as tar:
        names = tar.getnames()
        wanted = [UNIREF_A3M] + ([ENV_A3M] if use_env else [])
        a3m_texts = {}
        for name in wanted:
            if name not in names:
                if name == UNIREF_A3M:
                    raise FileNotFoundError(f"{acc} {a3m_tar_gz} has no {UNIREF_A3M}, only {names}")
                logging.warning(f"{acc} no {name} in the archive, continuing with UniRef hits only")
                continue
            extracted = tar.extractfile(name)
            if extracted is None:
                raise FileNotFoundError(f"{name} is not a regular file in {a3m_tar_gz}")
            with extracted as handle:
                a3m_texts[name] = handle.read().decode()

    all_rows: list[dict] = []
    query_seq = ""
    n_excluded_total = 0

    for name, text in a3m_texts.items():
        source = "uniref" if name == UNIREF_A3M else "env"
        rows, this_query_seq, n_excluded = _rows_from_a3m(text, source, e_value_cutoff, logging)
        if query_seq and this_query_seq != query_seq:
            raise ValueError(f"{acc} the a3m files in {a3m_tar_gz} were built from different queries")
        query_seq = this_query_seq
        all_rows.extend(rows)
        n_excluded_total += n_excluded

    # The two databases are searched separately and can both return the same UniRef entry, so a
    # merged alignment has to be deduplicated. Keeping the first occurrence keeps the UniRef hit
    # over the environmental one, which is the better-annotated of the two.
    seen: set[str] = set()
    deduplicated: list[dict] = []
    for row in all_rows:
        key = row["subject_align_seq"]
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(row)

    n_duplicates = len(all_rows) - len(deduplicated)

    # Best hit first, matching BLAST's ordering, so that any downstream code that treats earlier
    # rows as better behaves the way it did on the NCBI path.
    deduplicated.sort(key=lambda row: row["FASTA_expectation"])

    query_row = {
        "description": f"{acc}_colabfold_query_sequence",
        "organism": "no_organism",
        "FASTA_expectation": 0.0,
        "FASTA_identity": 1.0,
        "query_align_seq": query_seq,
        "subject_align_seq": query_seq,
        "match_markup_seq": query_seq,
        "query_start": 1,
        "query_end": len(query_seq),
        "subject_start": 1,
        "subject_end": len(query_seq),
    }

    rows_out = [query_row] + deduplicated
    for hit_num, row in enumerate(rows_out):
        row["hit_num"] = hit_num

    blast_csv_file = Path(str(blast_csv_tar)[: -len(".tar.gz")])
    blast_csv_file.parent.mkdir(parents=True, exist_ok=True)

    with open(blast_csv_file, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CSV_COLUMNS,
            extrasaction="ignore",
            delimiter=",",
            quotechar='"',
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
            doublequote=True,
        )
        writer.writeheader()
        writer.writerows(rows_out)

    with tarfile.open(blast_csv_tar, mode="w:gz") as tar:
        tar.add(blast_csv_file, arcname=blast_csv_file.name)
    blast_csv_file.unlink()

    logging.info(
        f"{acc} parse_a3m_to_csv finished: {len(deduplicated)} homologues "
        f"({n_duplicates} duplicates merged, {n_excluded_total} excluded by e-value cutoff "
        f"{e_value_cutoff}). Output = {blast_csv_tar}"
    )
    return acc, True, "no errors"


def parse_a3m_to_csv_mult_prot(
    paths: ArtefactPaths,
    df_set,
    e_value_cutoff: float,
    use_env: bool,
    logging,
) -> None:
    """Run parse_a3m_to_csv for a set of proteins.

    Parameters
    ----------
    paths : ArtefactPaths
        Locations of the pipeline's input and output files.
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process.
    e_value_cutoff : float
        Hits above this expectation value are dropped, as on the NCBI path.
    use_env : bool
        Include the environmental database hits alongside the UniRef ones.
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    logging.info("~~~~~~~~~~~~   starting parse_a3m_to_csv_mult_prot   ~~~~~~~~~~~~")

    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        parse_a3m_to_csv(
            acc,
            paths.colabfold_a3m_tar(database, acc),
            paths.homologue_csv_tar(database, acc),
            e_value_cutoff,
            use_env,
            logging,
        )

    logging.info("~~~~~~~~~~~~   finished parse_a3m_to_csv_mult_prot   ~~~~~~~~~~~~")
