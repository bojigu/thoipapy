"""What does THOIPA gain by taking its homologues from the ColabFold MSA server?

    THOIPA_COLABFOLD_CONTACT_EMAIL=you@example.com python scripts/compare_colabfold_alignment_depth.py

Every number in the paper rests on alignments built by a single blastp pass against NCBI ``nr`` in
May 2020. Two things are wrong with that today. The archive is six years old, and a remote ``nr``
query now queues for hours, so it cannot be rerun on demand. The deeper problem is the search
itself: one blastp pass is the weakest link in the pipeline, because conservation and coevolution
features are computed from the alignment and cannot be better than it is.

The ColabFold MSA server runs MMseqs2 with three profile iterations against UniRef30, then expands
each matched cluster back to its members, and optionally repeats the search against an
environmental database assembled from BFD, MGnify, MetaEuk and SMAG. That is a different class of
search from a single blastp pass, and the gain is not subtle.

This script measures the gain on one axis only -- alignment depth, the number of unique TMD
homologues that survive the existing gap and identity filters. That is the raw input to every
conservation and coevolution feature. It deliberately does not retrain anything: whether depth
buys accuracy is a separate question, and answering it means rebuilding the whole feature set,
which ``compare_blast_database_depth.py`` already lays out the shape of.

The comparison is confounded with time and cannot be otherwise. The tracked alignments are from
2020 and a search run today sees six years of extra sequence, so the difference measured here is
"iterative profile search now" against "single blastp pass then", not a clean test of the search
alone. NCBI does not archive old ``nr`` releases, so the 2020 snapshot cannot be recovered and the
clean test cannot be run by anyone.

Downloads are written outside the repository and are deliberately not DVC-tracked. The published
nr artefacts stay the reference until something is shown to be better than them.
"""

import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT))

from thoipapy.artefacts import ArtefactPaths  # noqa: E402
from thoipapy.homologues.colabfold_download import download_homologues_from_colabfold  # noqa: E402
from thoipapy.homologues.colabfold_parser import parse_a3m_to_csv  # noqa: E402
from thoipapy.homologues.NCBI_parser import extract_filtered_csv_homologues_to_alignments  # noqa: E402

NR_DATA = REPO_ROOT / "data"
_OUTPUT_DIR_ENV_VAR = "THOIPA_COLABFOLD_DATA"

SETNAME = "set08"

# From thoipapy/setting/run_settings_example.csv. Hardcoded rather than read from the spreadsheet,
# because a comparison whose filters can drift is not a comparison.
E_VALUE_CUTOFF = 100
MIN_IDENTITY_OF_TMD_SEQ = 0.2
MAX_N_GAPS_IN_TMD_QUERY_SEQ = 0
MAX_N_GAPS_IN_TMD_SUBJECT_SEQ = 5
NUM_OF_SUR_RESIDUES = 20

COLABFOLD_MODE = "env"


def output_dir() -> Path:
    """Where to put the downloaded alignments, from the environment.

    No default. These are large, they are not the reference dataset, and their location is a
    property of the machine rather than of this repository.
    """
    value = os.environ.get(_OUTPUT_DIR_ENV_VAR, "").strip()
    if not value:
        raise SystemExit(
            f"{_OUTPUT_DIR_ENV_VAR} is not set. Point it at a directory outside the repository "
            f"for the downloaded alignments to be written to."
        )
    return Path(value)


def n_seqs_in_fasta(path: Path) -> int:
    """Count records in a fasta file, or 0 if it was never written."""
    if not path.is_file():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith(">"))


def main() -> None:
    log = logging.getLogger("compare_colabfold_alignment_depth")
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

    df_set = pd.read_csv(NR_DATA / "Input_data" / f"{SETNAME}_train_processed.csv")

    paths = ArtefactPaths(
        data_dir=output_dir(),
        setname=SETNAME,
        num_of_sur_residues=NUM_OF_SUR_RESIDUES,
        max_n_gaps_in_TMD_subject_seq=MAX_N_GAPS_IN_TMD_SUBJECT_SEQ,
        homologue_source="colabfold",
    )
    nr_paths = ArtefactPaths(
        data_dir=NR_DATA,
        setname=SETNAME,
        num_of_sur_residues=NUM_OF_SUR_RESIDUES,
        max_n_gaps_in_TMD_subject_seq=MAX_N_GAPS_IN_TMD_SUBJECT_SEQ,
    )
    suffix = paths.surr_gaps_suffix

    rows = []
    for i in df_set.index:
        row = df_set.loc[i]
        acc, database = row["acc"], row["database"]

        a3m_tar = paths.colabfold_a3m_tar(database, acc)
        details_txt = Path(str(a3m_tar)[: -len(".a3m.tar.gz")] + "_details.txt")

        started = time.time()
        if not a3m_tar.is_file():
            download_homologues_from_colabfold(acc, row["TMD_seq_pl_surr"], a3m_tar, details_txt, COLABFOLD_MODE, log)
        seconds = time.time() - started

        csv_tar = paths.homologue_csv_tar(database, acc)
        csv_tar.parent.mkdir(parents=True, exist_ok=True)
        parse_a3m_to_csv(acc, a3m_tar, csv_tar, E_VALUE_CUTOFF, COLABFOLD_MODE == "env", log)

        alignments_dir = paths.alignment_dir(database)
        alignments_dir.mkdir(parents=True, exist_ok=True)
        summary = extract_filtered_csv_homologues_to_alignments(
            MAX_N_GAPS_IN_TMD_QUERY_SEQ,
            MAX_N_GAPS_IN_TMD_SUBJECT_SEQ,
            MIN_IDENTITY_OF_TMD_SEQ,
            acc,
            row["TMD_len"],
            alignments_dir / f"{acc}.{suffix}.redundant.fas",
            alignments_dir / f"{acc}.{suffix}.uniq.for_PSSM_FREECONTACT.txt",
            alignments_dir / f"{acc}.{paths.surr_suffix}.gaps0.uniq.for_LIPS.txt",
            alignments_dir / f"{acc}.surr5.gaps{MAX_N_GAPS_IN_TMD_SUBJECT_SEQ}.uniq.for_LIPO.txt",
            csv_tar,
            row["TMD_seq"],
            row["TMD_seq_pl_surr5"],
            log,
        )

        tracked = nr_paths.alignment_dir(database) / f"{acc}.{suffix}.uniq.for_PSSM_FREECONTACT.fas"
        n_nr = n_seqs_in_fasta(tracked)
        n_colabfold = summary["n_uniq_TMD_seqs_for_PSSM_FREECONTACT"]

        rows.append(
            {
                "acc": acc,
                "database": database,
                "TMD_len": row["TMD_len"],
                "n_uniq_TMD_seqs_nr_2020": n_nr,
                "n_uniq_TMD_seqs_colabfold": n_colabfold,
                "n_hits_colabfold": summary["n_total_BLAST_hits"],
                "download_seconds": round(seconds, 1),
            }
        )
        print(f"{acc:10s} {database:8s} nr={n_nr:6d}  colabfold={n_colabfold:6d}")

    df = pd.DataFrame(rows)
    # Proteins whose tracked alignment is missing cannot be compared, only reported.
    comparable = df.loc[df["n_uniq_TMD_seqs_nr_2020"] > 0].copy()
    comparable["fold_deeper"] = comparable["n_uniq_TMD_seqs_colabfold"] / comparable["n_uniq_TMD_seqs_nr_2020"]

    out_csv = output_dir() / f"{SETNAME}_colabfold_vs_nr_alignment_depth.csv"
    df.to_csv(out_csv, index=False)

    print(f"\n{len(comparable)} of {len(df)} proteins have a tracked nr alignment to compare against")
    print(f"median fold deeper : {comparable['fold_deeper'].median():.2f}")
    print(f"range              : {comparable['fold_deeper'].min():.2f} to {comparable['fold_deeper'].max():.2f}")
    print("\nby database:")
    print(comparable.groupby("database")["fold_deeper"].agg(["count", "median"]).to_string())
    print(f"\nwritten to {out_csv}")


if __name__ == "__main__":
    main()
