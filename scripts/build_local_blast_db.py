"""Download a UniProt release and build a local BLAST database from it.

Predictions query NCBI by default, which is fine until NCBI queues you. A probe of a
three-residue query has returned RTOE=11165 -- their own estimate of 3.1 hours before results are
ready. A webserver job with a one-hour ceiling cannot survive that, and a test suite must not
depend on it at all. Point THOIPA_LOCAL_BLAST_DB at the output of this script and the
predictor searches locally instead.

    python scripts/build_local_blast_db.py

Two databases are worth building, for different jobs:

swissprot   575k sequences, ~90 MB download, ~340 MB built, seconds to build.
            Use for tests, CI and development. NOT for production: it yields roughly a quarter
            of the alignment depth of nr (26 hits vs 76 unique TMD homologues for 1xioA4,
            18 vs 101 for 4hksA1), and conservation and coevolution features degrade with
            shallow alignments.

uniref90    121M clusters, ~30 GB download, expect 60-80 GB built.
            The production option. UniRef clusters at 90% identity, which is close to what the
            pipeline discards anyway -- it collapses exact-duplicate TMD sequences and runs
            CD-HIT -- so little of value is lost.

Deliberately not filtered to membrane proteins. THOIPA needs homologues of the query protein
itself, so what matters is taxonomic depth within the query's own family, not enrichment for
membrane proteins. Filtering on a transmembrane annotation would also silently drop genuine
homologues that simply lack the annotation, and alignment depth is the thing you cannot afford
to lose.
"""

import subprocess
import sys
import urllib.request
from pathlib import Path

BLAST_DB_DIR = Path("/mnt/shared/blastdb")

DATABASES = {
    "swissprot": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz",
        "fasta": "uniprot_sprot.fasta",
        "title": "UniProt Swiss-Prot",
    },
    "uniref90": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz",
        "fasta": "uniref90.fasta",
        "title": "UniProt UniRef90",
    },
}

# Which to build when the script is run. uniref90 is the default because it is what a real
# server needs; switch to "swissprot" for a fast, small database for tests and development.
DATABASE_TO_BUILD = "uniref90"


def download(url: str, destination: Path) -> None:
    """Fetch url to destination, skipping the download if the file is already there."""
    if destination.is_file():
        print(f"  already downloaded: {destination}")
        return
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, destination)
    print(f"  wrote {destination} ({destination.stat().st_size / 1e6:.0f} MB)")


def decompress(gzipped: Path, destination: Path) -> None:
    """Expand a .gz file, skipping if the output is already there."""
    if destination.is_file():
        print(f"  already decompressed: {destination}")
        return
    print(f"  decompressing {gzipped}")
    import gzip
    import shutil

    with gzip.open(gzipped, "rb") as src, open(destination, "wb") as dst:
        shutil.copyfileobj(src, dst)


def build_blast_database(fasta: Path, db_name: str, title: str) -> None:
    """Run makeblastdb over a fasta file."""
    print(f"  running makeblastdb -> {db_name}")
    result = subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "prot",
            "-out",
            str(BLAST_DB_DIR / db_name),
            "-title",
            title,
            "-parse_seqids",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {result.stderr.strip()}")
    print(result.stdout.strip())


def main() -> None:
    """Download and build the configured database."""
    if not shutil_which("makeblastdb"):
        sys.exit("makeblastdb not found. Install it with: mamba install -c bioconda blast")

    spec = DATABASES[DATABASE_TO_BUILD]
    BLAST_DB_DIR.mkdir(parents=True, exist_ok=True)

    gzipped = BLAST_DB_DIR / Path(spec["url"]).name
    fasta = BLAST_DB_DIR / spec["fasta"]

    download(spec["url"], gzipped)
    decompress(gzipped, fasta)
    build_blast_database(fasta, DATABASE_TO_BUILD, spec["title"])

    print()
    print("Done. To use it:")
    print(f"  export THOIPA_LOCAL_BLAST_DB={BLAST_DB_DIR / DATABASE_TO_BUILD}")


def shutil_which(name: str):
    """Return the path to an executable, or None."""
    from shutil import which

    return which(name)


if __name__ == "__main__":
    main()
