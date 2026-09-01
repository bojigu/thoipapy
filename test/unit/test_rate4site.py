from pathlib import Path

import pandas as pd
import pytest
from thoipapy.features.rate4site import rate4site_calculation, read_fasta_records
from thoipapy.utils import LogOnlyToConsole, SurroundingSequence, safe_filename_component

from test.helpers.helpers import TestProtein

# The committed alignment's query defline is ">1xioA4_01-2211oclock_orig_seq". The hyphens in it
# are why this fixture reproduces the live bug: the gap-stripping that built the cd-hit input used
# to run over the deflines too, so the query came back from cd-hit under a different name, failed
# the membership test, and was dropped from the alignment that rate4site scored.
ALIGNMENT_ACC = "1xioA4_01-2211oclock"


COMMITTED_ALIGNMENT = Path(__file__).parents[1] / "test_inputs/tm_homologues/homologues_uniq_for_pssm_freecontact.fas"


def alignment_with_query_named(acc: str, out_dir: Path) -> Path:
    """Copy the committed alignment, renaming its query defline to match the given accession.

    The alignment is written by NCBI_parser.save_fasta_from_array, which names the query record
    after the accession, so the two always agree in the pipeline. A test that varies the accession
    has to vary the alignment with it.
    """
    renamed = out_dir / f"{acc}.for_LIPO.fas"
    renamed.write_text(
        "".join(
            f">{acc}_orig_seq\n" if line.startswith(">") and line.rstrip().endswith("_orig_seq") else line
            for line in COMMITTED_ALIGNMENT.read_text().splitlines(keepends=True)
        )
    )
    return renamed


def run_rate4site_on_the_committed_alignment(acc: str, out_dir: Path) -> pd.DataFrame:
    """Run rate4site over the committed homologue alignment and return the parsed features."""
    fasta_uniq_TMD_seqs_surr5_for_LIPO = alignment_with_query_named(acc, out_dir)
    rate4site_csv = out_dir / "test_rate4site.csv"
    rate4site_csv.parent.mkdir(parents=True, exist_ok=True)
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    surrounding_sequence = SurroundingSequence(tp.tmd_start, tp.tmd_end, len(tp.full_seq), num_of_sur_residues=5)

    rate4site_calculation(
        tp.tmd_seq,
        acc,
        fasta_uniq_TMD_seqs_surr5_for_LIPO,
        rate4site_csv,
        surrounding_sequence.n_term_offset,
        LogOnlyToConsole(),
    )

    assert rate4site_csv.is_file()
    return pd.read_csv(rate4site_csv)


@pytest.mark.requires_tool("rate4site", "cd-hit")
def test_rate4site_scores_the_submitted_tmd(tmp_path):
    """The residue names must be the query's, not those of whichever homologue headed the alignment.

    This used to assert "GFLMSTQVVVITTGLIADL" -- homologue ">1", which differs from the query at
    positions 8 and 13 -- and treated that as expected behaviour.
    """
    tp = TestProtein()
    tp.with_1xioA4()
    df = run_rate4site_on_the_committed_alignment(ALIGNMENT_ACC, tmp_path)

    assert "rate4site" in df.columns
    assert "".join(df["residue_name"].to_list()) == tp.tmd_seq
    assert list(df["residue_num"]) == list(range(1, len(tp.tmd_seq) + 1))
    assert df["rate4site"].notna().all()


@pytest.mark.requires_tool("rate4site", "cd-hit")
def test_the_query_heads_the_alignment_handed_to_rate4site(tmp_path):
    """The invariant the fix establishes, checked on the file rate4site actually reads."""
    run_rate4site_on_the_committed_alignment(ALIGNMENT_ACC, tmp_path)

    records = read_fasta_records(tmp_path / f"{ALIGNMENT_ACC}.rate4site_input.fas")
    assert records[0][0] == f"{ALIGNMENT_ACC}_orig_seq"


@pytest.mark.requires_tool("rate4site", "cd-hit")
def test_rate4site_runs_with_a_name_derived_from_a_uniprot_fasta_header(tmp_path):
    """A full UniProt header is the most likely thing a user pastes, and it triggered both bugs.

    Sanitised it still contains hyphens ("Glycophorin-A"), which is the case that displaced the
    rate4site reference sequence.
    """
    tp = TestProtein()
    tp.with_1xioA4()
    acc = safe_filename_component("sp|P02724|GLPA_HUMAN Glycophorin-A OS=Homo sapiens OX=9606")
    assert "-" in acc, "the regression case needs a hyphen in the accession"

    df = run_rate4site_on_the_committed_alignment(acc, tmp_path)

    assert "".join(df["residue_name"].to_list()) == tp.tmd_seq


@pytest.mark.requires_tool("rate4site", "cd-hit")
def test_rate4site_leaves_no_temporary_files_in_the_working_directory(tmp_path):
    """rate4site writes r4s.res, r4sOrig.res and TheTree.txt into the working directory."""
    working_directory_before = {p.name for p in Path.cwd().iterdir()}

    run_rate4site_on_the_committed_alignment(ALIGNMENT_ACC, tmp_path)

    for temporary in ["r4s.res", "r4sOrig.res", "TheTree.txt"]:
        assert not (tmp_path / temporary).exists(), f"{temporary} was left in the output directory"
        if temporary not in working_directory_before:
            assert not (Path.cwd() / temporary).exists(), f"{temporary} was written to the working directory"
