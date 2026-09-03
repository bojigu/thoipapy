"""Unit tests for the a3m -> homologue csv conversion.

The reconstruction of a pairwise alignment from an a3m record is the only place in the ColabFold
path where THOIPA does alignment arithmetic of its own, so it is the part worth pinning down. The
fixtures below are hand-built rather than server-captured, so that each one isolates a single
property: an exact hit, an insertion, partial coverage, and the trailing NUL byte MMseqs2 leaves
on the last record of an a3m.
"""

import tarfile
from pathlib import Path

import pandas as pd
import pytest
from thoipapy.homologues.colabfold_parser import (
    pairwise_alignment_from_a3m_record,
    parse_a3m_to_csv,
    read_a3m,
)
from thoipapy.utils import LogOnlyToConsole

QUERY = "AAACCCGGG"


def test_read_a3m_returns_query_first():
    records = read_a3m(f">101\n{QUERY}\n>hit1\tx\nAAACCCGGT\n")
    assert records[0] == ("101", QUERY)
    assert records[1] == ("hit1\tx", "AAACCCGGT")


def test_read_a3m_strips_the_trailing_nul_mmseqs2_leaves_behind():
    """The last record of an a3m written by result2msa carries MMseqs2's NUL terminator.

    Left in, it makes that record one column longer than the query and the reconstruction walks
    off the end of the query sequence.
    """
    records = read_a3m(f">101\n{QUERY}\n>hit1\tx\nAAACCCGGT\x00\n")
    assert records[-1][1] == "AAACCCGGT"


def test_read_a3m_joins_wrapped_sequence_lines():
    records = read_a3m(f">101\n{QUERY}\n>hit1\tx\nAAACC\nCGGT\n")
    assert records[1][1] == "AAACCCGGT"


def test_exact_hit_gives_identical_alignment_and_full_markup():
    query_aln, subject_aln, markup = pairwise_alignment_from_a3m_record(QUERY, QUERY)
    assert query_aln == QUERY
    assert subject_aln == QUERY
    assert markup == QUERY


def test_mismatch_is_a_space_in_the_markup():
    _, _, markup = pairwise_alignment_from_a3m_record(QUERY, "AAACCCGGT")
    assert markup == "AAACCCGG "


def test_lowercase_is_an_insertion_that_gaps_the_query():
    """A lowercase column is a residue the subject has and the query does not."""
    query_aln, subject_aln, markup = pairwise_alignment_from_a3m_record(QUERY, "AAAcccCCCGGG")

    assert query_aln == "AAA---CCCGGG"
    assert subject_aln == "AAACCCCCCGGG"
    # An insertion is never an identity, so the markup is blank across it.
    assert markup == "AAA   CCCGGG"
    assert len(query_aln) == len(subject_aln) == len(markup)


def test_dash_is_a_query_residue_the_subject_does_not_cover():
    query_aln, subject_aln, markup = pairwise_alignment_from_a3m_record(QUERY, "---CCCGGG")
    assert query_aln == QUERY
    assert subject_aln == "---CCCGGG"
    assert markup == "   CCCGGG"


def test_record_with_too_many_match_states_is_rejected():
    with pytest.raises(ValueError, match="more match-state columns"):
        pairwise_alignment_from_a3m_record(QUERY, QUERY + "A")


def test_record_with_too_few_match_states_is_rejected():
    with pytest.raises(ValueError, match="of the query's"):
        pairwise_alignment_from_a3m_record(QUERY, "AAACCC")


def _write_a3m_archive(tmp_path: Path, uniref_text: str, env_text: str | None = None) -> Path:
    """Build an archive shaped like the one the ColabFold MSA server returns."""
    archive = tmp_path / "query.a3m.tar.gz"
    (tmp_path / "uniref.a3m").write_text(uniref_text)
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(tmp_path / "uniref.a3m", arcname="uniref.a3m")
        if env_text is not None:
            (tmp_path / "env.a3m").write_text(env_text)
            tar.add(tmp_path / "env.a3m", arcname="bfd.mgnify30.metaeuk30.smag30.a3m")
    return archive


def _header(target, fident, evalue, qstart, qend, tstart=1, tend=9, tlen=50):
    return f"{target}\t42\t{fident}\t{evalue}\t{qstart}\t{qend}\t{len(QUERY)}\t{tstart}\t{tend}\t{tlen}"


def _read_csv_from_tar(csv_tar: Path) -> pd.DataFrame:
    with tarfile.open(csv_tar, "r:gz") as tar:
        member = tar.getnames()[0]
        extracted = tar.extractfile(member)
        assert extracted is not None
        with extracted as handle:
            return pd.read_csv(handle)


def test_parse_a3m_to_csv_writes_the_query_as_hit_zero(tmp_path):
    """Downstream takes row 0 as the ungapped reference length, so it must be the query."""
    a3m = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
    archive = _write_a3m_archive(tmp_path, a3m)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=False, logging=LogOnlyToConsole())

    df = _read_csv_from_tar(csv_tar)
    assert df.loc[0, "hit_num"] == 0
    assert df.loc[0, "description"] == "testacc_colabfold_query_sequence"
    assert df.loc[0, "subject_align_seq"] == QUERY
    assert df.loc[0, "query_align_seq"] == QUERY
    assert df.loc[1, "description"] == "UniRef100_X uniref"


def test_parse_a3m_to_csv_produces_the_columns_the_ncbi_path_produces(tmp_path):
    """The csv is read by code written for BLAST output, so the two must be interchangeable."""
    a3m = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
    archive = _write_a3m_archive(tmp_path, a3m)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=False, logging=LogOnlyToConsole())

    df = _read_csv_from_tar(csv_tar)
    assert set(df.columns) == {
        "FASTA_expectation",
        "FASTA_identity",
        "description",
        "hit_num",
        "match_markup_seq",
        "organism",
        "query_align_seq",
        "query_end",
        "query_start",
        "subject_align_seq",
        "subject_end",
        "subject_start",
    }


def test_parse_a3m_to_csv_converts_coordinates_to_one_based(tmp_path):
    """a3m coordinates are 0-based and inclusive; the csv is read as BLAST's 1-based."""
    a3m = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.667, '1.0E-10', 3, 8, tstart=10, tend=15)}\n---CCCGGG\n"
    archive = _write_a3m_archive(tmp_path, a3m)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=False, logging=LogOnlyToConsole())

    hit = _read_csv_from_tar(csv_tar).iloc[1]
    assert hit["query_start"] == 4
    assert hit["query_end"] == 9
    assert hit["subject_start"] == 11
    assert hit["subject_end"] == 16


def test_parse_a3m_to_csv_drops_hits_above_the_e_value_cutoff(tmp_path):
    a3m = (
        f">101\n{QUERY}\n"
        f">{_header('UniRef100_GOOD', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
        f">{_header('UniRef100_BAD', 0.333, '9.0E+00', 0, 8)}\nAAAGGGTTT\n"
    )
    archive = _write_a3m_archive(tmp_path, a3m)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=1.0, use_env=False, logging=LogOnlyToConsole())

    descriptions = _read_csv_from_tar(csv_tar)["description"].tolist()
    assert "UniRef100_GOOD uniref" in descriptions
    assert "UniRef100_BAD uniref" not in descriptions


def test_parse_a3m_to_csv_merges_the_environmental_database(tmp_path):
    uniref = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
    env = f">101\n{QUERY}\n>{_header('ENV_Y', 0.778, '1.0E-05', 0, 8)}\nAAACCCTTT\n"
    archive = _write_a3m_archive(tmp_path, uniref, env)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=True, logging=LogOnlyToConsole())

    descriptions = _read_csv_from_tar(csv_tar)["description"].tolist()
    assert "UniRef100_X uniref" in descriptions
    assert "ENV_Y env" in descriptions


def test_parse_a3m_to_csv_ignores_the_environmental_database_when_asked_to(tmp_path):
    uniref = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
    env = f">101\n{QUERY}\n>{_header('ENV_Y', 0.778, '1.0E-05', 0, 8)}\nAAACCCTTT\n"
    archive = _write_a3m_archive(tmp_path, uniref, env)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=False, logging=LogOnlyToConsole())

    descriptions = _read_csv_from_tar(csv_tar)["description"].tolist()
    assert "ENV_Y env" not in descriptions


def test_parse_a3m_to_csv_deduplicates_hits_found_in_both_databases(tmp_path):
    """The two databases are searched separately and can return the same sequence."""
    uniref = f">101\n{QUERY}\n>{_header('UniRef100_X', 0.889, '1.0E-10', 0, 8)}\nAAACCCGGT\n"
    env = f">101\n{QUERY}\n>{_header('ENV_SAME', 0.889, '1.0E-08', 0, 8)}\nAAACCCGGT\n"
    archive = _write_a3m_archive(tmp_path, uniref, env)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=True, logging=LogOnlyToConsole())

    df = _read_csv_from_tar(csv_tar)
    # the query, plus one of the two identical subject sequences
    assert len(df) == 2
    assert df["description"].tolist() == ["testacc_colabfold_query_sequence", "UniRef100_X uniref"]


def test_parse_a3m_to_csv_orders_hits_by_e_value(tmp_path):
    a3m = (
        f">101\n{QUERY}\n"
        f">{_header('UniRef100_WEAK', 0.444, '1.0E-02', 0, 8)}\nAAACCCTTT\n"
        f">{_header('UniRef100_STRONG', 0.889, '1.0E-20', 0, 8)}\nAAACCCGGT\n"
    )
    archive = _write_a3m_archive(tmp_path, a3m)
    csv_tar = tmp_path / "out.csv.tar.gz"

    parse_a3m_to_csv("testacc", archive, csv_tar, e_value_cutoff=100, use_env=False, logging=LogOnlyToConsole())

    descriptions = _read_csv_from_tar(csv_tar)["description"].tolist()
    assert descriptions[1] == "UniRef100_STRONG uniref"
    assert descriptions[2] == "UniRef100_WEAK uniref"


def test_parse_a3m_to_csv_reports_a_missing_archive_rather_than_raising(tmp_path):
    acc, worked, message = parse_a3m_to_csv(
        "testacc", tmp_path / "absent.a3m.tar.gz", tmp_path / "out.csv.tar.gz", 100, False, LogOnlyToConsole()
    )
    assert acc == "testacc"
    assert worked is False
    assert "not found" in message


def test_parse_a3m_to_csv_rejects_an_archive_without_uniref(tmp_path):
    archive = tmp_path / "query.a3m.tar.gz"
    (tmp_path / "other.a3m").write_text(f">101\n{QUERY}\n")
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(tmp_path / "other.a3m", arcname="something_else.a3m")

    with pytest.raises(FileNotFoundError, match="uniref.a3m"):
        parse_a3m_to_csv("testacc", archive, tmp_path / "out.csv.tar.gz", 100, False, LogOnlyToConsole())
