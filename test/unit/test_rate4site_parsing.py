"""Tests for the rate4site output parser and the alignment it is given.

These need no external binaries, so they cover the fix on a clean checkout. Every other rate4site
test is marked requires_tool and skips when rate4site and cd-hit are absent, which is what let the
original defect survive: the suite passed, having tested nothing.
"""

from pathlib import Path

import pytest
from thoipapy.features.rate4site import (
    count_rate4site_header_lines,
    parse_rate4site_output,
    query_record_from_alignment,
    read_fasta_records,
    truncate_alignment,
    write_cdhit_input,
    write_rate4site_input,
)

FIXTURES = Path(__file__).parents[1] / "test_inputs/rate4site"


@pytest.fixture
def rate4site_output(tmp_path):
    """Copy a committed rate4site output into a temporary directory before parsing it.

    parse_rate4site_output writes a .orig.csv beside its input, so parsing a fixture in place
    would litter test_inputs/.
    """

    def copy(name: str) -> Path:
        destination = tmp_path / name
        destination.write_text((FIXTURES / name).read_text())
        return destination

    return copy


# The 1xioA4 test protein: TMD as submitted, and the same TMD plus 5 residues each side, which is
# what the LIPO alignment holds and what rate4site scores.
TMD_SEQ = "GFLMSTQIVVITSGLIADL"
QUERY_SURR5 = "DWTLI" + TMD_SEQ + "SERDW"
N_TERM_OFFSET = 5
ACC = "1xioA4"


def test_the_header_length_is_measured_not_assumed(rate4site_output):
    """It was hardcoded to 13. A banner one line off would eat the first data row."""
    assert count_rate4site_header_lines(rate4site_output("matching_reference.rate4site_orig_output.txt")) == 13


def test_residue_names_come_from_the_query(rate4site_output):
    df = parse_rate4site_output(
        rate4site_output("matching_reference.rate4site_orig_output.txt"), TMD_SEQ, N_TERM_OFFSET, ACC
    )

    assert "".join(df["residue_name"].tolist()) == TMD_SEQ
    assert list(df.index) == list(range(1, len(TMD_SEQ) + 1))
    assert df["rate4site"].notna().all()


def test_a_reference_that_is_not_the_query_is_refused(rate4site_output):
    """The live bug. rate4site scored a homologue, so its scores belong to another sequence.

    Reporting them under the query's residue names would be a confidently wrong feature; the merge
    in combine_all_features used to drop the offending position and complain about the length two
    stages later.
    """
    with pytest.raises(ValueError) as excinfo:
        parse_rate4site_output(
            rate4site_output("diverging_reference.rate4site_orig_output.txt"), TMD_SEQ, N_TERM_OFFSET, ACC
        )

    message = str(excinfo.value)
    assert "pos 8: query I, rate4site V" in message
    assert "pos 13: query S, rate4site T" in message


def test_a_short_rate4site_output_is_refused_rather_than_becoming_nan(rate4site_output):
    with pytest.raises(ValueError, match="produced no score for TMD residue"):
        parse_rate4site_output(rate4site_output("truncated.rate4site_orig_output.txt"), TMD_SEQ, N_TERM_OFFSET, ACC)


def write_alignment(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{header}\n{seq}\n" for header, seq in records))
    return path


def test_the_query_record_is_identified_and_checked(tmp_path):
    records = [(f"{ACC}_orig_seq", QUERY_SURR5), ("0", QUERY_SURR5), ("1", QUERY_SURR5)]
    header, seq = query_record_from_alignment(records, ACC, TMD_SEQ, N_TERM_OFFSET, "alignment.fas")

    assert header == f"{ACC}_orig_seq"
    assert seq == QUERY_SURR5


def test_an_alignment_whose_first_record_is_not_the_query_is_refused():
    records = [("0", QUERY_SURR5), (f"{ACC}_orig_seq", QUERY_SURR5)]
    with pytest.raises(ValueError, match="is '0', not the query"):
        query_record_from_alignment(records, ACC, TMD_SEQ, N_TERM_OFFSET, "alignment.fas")


def test_a_gapped_query_record_is_refused():
    """rate4site numbers its output by ungapped position in the reference."""
    gapped = QUERY_SURR5[:10] + "-" + QUERY_SURR5[11:]
    with pytest.raises(ValueError, match="contains gaps"):
        query_record_from_alignment([(f"{ACC}_orig_seq", gapped)], ACC, TMD_SEQ, N_TERM_OFFSET, "alignment.fas")


def test_a_query_record_whose_tmd_is_not_at_the_expected_offset_is_refused():
    """The assumption the whole index shift rests on, and which was never checked."""
    with pytest.raises(ValueError, match="not at offset 5"):
        query_record_from_alignment(
            [(f"{ACC}_orig_seq", "AAA" + TMD_SEQ + "SERDW")], ACC, TMD_SEQ, N_TERM_OFFSET, "alignment.fas"
        )


def test_the_query_heads_the_rate4site_input_even_when_cdhit_drops_it(tmp_path):
    """cd-hit removes redundancy, and the query is exactly the kind of sequence it removes.

    Its absence used to be invisible: rate4site simply scored whatever record was left first.
    """
    records = [(f"{ACC}_orig_seq", QUERY_SURR5), ("0", "A" * len(QUERY_SURR5)), ("1", "C" * len(QUERY_SURR5))]
    rate4site_input = tmp_path / "rate4site_input.fas"

    write_rate4site_input(records, f"{ACC}_orig_seq", cdhit_cluster_reps=["0", "1"], rate4site_input=rate4site_input)
    written = read_fasta_records(rate4site_input)

    assert written[0] == (f"{ACC}_orig_seq", QUERY_SURR5)
    assert [header for header, _ in written] == [f"{ACC}_orig_seq", "0", "1"]


def test_the_query_is_not_written_twice_when_cdhit_keeps_it(tmp_path):
    records = [(f"{ACC}_orig_seq", QUERY_SURR5), ("0", "A" * len(QUERY_SURR5))]
    rate4site_input = tmp_path / "rate4site_input.fas"

    write_rate4site_input(
        records, f"{ACC}_orig_seq", cdhit_cluster_reps=[f"{ACC}_orig_seq", "0"], rate4site_input=rate4site_input
    )

    assert [header for header, _ in read_fasta_records(rate4site_input)] == [f"{ACC}_orig_seq", "0"]


def test_gap_stripping_does_not_touch_the_deflines(tmp_path):
    """The defect itself: a hyphen in the accession was removed from the header as well.

    The query is then the one record whose name no longer matches on the way back from cd-hit, so
    it is dropped from the alignment and rate4site scores a homologue instead.
    """
    hyphenated = "Serine-rich_protein"
    records = [(f"{hyphenated}_orig_seq", QUERY_SURR5), ("0", "A" * 10 + "-" * 5 + "A" * 14)]
    cdhit_input = tmp_path / "cdhit_input.fas"

    write_cdhit_input(records, cdhit_input)
    written = read_fasta_records(cdhit_input)

    assert written[0][0] == f"{hyphenated}_orig_seq", "the defline must keep its hyphens"
    assert written[0][1] == QUERY_SURR5, "the query is gap-free and must be unchanged"
    assert written[1][1] == "A" * 24, "gaps must still be stripped from the sequences"


def test_truncating_an_alignment_keeps_whole_records(tmp_path):
    """Truncation used to cut at a fixed line count, which could leave a defline with no sequence."""
    alignment = tmp_path / "cdhit_output.fas"
    records = [(str(n), f"{'ACDEF' * 5}{n % 10}") for n in range(50)]
    alignment.write_text("".join(f">{header}\n{seq}\n" for header, seq in records))

    kept = truncate_alignment(alignment, max_n_sequences=7)

    assert kept == [str(n) for n in range(7)]
    written = read_fasta_records(alignment)
    assert len(written) == 7
    assert all(seq for _header, seq in written), "every record kept must have its sequence"


def test_an_accession_rate4site_would_read_as_an_option_is_refused():
    with pytest.raises(ValueError, match="command-line option"):
        query_record_from_alignment([("-a_orig_seq", QUERY_SURR5)], "-a", TMD_SEQ, N_TERM_OFFSET, "alignment.fas")
