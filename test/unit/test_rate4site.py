import pandas as pd
from pathlib import Path
from test.helpers.helpers import TestProtein
from thoipapy.features.rate4site import rate4site_calculation
from thoipapy.utils import LogOnlyToConsole, SurroundingSequence


def test_rate4site():
    fasta_uniq_TMD_seqs_surr5_for_LIPO = Path(__file__).parents[1] / "test_inputs/tm_homologues/homologues_uniq_for_pssm_freecontact.fas"
    rate4site_csv = Path(__file__).parents[1] / "test_outputs/test_rate4site/test_rate4site.csv"
    if not rate4site_csv.parent.is_dir():
        rate4site_csv.parent.mkdir(parents=True)
    tp: TestProtein = TestProtein()
    tp.with_1xioA4()
    surrounding_sequence = SurroundingSequence(tp.tmd_start, tp.tmd_end, len(tp.full_seq), num_of_sur_residues=5)
    logging = LogOnlyToConsole()

    rate4site_calculation(tp.tmd_seq, tp.acc, fasta_uniq_TMD_seqs_surr5_for_LIPO, rate4site_csv, surrounding_sequence.n_term_offset, logging)

    assert rate4site_csv.is_file()
    with open(rate4site_csv, "r") as f:
        lines = f.readlines()
    assert "rate4site" in lines[0]
    df = pd.read_csv(rate4site_csv)
    seq_in_rate4site_output: str = "".join(df["residue_name"].to_list())
    assert seq_in_rate4site_output == "GFLMSTQVVVITTGLIADL"
    # note that the sequence in the rate4site output does not match the canonical TMD sequence
    assert tp.tmd_seq is not seq_in_rate4site_output

