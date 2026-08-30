"""Pin ArtefactPaths to the exact paths the pre-refactor code produced.

Each expected value below was taken from the inline path expression it replaces, so migrating a
call site to ArtefactPaths is verifiably a no-op. The literals are written out longhand on purpose:
deriving them from the same code under test would prove nothing.
"""
from pathlib import Path

import pytest

from thoipapy.artefacts import ArtefactPaths

DATA_DIR = Path("/tmp/thoipapy_test_data")
SETNAME = "set08"
DATABASE = "crystal"
ACC = "1xioA4"


@pytest.fixture
def paths():
    return ArtefactPaths(
        data_dir=DATA_DIR,
        setname=SETNAME,
        num_of_sur_residues=20,
        max_n_gaps_in_TMD_subject_seq=5,
    )


def test_suffixes(paths):
    assert paths.surr_gaps_suffix == "surr20.gaps5"
    assert paths.surr_suffix == "surr20"


@pytest.mark.parametrize("method,args,expected", [
    # per-protein features
    ("pssm_csv", (DATABASE, ACC), "features/pssm/crystal/1xioA4.surr20.gaps5.pssm.csv"),
    ("entropy_csv", (DATABASE, ACC), "features/entropy/crystal/1xioA4.surr20.gaps5.uniq.entropy.csv"),
    ("rate4site_csv", (DATABASE, ACC), "features/rate4site/crystal/1xioA4_rate4site.csv"),
    ("freecontact_csv", (DATABASE, ACC), "features/coevolution/crystal/1xioA4.surr20.gaps5.freecontact.csv"),
    ("freecontact_parsed_csv", (DATABASE, ACC), "features/coevolution/crystal/1xioA4.surr20.gaps5.freecontact_parsed.csv"),
    ("combined_features_csv", (DATABASE, ACC), "features/combined/crystal/1xioA4.surr20.gaps5.combined_features.csv"),
    ("homologue_xml_tar", (DATABASE, ACC), "homologues/xml/crystal/1xioA4.surr20.BLAST.xml.tar.gz"),
    ("homologue_csv_tar", (DATABASE, ACC), "homologues/ncbi/crystal/1xioA4.surr20.BLAST.csv.tar.gz"),
    ("alignment_dir", (DATABASE,), "homologues/alignments/crystal"),
    ("structure_dir", (DATABASE,), "features/structure/crystal"),
    ("other_predictors_dir", (DATABASE,), "Predictions/other_predictors/crystal"),
    ("protein_dir", (DATABASE,), "Proteins/crystal"),
    # training table
    ("train_data_orig_csv", (), "results/set08/train_data/01_train_data_orig.csv"),
    ("feat_imp_MDI_before_feature_seln_csv", (), "results/set08/train_data/01_feat_imp_MDI_before_feature_seln.csv"),
    ("tuned_ensemble_parameters_before_feature_seln_csv", (), "results/set08/train_data/01_tuned_ensemble_parameters_before_feature_seln.csv"),
    ("train_data_excl_duplicates_csv", (), "results/set08/train_data/02_train_data_excl_duplicates.csv"),
    ("train_data_after_first_feature_seln_csv", (), "results/set08/train_data/03_train_data_after_first_feature_seln.csv"),
    ("tuned_ensemble_parameters_csv", (), "results/set08/train_data/04_tuned_ensemble_parameters.csv"),
    ("hetero_contact_residues_csv", (), "results/set08/train_data/00_hetero_contact_residues.csv"),
    ("top_features_anova_csv", (), "results/set08/feat_imp/top_features_anova.csv"),
    ("top_features_rfe_csv", (), "results/set08/feat_imp/top_features_rfe.csv"),
    # model and clusters
    ("ml_model_lpkl", (), "results/set08/set08_ML_model.lpkl"),
    ("protein_set_full_seq_fasta", (), "results/set08/clusters/set08_full_seqs.fas"),
    ("protein_set_tmd_seq_fasta", (), "results/set08/clusters/set08_tmd_seqs.fas"),
    ("cdhit_cluster_txt", (), "results/set08/clusters/set08.fas.1.clstr.sorted.txt"),
    # input tables
    ("processed_set_csv", ("set08_train",), "Input_data/set08_train_processed.csv"),
    ("logging_dir", (), "Logging"),
])
def test_path_matches_pre_refactor_spelling(paths, method, args, expected):
    actual = getattr(paths, method)(*args)
    assert actual == DATA_DIR / expected, f"{method}{args} produced {actual}"


def test_suffix_follows_the_settings_rather_than_being_hardcoded():
    """Changing the settings must change the filenames, for readers as well as writers.

    Five reader sites used to hardcode surr20.gaps5 while every writer built it from the settings,
    so a non-default setting made the pipeline write one filename and read another.
    """
    non_default = ArtefactPaths(
        data_dir=DATA_DIR, setname=SETNAME,
        num_of_sur_residues=30, max_n_gaps_in_TMD_subject_seq=2,
    )
    assert non_default.surr_gaps_suffix == "surr30.gaps2"
    assert non_default.combined_features_csv(DATABASE, ACC).name == "1xioA4.surr30.gaps2.combined_features.csv"
    assert non_default.homologue_xml_tar(DATABASE, ACC).name == "1xioA4.surr30.BLAST.xml.tar.gz"


def test_from_settings_reads_the_expected_keys():
    s = {
        "data_dir": DATA_DIR,
        "setname": "set07",
        "num_of_sur_residues": 20,
        "max_n_gaps_in_TMD_subject_seq": 5,
    }
    paths = ArtefactPaths.from_settings(s)
    assert paths.setname == "set07"
    assert paths.train_data_orig_csv() == DATA_DIR / "results/set07/train_data/01_train_data_orig.csv"


def test_is_frozen(paths):
    """A mutable path layer is how two spellings of the same file appeared in the first place."""
    with pytest.raises(Exception):
        paths.setname = "set07"
