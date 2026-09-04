"""Rebuild the feature set from ColabFold alignments, so it can be compared with the tracked one.

    THOIPA_COLABFOLD_CONTACT_EMAIL=you@example.com \
    THOIPA_A3M_DIR=/path/to/downloaded/a3m/data \
    THOIPA_COLABFOLD_DATA=/path/for/the/rebuild \
    python scripts/build_colabfold_dataset.py

Changing the homologue source changes far more than the alignments. Conservation, entropy,
rate4site, the PSSM, the lipophilicity derived from it, the LIPS scores and every FreeContact
coevolution feature are all computed from the alignment, so a fair comparison has to rebuild all
of them and retrain, not just count sequences.

This builds a parallel data directory to do that. Everything that does not depend on the
alignment is copied from the tracked directory -- structure-derived residue depth, the
experimental interfaces, the CD-HIT clusters that define the cross-validation folds, the protein
sets -- and everything that does is deleted and recomputed. What comes out has the same shape as
the published table, the same residues and the same labels, so the two can be scored by identical
code with only the feature values differing. ``compare_homologue_sources.py`` does that scoring.

The repository's own data/ directory is read, never written. The rebuild lives outside it and is
deliberately not DVC-tracked: the published artefacts stay the reference.

Alignments have to be downloaded first, which ``compare_colabfold_alignment_depth.py`` does as a
side effect of measuring depth. Point THOIPA_A3M_DIR at the directory it wrote.
"""

import os
import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT))

from thoipapy.paths import DATA_DIR, RUN_SETTINGS_CSV  # noqa: E402

# Which sets to rebuild. set08 trains, set07 is the blind test, and the comparison needs both.
SET_NUMBERS = (8, 7)

# Must match the mode the alignments were downloaded with, because it decides whether the
# environmental hits in each archive are read.
COLABFOLD_MODE = "env-nofilter"

# Recomputed from the new alignments. Everything else under data/ is alignment-independent.
ALIGNMENT_DERIVED_FEATURE_DIRS = (
    "pssm",
    "entropy",
    "rate4site",
    "coevolution",
    "lipophilicity",
    "lips_score",
    "motifs",
    "relative_position",
    "combined",
)


def _dir_from_env(variable: str, purpose: str) -> Path:
    """Read a directory from the environment, refusing to guess.

    No defaults. These directories are large, they are not the reference dataset, and where they
    live is a property of the machine rather than of this repository.
    """
    value = os.environ.get(variable, "").strip()
    if not value:
        raise SystemExit(f"{variable} is not set. Point it at {purpose}.")
    return Path(value)


def prepare_data_dir(output_dir: Path, a3m_dir: Path) -> None:
    """Copy the tracked data directory, then clear everything derived from the old alignments."""
    if output_dir.exists():
        print(f"removing previous {output_dir}")
        shutil.rmtree(output_dir)

    print(f"copying {DATA_DIR} -> {output_dir}")
    shutil.copytree(DATA_DIR, output_dir, symlinks=True)

    # Deleted rather than left to be overwritten, so that a stage which quietly fails shows up as
    # a missing file instead of silently contributing the old source's numbers to the comparison.
    for subdir in ("ncbi", "xml", "alignments", "colabfold", "a3m"):
        target = output_dir / "homologues" / subdir
        if target.exists():
            shutil.rmtree(target)
            print(f"  cleared homologues/{subdir}")

    for name in ALIGNMENT_DERIVED_FEATURE_DIRS:
        target = output_dir / "features" / name
        if target.exists():
            shutil.rmtree(target)
            print(f"  cleared features/{name}")

    # The training tables and the feature selection built from them. clusters/ is deliberately
    # kept: CD-HIT clusters come from the query sequences, so they are the same whatever the
    # homologue source, and reusing them keeps the cross-validation folds identical.
    for set_number in SET_NUMBERS:
        train_data = output_dir / "results" / f"set{set_number:02d}" / "train_data"
        if train_data.exists():
            shutil.rmtree(train_data)
            print(f"  cleared results/set{set_number:02d}/train_data")

    source = a3m_dir / "homologues" / "a3m"
    if not source.is_dir():
        raise SystemExit(f"no alignments at {source}. Run compare_colabfold_alignment_depth.py first.")
    shutil.copytree(source, output_dir / "homologues" / "a3m")
    print(f"  {len(list((output_dir / 'homologues' / 'a3m').rglob('*.a3m.tar.gz')))} alignments in place")


def settings(output_dir: Path) -> dict:
    """The published settings, with the homologue source switched and the rebuild stages on."""
    import thoipapy.common

    s = thoipapy.common.create_settingdict(RUN_SETTINGS_CSV)
    s["data_dir"] = output_dir
    s["homologue_source"] = "colabfold"
    s["colabfold_mode"] = COLABFOLD_MODE

    # The alignments are already downloaded, so the search stage stays off.
    s["run_retrieve_homologues_from_colabfold"] = False

    for stage in (
        "run_parse_colabfold_a3m_into_csv",
        "parse_csv_homologues_to_alignment",
        "pssm_calculation",
        "entropy_calculation",
        "rate4site_calculation",
        "coevolution_calculation",
        "clac_relative_position",
        "calc_lipo_from_pssm",
        "lips_score_calculation",
        "motifs_from_seq",
        "combine_feature_into_train_data",
        "run_feature_selection",
    ):
        s[stage] = True

    # Off: the NCBI path, and everything past feature selection. The comparison fits its own
    # forests from the tuned parameters, so training and validation are not needed here.
    for stage in (
        "run_retrieve_NCBI_homologues_with_blastp",
        "run_parse_homologues_xml_into_csv",
        "download_10_homologues_from_ncbi",
        "train_machine_learning_model",
        "calc_feature_importances",
        "conduct_ttest",
        "run_validation",
        "run_testset_trainset_validation",
        "create_identity_matrix_from_set_seqs",
    ):
        s[stage] = False

    # Reproduces the published 850-residue training table. Asserted rather than assigned, so that
    # a change to the shipped settings shows up here instead of quietly producing 887 residues.
    if s["remove_crystal_hetero"] is not True:
        raise SystemExit(f"remove_crystal_hetero is {s['remove_crystal_hetero']!r}, expected True")
    return s


def main() -> None:
    output_dir = _dir_from_env("THOIPA_COLABFOLD_DATA", "a directory outside the repository for the rebuild")
    a3m_dir = _dir_from_env("THOIPA_A3M_DIR", "the directory holding the downloaded a3m archives")

    prepare_data_dir(output_dir, a3m_dir)

    import thoipapy.run

    s = settings(output_dir)
    for set_number in SET_NUMBERS:
        started = time.time()
        print(f"\n{'=' * 80}\nrebuilding set{set_number:02d}\n{'=' * 80}", flush=True)
        thoipapy.run.run_one_set(s, set_number)
        print(f"set{set_number:02d} finished in {(time.time() - started) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
