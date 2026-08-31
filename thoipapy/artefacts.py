"""Single source of truth for the paths of every file the pipeline reads and writes.

Before this, the same logical artefact was spelled out inline wherever it was needed -- 297
references to ``s["data_dir"]`` / ``s["sets_dir"]`` / ``s["base_dir"]``, across three different
idioms (``os.path.join``, ``Path`` division and ``str.format``), producing 155 distinct path
templates for a much smaller number of actual artefacts. The average artefact was spelled 1.7
times and the busiest ones up to nine times.

That was not merely untidy. The feature files carry the settings in their names
(``{acc}.surr{num_of_sur_residues}.gaps{max_n_gaps_in_TMD_subject_seq}.combined_features.csv``),
and while every writer built that suffix from the settings, five readers hardcoded
``surr20.gaps5``. Changing either setting therefore made the pipeline write to one filename and
read from another: at best a FileNotFoundError, at worst silently validating a model against stale
features left over from an earlier run at the old settings. Computing the suffix in exactly one
place removes that whole class of bug by construction.

The other benefit is that artefacts become greppable. ``grep -rn "combined_features_csv"`` finds
every producer and consumer; ``grep -rn "combined_features"`` did not, because some call sites
split the name across lines with ``.format``.

Usage::

    paths = ArtefactPaths.from_settings(s)
    df = pd.read_csv(paths.train_data_orig_csv())
"""

from dataclasses import dataclass
from pathlib import Path

from thoipapy.paths import BASE_DIR, SETS_DIR


@dataclass(frozen=True)
class ArtefactPaths:
    """Every pipeline path, derived from the data directory and the settings that appear in filenames.

    Frozen, because a path layer that can be mutated mid-run is how the pipeline ended up with two
    spellings of the same file in the first place.

    Attributes
    ----------
    data_dir : Path
        Root of the DVC-tracked data directory.
    setname : str
        Protein set being processed, e.g. "set08".
    num_of_sur_residues : int
        Residues either side of the TMD included in the homologue search. Appears in filenames.
    max_n_gaps_in_TMD_subject_seq : int
        Gap limit used when filtering homologues. Appears in filenames.
    """

    data_dir: Path
    setname: str
    num_of_sur_residues: int
    max_n_gaps_in_TMD_subject_seq: int

    @classmethod
    def from_settings(cls, s: dict) -> "ArtefactPaths":
        """Build from a settings dict, which must already carry setname."""
        return cls(
            data_dir=Path(s["data_dir"]),
            setname=s["setname"],
            num_of_sur_residues=int(s["num_of_sur_residues"]),
            max_n_gaps_in_TMD_subject_seq=int(s["max_n_gaps_in_TMD_subject_seq"]),
        )

    # ------------------------------------------------------------------ suffixes

    @property
    def surr_gaps_suffix(self) -> str:
        """The ``surr20.gaps5``-style fragment that appears in per-protein feature filenames."""
        return f"surr{self.num_of_sur_residues}.gaps{self.max_n_gaps_in_TMD_subject_seq}"

    @property
    def surr_suffix(self) -> str:
        """The ``surr20``-style fragment used by the homologue files, which have no gap filter."""
        return f"surr{self.num_of_sur_residues}"

    @property
    def set_number(self) -> int:
        """The numeric part of setname, e.g. 8 for "set08"."""
        return int(self.setname.removeprefix("set"))

    # ------------------------------------------------------------------ roots

    @property
    def results_dir(self) -> Path:
        return self.data_dir / "results" / self.setname

    @property
    def train_data_dir(self) -> Path:
        return self.results_dir / "train_data"

    @property
    def feat_imp_dir(self) -> Path:
        return self.results_dir / "feat_imp"

    @property
    def crossvalidation_dir(self) -> Path:
        return self.results_dir / "crossvalidation"

    @property
    def clusters_dir(self) -> Path:
        return self.results_dir / "clusters"

    @property
    def predictions_dir(self) -> Path:
        return self.results_dir / "predictions"

    @property
    def sets_dir(self) -> Path:
        return SETS_DIR

    @property
    def base_dir(self) -> Path:
        return BASE_DIR

    @property
    def protein_names_csv(self) -> Path:
        return BASE_DIR / "protein_names.csv"

    # ------------------------------------------------------- per-protein features

    def homologue_xml_tar(self, database: str, acc: str) -> Path:
        return self.data_dir / "homologues" / "xml" / database / f"{acc}.{self.surr_suffix}.BLAST.xml.tar.gz"

    def homologue_csv_tar(self, database: str, acc: str) -> Path:
        return self.data_dir / "homologues" / "ncbi" / database / f"{acc}.{self.surr_suffix}.BLAST.csv.tar.gz"

    def alignment_dir(self, database: str) -> Path:
        return self.data_dir / "homologues" / "alignments" / database

    def pssm_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "pssm" / database / f"{acc}.{self.surr_gaps_suffix}.pssm.csv"

    def entropy_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "entropy" / database / f"{acc}.{self.surr_gaps_suffix}.uniq.entropy.csv"

    def rate4site_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "rate4site" / database / f"{acc}_rate4site.csv"

    def freecontact_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "coevolution" / database / f"{acc}.{self.surr_gaps_suffix}.freecontact.csv"

    def freecontact_parsed_csv(self, database: str, acc: str) -> Path:
        return (
            self.data_dir
            / "features"
            / "coevolution"
            / database
            / f"{acc}.{self.surr_gaps_suffix}.freecontact_parsed.csv"
        )

    def combined_features_csv(self, database: str, acc: str, nohetero: bool = False) -> Path:
        """The per-protein feature table. Five readers used to hardcode surr20.gaps5 here.

        Parameters
        ----------
        nohetero : bool
            Select the variant with crystal hetero-contact residues already removed. Only crystal
            structures have one; see remove_crystal_hetero in the settings.
        """
        infix = ".nohetero" if nohetero else ""
        return (
            self.data_dir
            / "features"
            / "combined"
            / database
            / f"{acc}{infix}.{self.surr_gaps_suffix}.combined_features.csv"
        )

    # --- alignments -----------------------------------------------------------

    def uniq_tmd_seqs_for_pssm_freecontact_txt(self, database: str, acc: str) -> Path:
        return self.alignment_dir(database) / f"{acc}.{self.surr_gaps_suffix}.uniq.for_PSSM_FREECONTACT.txt"

    def uniq_tmd_seqs_surr5_for_lipo_txt(self, database: str, acc: str) -> Path:
        """Always surr5, regardless of num_of_sur_residues: lipophilicity uses a fixed 5-residue window."""
        return self.alignment_dir(database) / f"{acc}.surr5.gaps{self.max_n_gaps_in_TMD_subject_seq}.uniq.for_LIPO.txt"

    def uniq_tmd_seqs_no_gaps_for_lips_txt(self, database: str, acc: str) -> Path:
        """Always gaps0: LIPS requires an ungapped alignment."""
        return self.alignment_dir(database) / f"{acc}.{self.surr_suffix}.gaps0.uniq.for_LIPS.txt"

    def alignment_summary_csv(self, database: str, acc: str) -> Path:
        return self.alignment_dir(database) / f"{acc}.{self.surr_gaps_suffix}.alignment_summary.csv"

    def lips_output_csv(self, database: str, acc: str) -> Path:
        return self.alignment_dir(database) / f"{acc}.{self.surr_suffix}.LIPS_output.csv"

    # --- remaining per-protein features ---------------------------------------

    def pssm_surr5_csv(self, database: str, acc: str) -> Path:
        return (
            self.data_dir
            / "features"
            / "pssm"
            / database
            / f"{acc}.surr5.gaps{self.max_n_gaps_in_TMD_subject_seq}.pssm.csv"
        )

    def lipo_csv(self, database: str, acc: str, scalename: str) -> Path:
        return self.data_dir / "features" / "lipophilicity" / database / f"{acc}_{scalename}_lipo.csv"

    def relative_position_csv(self, database: str, acc: str, surres: str) -> Path:
        return self.data_dir / "features" / "relative_position" / database / f"{acc}.relative_position{surres}.csv"

    def lips_score_parsed_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "lips_score" / database / f"{acc}.{self.surr_suffix}.LIPS_score_parsed.csv"

    def lips_mem_output(self, database: str, acc: str, surres: str) -> Path:
        return self.data_dir / "features" / "lips_score" / database / f"{acc}.mem.lips.output{surres}"

    def motifs_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "motifs" / database / f"{acc}.motifs.csv"

    def merged_predictions_csv(self, database: str, acc: str) -> Path:
        return self.predictions_dir / "merged" / f"{database}.{acc}.merged.csv"

    def full_seq_fasta(self, database: str, acc: str) -> Path:
        return self.protein_dir(database) / f"{acc}.fasta"

    def blast_xml(self, database: str, acc: str) -> Path:
        return self.data_dir / "homologues" / "xml" / database / f"{acc}.{self.surr_suffix}.BLAST.xml"

    def structure_dir(self, database: str) -> Path:
        return self.data_dir / "features" / "structure" / database

    def other_predictors_dir(self, database: str) -> Path:
        return self.data_dir / "Predictions" / "other_predictors" / database

    def protein_dir(self, database: str) -> Path:
        return self.data_dir / "Proteins" / database

    # ------------------------------------------------------------- training table

    def train_data_orig_csv(self) -> Path:
        return self.train_data_dir / "01_train_data_orig.csv"

    def feat_imp_MDI_before_feature_seln_csv(self) -> Path:
        return self.train_data_dir / "01_feat_imp_MDI_before_feature_seln.csv"

    def tuned_ensemble_parameters_before_feature_seln_csv(self) -> Path:
        return self.train_data_dir / "01_tuned_ensemble_parameters_before_feature_seln.csv"

    def train_data_excl_duplicates_csv(self) -> Path:
        return self.train_data_dir / "02_train_data_excl_duplicates.csv"

    def train_data_after_first_feature_seln_csv(self) -> Path:
        return self.train_data_dir / "03_train_data_after_first_feature_seln.csv"

    def tuned_ensemble_parameters_csv(self) -> Path:
        return self.train_data_dir / "04_tuned_ensemble_parameters.csv"

    def hetero_contact_residues_csv(self) -> Path:
        return self.train_data_dir / "00_hetero_contact_residues.csv"

    def top_features_anova_csv(self) -> Path:
        return self.feat_imp_dir / "top_features_anova.csv"

    def top_features_rfe_csv(self) -> Path:
        return self.feat_imp_dir / "top_features_rfe.csv"

    # -------------------------------------------------------------------- model

    def ml_model_lpkl(self) -> Path:
        return self.results_dir / f"{self.setname}_ML_model.lpkl"

    # ------------------------------------------------------------------ clusters

    def protein_set_full_seq_fasta(self) -> Path:
        return self.clusters_dir / f"{self.setname}_full_seqs.fas"

    def protein_set_tmd_seq_fasta(self) -> Path:
        return self.clusters_dir / f"{self.setname}_tmd_seqs.fas"

    def cdhit_cluster_txt(self) -> Path:
        return self.clusters_dir / f"{self.setname}.fas.1.clstr.sorted.txt"

    # -------------------------------------------------------------- input tables

    def processed_set_csv(self, set_basename: str) -> Path:
        """Input_data/<setname>_processed.csv, written by process_set_protein_seqs."""
        return self.data_dir / "Input_data" / f"{set_basename}_processed.csv"

    # --- clustering, feature importance and t-test outputs ---------------------

    def sim_matrix_alignments_txt(self) -> Path:
        return self.clusters_dir / f"{self.setname}_sim_matrix_alignments.txt"

    def sim_matrix_xlsx(self) -> Path:
        return self.clusters_dir / f"{self.setname}_sim_matrix.xlsx"

    def feat_imp_mean_decrease_accuracy_xlsx(self) -> Path:
        return self.feat_imp_dir / "feat_imp_mean_decrease_accuracy.xlsx"

    def feat_imp_temp_bocurve_csv(self) -> Path:
        return self.feat_imp_dir / "feat_imp_temp_THOIPA.best_overlap_data.csv"

    def feat_imp_temp_bocurve_xlsx(self) -> Path:
        return self.feat_imp_dir / "feat_imp_temp_bocurve_data.xlsx"

    def variable_importance_png(self) -> Path:
        return self.feat_imp_dir / "FigS17_BZ13_feature_importance.png"

    def variable_importance_xlsx(self) -> Path:
        return self.feat_imp_dir / "FigS17_BZ13_feature_importance.xlsx"

    def results_remove_dup_feat_with_low_MDI_csv(self) -> Path:
        return self.feat_imp_dir / "results_remove_dup_feat_with_low_MDI.csv"

    def ttest_dir(self) -> Path:
        return self.results_dir / "ttest"

    def homohetero_contact_csv(self, database: str, acc: str) -> Path:
        return self.data_dir / "features" / "structure" / database / f"{acc}.homohetero.bind.closedist.csv"

    def experimental_interface_csv(self, database: str, acc: str, inter_pair_max: int) -> Path:
        return self.structure_dir(database) / f"{acc}.{inter_pair_max}pairmax.bind.closedist.csv"

    def etra_experimental_csv(self, acc: str) -> Path:
        return self.base_dir / "ETRA_data" / "Average_with_interface" / f"{acc}_mul_scan_average_data.csv"

    def logging_dir(self) -> Path:
        return self.data_dir / "Logging"
