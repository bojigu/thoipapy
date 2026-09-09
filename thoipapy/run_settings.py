"""Typed representation of the pipeline stage switches.

Everything below the orchestrator now takes explicit parameters, which is clearer than passing a
dictionary around. ``run_one_set`` is the one place where that does not work: it reads 31 on/off
flags, and 31 boolean arguments would be worse than what it replaced. Those flags are also a
different kind of thing from the rest of the settings -- they are not parameters of a calculation,
they are a description of which stages to run.

So they get a type of their own. The gain over the raw dictionary is that a typo is caught when
the settings are loaded rather than an hour into a run, a missing flag has a documented default
instead of raising KeyError somewhere in the middle of the pipeline, and the set of stages is
discoverable by reading one class instead of grepping for ``if s[``.

The boolean coercion matters more than it looks. The ``type`` column in the settings CSV
only exists on one of the three sheets and is blank on eight rows of that one, so values used to
arrive as the string ``"TRUE"``. That worked by accident -- a non-empty string is truthy -- but so
is ``"FALSE"``, which meant a stage switched off in the spreadsheet could still run.
"""

from dataclasses import dataclass, fields
from typing import Any

TRUE_STRINGS = {"true", "t", "yes", "y", "1", "wahr"}
FALSE_STRINGS = {"false", "f", "no", "n", "0", "falsch"}


def _as_bool(value: Any, name: str) -> bool:
    """Coerce a settings value to a real bool, refusing anything ambiguous.

    Raises
    ------
    ValueError
        If the value cannot be read as a boolean. Better to fail at load time than to let
        "FALSE" run a stage because non-empty strings are truthy.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, int | float) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in TRUE_STRINGS:
            return True
        if lowered in FALSE_STRINGS:
            return False
    raise ValueError(
        f"settings flag {name!r} has the value {value!r}, which is not recognisably true or false. "
        f"Use TRUE or FALSE in the settings CSV."
    )


@dataclass(frozen=True)
class RunSettings:
    """Which pipeline stages thoipapy.run.run_one_set should execute.

    Every field defaults to False, so a flag missing from the spreadsheet means "do not run this
    stage" rather than KeyError partway through a run.
    """

    # homologues
    run_retrieve_NCBI_homologues_with_blastp: bool = False
    download_10_homologues_from_ncbi: bool = False
    run_parse_homologues_xml_into_csv: bool = False
    # The ColabFold MSA server is an alternative to the two stages above, not an addition to
    # them: it replaces the blastp search and the xml parsing with a search and an a3m parse.
    # parse_csv_homologues_to_alignment then runs unchanged on whichever csv files were written.
    run_retrieve_homologues_from_colabfold: bool = False
    run_parse_colabfold_a3m_into_csv: bool = False
    parse_csv_homologues_to_alignment: bool = False

    # per-residue features
    pssm_calculation: bool = False
    entropy_calculation: bool = False
    rate4site_calculation: bool = False
    coevolution_calculation: bool = False
    clac_relative_position: bool = False  # spelling preserved: it is the spreadsheet column name
    calc_lipo_from_pssm: bool = False
    lips_score_calculation: bool = False
    motifs_from_seq: bool = False

    # assembling the training table
    combine_feature_into_train_data: bool = False
    remove_crystal_hetero: bool = False
    generate_randomised_interfaces: bool = False
    add_PREDDIMER_TMDOCK_to_combined_features: bool = False
    calc_PREDDIMER_TMDOCK_closedist: bool = False

    # model
    run_feature_selection: bool = False
    calc_feature_importances: bool = False
    conduct_ttest: bool = False
    train_machine_learning_model: bool = False

    # validation
    run_validation: bool = False
    run_testset_trainset_validation: bool = False

    # structures and experimental data
    Atom_Close_Dist: bool = False
    calc_NMR_closedist: bool = False
    Get_Tmd_Homodimers: bool = False

    # clustering
    create_identity_matrix_from_set_seqs: bool = False

    # publication figures (thoipapy/paper_figures, not maintained)
    compare_selected_predictors: bool = False
    create_merged_heatmap_for_trainset_and_testset: bool = False
    retrospective_coevolution: bool = False
    plot_coev_vs_res_dist: bool = False

    @classmethod
    def from_settings(cls, s: dict) -> "RunSettings":
        """Build from the settings dict, coercing each flag and rejecting unrecognised values.

        Parameters
        ----------
        s : dict
            Settings dictionary from thoipapy.common.create_settingdict.

        Returns
        -------
        RunSettings
            One field per pipeline stage. Flags absent from the spreadsheet default to False.
        """
        kwargs = {}
        for field in fields(cls):
            if field.name in s:
                kwargs[field.name] = _as_bool(s[field.name], field.name)
        return cls(**kwargs)

    def enabled_stages(self) -> list[str]:
        """Names of the stages that are switched on, for logging what a run will actually do."""
        return [f.name for f in fields(self) if getattr(self, f.name)]
