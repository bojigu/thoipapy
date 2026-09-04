"""Which homologue source should THOIPA use? Depth, accuracy and cost for all three.

    THOIPA_UNIREF_DATA=/path/to/uniref90/rebuild \
    THOIPA_COLABFOLD_DATA=/path/to/colabfold/rebuild \
    python scripts/compare_homologue_sources.py

``compare_blast_database_depth.py`` asked whether a local UniRef90 database could stand in for the
``nr`` archive the paper used. This asks the question that actually has to be answered now, which
is a three-way one, because ``nr`` is no longer an option for anybody: NCBI does not archive old
releases, so the May 2020 snapshot behind the published numbers cannot be recovered, and a live
remote ``nr`` query has been measured queueing for hours. It is the reference, not a candidate.

That leaves two real candidates, and they differ in what they cost rather than in what they are:

    UniRef90 local    a single blastp pass. Tens of GB of server disk, rebuilt every release,
                      but it depends on nothing outside the machine.
    ColabFold         MMseqs2, three profile iterations against UniRef30 plus an environmental
                      database, run by somebody else's server. No disk, no maintenance, but a
                      dependency on a free academic service with a fair-use budget.

Both rebuilds must already exist; ``build_colabfold_dataset.py`` makes the ColabFold one. Every
dataset is read, none is written.

Accuracy is reported per CD-HIT cluster rather than per protein. Several proteins in set08 are
homologues of each other, so treating them as independent would count the same evidence more than
once and make the confidence interval too narrow.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compare_blast_database_depth import (  # noqa: E402
    NR_DATA,
    TEST_SET,
    TRAIN_SET,
    alignment_depth,
    bootstrap_ci,
    cluster_means,
    leave_one_out,
    selected_features,
    train_on_trainset_test_on_testset,
    training_table,
)

NR_LABEL = "nr 2020"
UNIREF_LABEL = "UniRef90 local"
COLABFOLD_LABEL = "ColabFold"


def _dir_from_env(variable: str, purpose: str) -> Path:
    value = os.environ.get(variable, "").strip()
    if not value:
        raise SystemExit(f"{variable} is not set. Point it at {purpose}.")
    directory = Path(value)
    if not directory.is_dir():
        raise SystemExit(f"{variable} points at {directory}, which is not a directory.")
    return directory


def report_depth(sources: dict[str, Path]) -> None:
    print("=" * 92)
    print("1. ALIGNMENT DEPTH: unique TMD homologues surviving the identity and gap filters")
    print("=" * 92)

    rows = []
    depths: dict[str, dict[str, pd.Series]] = {}
    for label, data_dir in sources.items():
        depths[label] = {}
        row: dict[str, str | int] = {"source": label}
        for setname in (TRAIN_SET, TEST_SET):
            series = alignment_depth(data_dir, setname)
            depths[label][setname] = series
            row[f"{setname} median"] = int(series.median())
            row[f"{setname} total"] = int(series.sum())
        rows.append(row)
    print(pd.DataFrame(rows).to_string(index=False))

    baseline = depths[UNIREF_LABEL][TRAIN_SET]
    candidate = depths[COLABFOLD_LABEL][TRAIN_SET]
    common = baseline.index.intersection(candidate.index)
    ratio = candidate[common] / baseline[common].replace(0, np.nan)
    print(
        f"\n  {COLABFOLD_LABEL} against {UNIREF_LABEL} on {TRAIN_SET}: median {ratio.median():.2f}x, "
        f"deeper on {int((candidate[common] > baseline[common]).sum())} of {len(common)} proteins"
    )


def report_accuracy(sources: dict[str, Path]) -> None:
    print("\n" + "=" * 92)
    print("2. PREDICTIVE ACCURACY: identical code, identical folds, only the homologue source differs")
    print("=" * 92)

    loo: dict[str, tuple[pd.Series, pd.Series]] = {}
    blind: dict[str, pd.Series] = {}
    clusters: list = []

    for label, data_dir in sources.items():
        roc, pr, clusters = leave_one_out(data_dir, TRAIN_SET)
        loo[label] = (roc, pr)
        line = f"  {label:<16} {TRAIN_SET} leave-one-out ROC {roc.mean():.3f}  PR {pr.mean():.3f}"
        try:
            blind_roc, blind_pr = train_on_trainset_test_on_testset(data_dir)
            blind[label] = blind_roc
            line += f"   {TEST_SET} blind ROC {blind_roc.mean():.3f}  PR {blind_pr.mean():.3f}"
        except (FileNotFoundError, ValueError) as e:
            line += f"   {TEST_SET} blind unavailable ({type(e).__name__})"
        print(line)

    print(f"\n  paired differences over CD-HIT clusters, {TRAIN_SET} leave-one-out ROC AUC:")
    comparisons = (
        (COLABFOLD_LABEL, UNIREF_LABEL),
        (UNIREF_LABEL, NR_LABEL),
        (COLABFOLD_LABEL, NR_LABEL),
    )
    for candidate_label, baseline_label in comparisons:
        candidate, baseline = loo[candidate_label][0], loo[baseline_label][0]
        common = candidate.index.intersection(baseline.index)
        by_cluster = cluster_means(candidate[common] - baseline[common], clusters)
        low, high = bootstrap_ci(by_cluster)
        p = wilcoxon(by_cluster).pvalue
        print(
            f"    {candidate_label:<16} minus {baseline_label:<16} "
            f"{by_cluster.mean():+.3f} [95% CI {low:+.3f}, {high:+.3f}]  Wilcoxon p = {p:.3f}"
        )
    print(f"    ({len(by_cluster)} clusters, not {len(common)} proteins: several are homologues)")

    if COLABFOLD_LABEL in blind and UNIREF_LABEL in blind:
        candidate, baseline = blind[COLABFOLD_LABEL], blind[UNIREF_LABEL]
        shared = candidate.index.intersection(baseline.index)
        print(
            f"    {TEST_SET} blind: {COLABFOLD_LABEL} minus {UNIREF_LABEL} "
            f"{(candidate[shared] - baseline[shared]).mean():+.3f} "
            f"({len(shared)} proteins, no confidence interval, descriptive only)"
        )


def report_feature_agreement(sources: dict[str, Path]) -> None:
    print("\n" + "=" * 92)
    print("3. FEATURE AGREEMENT with the published nr features, per feature Spearman")
    print("=" * 92)

    published = training_table(NR_DATA, TRAIN_SET)
    features = selected_features(NR_DATA, TRAIN_SET)

    for label in (UNIREF_LABEL, COLABFOLD_LABEL):
        table = training_table(sources[label], TRAIN_SET)
        rows = published.index.intersection(table.index)
        correlations = pd.Series(
            {
                feature: spearmanr(published.loc[rows, feature], table.loc[rows, feature]).statistic
                for feature in features
                if feature in table.columns
                and published.loc[rows, feature].nunique() > 1
                and table.loc[rows, feature].nunique() > 1
            }
        )
        # DI and MI are the FreeContact coevolution features, the ones computed from the
        # alignment's covariation rather than from its column composition.
        coevolution = [f for f in correlations.index if f.startswith(("DI", "MI"))]
        print(
            f"  {label:<16} median rho {correlations.median():+.3f}   "
            f"coevolution features median rho {correlations[coevolution].median():+.3f}   "
            f"worst {correlations.idxmin()} {correlations.min():+.3f}"
        )

    print(
        "\n  A source that scores well here reproduces the published behaviour; it says nothing\n"
        "  about which source is more accurate. Read it alongside section 2, not instead of it."
    )


def main() -> None:
    sources = {
        NR_LABEL: NR_DATA,
        UNIREF_LABEL: _dir_from_env("THOIPA_UNIREF_DATA", "a data directory rebuilt from a local UniRef90 search"),
        COLABFOLD_LABEL: _dir_from_env("THOIPA_COLABFOLD_DATA", "the output of build_colabfold_dataset.py"),
    }

    report_depth(sources)
    report_accuracy(sources)
    report_feature_agreement(sources)

    print(
        "\nEvery dataset was read, none written. Both rebuilds are intentionally not DVC-tracked:\n"
        "the published nr artefacts remain the reference."
    )


if __name__ == "__main__":
    main()
