"""Does alignment depth on its own change THOIPA's accuracy?

    THOIPA_COLABFOLD_DATA=/path/to/env-nofilter/rebuild \
    THOIPA_COLABFOLD_ENV_DATA=/path/to/env/rebuild \
    python scripts/compare_alignment_depth_control.py

``compare_homologue_sources.py`` compares three homologue sources, and every one of those
comparisons is confounded. The 2020 ``nr`` archive differs from a 2026 search by six years of
sequence growth as well as by depth; a local UniRef90 blastp differs from an MMseqs2 profile search
by algorithm as well as by depth. Neither can separate "more sequences" from "different sequences".

This one can. Both arms come from the same server, the same three-iteration MMseqs2 search, the
same query and the same downstream filters. The only difference is the ``colabfold_mode``: with the
diversity filter on, MMseqs2 applies ``--qsc 0.8 --max-seq-id 0.95`` and thins near-identical hits;
with it off, they are kept. That is a pure depth manipulation of one search, so the comparison is
free of the time and algorithm confounds.

Two further variables are pinned so that only the feature values move:

* the feature set, taken from the published ``nr`` selection rather than each arm's own
* the hyperparameters, likewise

In ``compare_homologue_sources.py`` both are refit per arm, which is the right way to compare
sources as they would actually be deployed, but it means a difference there could come from a
different feature set rather than from the alignment. Pinning them answers the narrower question.

Reading the result: pinning to ``nr``'s selection biases the *nr* comparisons in nr's favour, since
those features were chosen on nr's values. It does not bias the depth comparison, because neither
arm had a hand in choosing them.
"""

import os
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402

from compare_blast_database_depth import (  # noqa: E402
    BIND_COLUMN,
    FOREST_SEEDS,
    NR_DATA,
    TRAIN_SET,
    alignment_depth,
    bootstrap_ci,
    cluster_means,
    forest,
    homologue_clusters,
    selected_features,
    training_table,
)

DEEP_LABEL = "ColabFold env-nofilter"
SHALLOW_LABEL = "ColabFold env"


def _dir_from_env(variable: str, purpose: str) -> Path:
    value = os.environ.get(variable, "").strip()
    if not value:
        raise SystemExit(f"{variable} is not set. Point it at {purpose}.")
    directory = Path(value)
    if not directory.is_dir():
        raise SystemExit(f"{variable} points at {directory}, which is not a directory.")
    return directory


def leave_one_out_with_pinned_model(data_dir: Path, features: list[str], tuned_csv: Path):
    """Leave-one-out over CD-HIT clusters, with the feature set and hyperparameters supplied.

    Identical to ``compare_blast_database_depth.leave_one_out`` except that it does not read the
    feature selection or the tuning from ``data_dir``, so every arm is scored by the same model
    and only the feature values differ. Folds come from the published clustering for the same
    reason.
    """
    df = training_table(data_dir, TRAIN_SET)
    missing = [f for f in features if f not in df.columns]
    if missing:
        raise SystemExit(f"{data_dir} is missing pinned features: {missing}")

    y, groups = df[BIND_COLUMN], df["acc_db"]
    homologues, clusters = homologue_clusters(NR_DATA, TRAIN_SET)

    roc, pr = {}, {}
    for acc_db in groups.unique():
        train = (~groups.isin(homologues.get(acc_db, {acc_db}))).to_numpy()
        test = (groups == acc_db).to_numpy()
        if test.sum() == 0 or y[test].nunique() < 2:
            continue
        X = df[features]
        r, p = [], []
        for seed in FOREST_SEEDS:
            fitted = forest(tuned_csv, seed).fit(X[train], y[train])
            prediction = fitted.predict_proba(X[test])[:, 1]
            r.append(roc_auc_score(y[test], prediction))
            p.append(average_precision_score(y[test], prediction))
        roc[acc_db], pr[acc_db] = float(sum(r) / len(r)), float(sum(p) / len(p))
    return pd.Series(roc), pd.Series(pr), clusters


def main() -> None:
    deep = _dir_from_env("THOIPA_COLABFOLD_DATA", "the env-nofilter rebuild")
    shallow = _dir_from_env("THOIPA_COLABFOLD_ENV_DATA", "the env rebuild")

    features = selected_features(NR_DATA, TRAIN_SET)
    tuned_csv = NR_DATA / f"results/{TRAIN_SET}/train_data/04_tuned_ensemble_parameters.csv"

    print("=" * 92)
    print("DEPTH, WITH EVERYTHING ELSE HELD CONSTANT")
    print("=" * 92)
    print("  same server, same search, same query, same filters, same folds; only colabfold_mode differs")
    print(f"  model pinned to the published nr selection: {len(features)} features, {tuned_csv.name}\n")

    depths = {}
    for label, data_dir in ((SHALLOW_LABEL, shallow), (DEEP_LABEL, deep)):
        d = alignment_depth(data_dir, TRAIN_SET)
        depths[label] = d
        print(f"  {label:<24} median depth {int(d.median()):>6}   total {int(d.sum()):>8,}")

    common = depths[SHALLOW_LABEL].index.intersection(depths[DEEP_LABEL].index)
    median_ratio = (depths[DEEP_LABEL][common] / depths[SHALLOW_LABEL][common]).median()
    total_ratio = depths[DEEP_LABEL][common].sum() / depths[SHALLOW_LABEL][common].sum()
    print(f"  depth manipulation: median {median_ratio:.2f}x, total {total_ratio:.2f}x\n")

    results = {}
    for label, data_dir in ((SHALLOW_LABEL, shallow), (DEEP_LABEL, deep)):
        roc, pr, clusters = leave_one_out_with_pinned_model(data_dir, features, tuned_csv)
        results[label] = (roc, pr, clusters)
        print(
            f"  {label:<24} {TRAIN_SET} LOO ROC {cluster_means(roc, clusters).mean():.4f}  "
            f"PR {cluster_means(pr, clusters).mean():.4f}   (per cluster)"
        )

    deep_roc, deep_pr, clusters = results[DEEP_LABEL]
    shallow_roc, shallow_pr, _ = results[SHALLOW_LABEL]
    shared = deep_roc.index.intersection(shallow_roc.index)

    print(f"\n  paired over {len(cluster_means(deep_roc[shared] - shallow_roc[shared], clusters))} clusters:")
    for metric, deep_s, shallow_s in (("ROC AUC", deep_roc, shallow_roc), ("PR AUC", deep_pr, shallow_pr)):
        by_cluster = cluster_means(deep_s[shared] - shallow_s[shared], clusters)
        low, high = bootstrap_ci(by_cluster)
        p = wilcoxon(by_cluster).pvalue
        print(
            f"    deeper minus shallower, {metric:8s} {by_cluster.mean():+.4f} "
            f"[95% CI {low:+.4f}, {high:+.4f}]  Wilcoxon p = {p:.3f}"
        )

    print(
        "\n  A difference indistinguishable from zero here is the cleanest evidence available that\n"
        "  depth alone does not drive accuracy, because nothing else varies between the two arms.\n"
        "  It does not rule out that a deeper alignment built by a different search, or filtered\n"
        "  differently, would help: depth and alignment composition are still varied together."
    )


if __name__ == "__main__":
    main()
