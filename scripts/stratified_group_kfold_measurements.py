"""Regenerate every number quoted in docs/stratified_group_kfold_plan.md.

The plan proposes replacing the set08/set07 train/blind-test division with repeated
StratifiedGroupKFold over set05, grouped by homology cluster and stratified by experimental
database. This script produces the evidence it cites, so the numbers can be checked rather than
taken on trust.

Printed sections are numbered to match the plan:

    2.1  data tables, cluster structure, clustering-threshold sweep
    2.3  blind test bootstrap intervals
    2.4  splitting-scheme comparison over 10 seeds
    2.5  the proposed design's own interval, and the paired instrument
    3.4  fold composition at k=5 and k=10, distinct groups per stratum
    3.5  shuffle behaviour and partition diversity
    6    per-fit timing
    4    the column-filter check

Reads only from data/results/ and thoipapy/setting/. Writes nothing. Takes about four minutes.

Requires the DVC-tracked data directory, so run `dvc pull` first. Run with the project
environment active:

    source activate thoipapy
    python scripts/stratified_group_kfold_measurements.py
"""

import time
from ast import literal_eval
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold, StratifiedKFold

REPO_ROOT = Path(__file__).parents[1]
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = DATA_DIR / "results"
MODEL_FEATURES_CSV = REPO_ROOT / "thoipapy" / "setting" / "model_features.csv"

# The classifier used for every measurement. Deliberately not the shipped model's parameters:
# a single fixed configuration keeps the comparisons between splitting schemes clean. Section 2.3
# of the plan notes what changes when the shipped parameters are used instead.
CLASSIFIER_PARAMS = {
    "n_estimators": 200,
    "min_samples_leaf": 4,
    "max_features": "sqrt",
    "bootstrap": False,
    "random_state": 0,
    "n_jobs": -1,
}

N_BOOTSTRAP_RESAMPLES = 4000
N_SEEDS = 10
N_FOLDS = 5


def classifier(**overrides):
    """Return the fixed ExtraTreesClassifier, optionally overriding a parameter."""
    return ExtraTreesClassifier(**{**CLASSIFIER_PARAMS, **overrides})


def sim_matrix_xlsx(setname):
    return RESULTS_DIR / setname / "clusters" / f"{setname}_sim_matrix.xlsx"


def train_data_orig_csv(setname):
    return RESULTS_DIR / setname / "train_data" / "01_train_data_orig.csv"


def selected_features():
    """The 26 features the shipped model uses, taken from the set08 selection artefact."""
    path = RESULTS_DIR / "set08" / "train_data" / "03_train_data_after_first_feature_seln.csv"
    return [c for c in pd.read_csv(path, index_col=0, nrows=1).columns if c != "interface"]


def homology_groups(setname):
    """Map each protein to an integer homology-cluster id, from the reduced_clusters sheet."""
    df = pd.read_excel(sim_matrix_xlsx(setname), sheet_name="reduced_clusters", index_col=0)
    clusters = df["reduced_clusters"].apply(literal_eval)
    # Cluster members are prefixed with the cd-hit numbering ("4-2j58A1-crystal"); strip it.
    members = clusters.apply(lambda cluster: ["-".join(name.split("-")[1:]) for name in cluster])
    protein_to_cluster = {}
    for cluster in members:
        for protein in cluster:
            protein_to_cluster[protein] = frozenset(cluster)
    cluster_ids = {cluster: i for i, cluster in enumerate(dict.fromkeys(protein_to_cluster.values()))}
    return {protein: cluster_ids[cluster] for protein, cluster in protein_to_cluster.items()}


def homology_groups_at_cutoff(setname, cutoff):
    """Re-cluster from the raw similarity matrix at a different identity cutoff.

    Reproduces the transitive closure in clustering/pairwise_aln_similarity_matrix.py, so the
    sensitivity of the grouping to its 15% threshold can be measured.
    """
    sim = pd.read_excel(sim_matrix_xlsx(setname), sheet_name="sim_matrix", index_col=0)
    names = ["-".join(name.split("-")[1:]) for name in sim.index]
    sim.index = names
    sim.columns = names
    parent = dict(zip(names, names))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for i, first in enumerate(names):
        for second in names[i + 1 :]:
            if sim.at[first, second] > cutoff:
                parent[find(first)] = find(second)

    components = {}
    for name in names:
        components.setdefault(find(name), []).append(name)
    return {name: i for i, members in enumerate(components.values()) for name in members}


def load_set(setname):
    """The stacked per-residue training table, with the database parsed out of acc_db."""
    df = pd.read_csv(train_data_orig_csv(setname), index_col=0)
    df["db"] = df["acc_db"].str.split("-").str[-1]
    return df


def out_of_fold_predictions(X, y, splitter, split_args):
    """Fit one model per fold and return the concatenated out-of-fold probabilities."""
    predictions = np.full(len(y), np.nan)
    for train_idx, test_idx in splitter.split(*split_args):
        fitted = classifier().fit(X.iloc[train_idx], y.iloc[train_idx])
        predictions[test_idx] = fitted.predict_proba(X.iloc[test_idx])[:, 1]
    return predictions


def _resample_indices(unit_labels, n_resamples, seed):
    """Yield index arrays from resampling whole units (clusters or proteins) with replacement."""
    rng = np.random.default_rng(seed)
    units = np.array(sorted(set(unit_labels)))
    positions = {unit: np.where(unit_labels == unit)[0] for unit in units}
    for _ in range(n_resamples):
        drawn = rng.choice(units, len(units), replace=True)
        yield np.concatenate([positions[unit] for unit in drawn])


def bootstrap_auc_ci(y, predictions, unit_labels, n_resamples=N_BOOTSTRAP_RESAMPLES, seed=0):
    """Percentile interval for the AUC, resampling whole units rather than residues.

    Residues within one TMD are not independent draws, so resampling rows would return an
    interval narrower than the data support.
    """
    scores = [
        roc_auc_score(y[idx], predictions[idx])
        for idx in _resample_indices(unit_labels, n_resamples, seed)
        if len(set(y[idx])) == 2
    ]
    return np.percentile(scores, [2.5, 97.5])


def bootstrap_auc_difference_ci(
    y, predictions_a, predictions_b, unit_labels, n_resamples=N_BOOTSTRAP_RESAMPLES, seed=0
):
    """Percentile interval for the AUC difference between two arms on identical folds.

    This is the estimator the plan nominates as its primary inferential claim: pairing removes
    the between-resample variation that dominates either arm's marginal interval.
    """
    differences = [
        roc_auc_score(y[idx], predictions_a[idx]) - roc_auc_score(y[idx], predictions_b[idx])
        for idx in _resample_indices(unit_labels, n_resamples, seed)
        if len(set(y[idx])) == 2
    ]
    return np.percentile(differences, [2.5, 97.5])


def partition_signature(groups, strata, seed, n_folds=N_FOLDS, shuffle=True):
    """A canonical, order-independent identity for one group-to-fold assignment."""
    splitter = StratifiedGroupKFold(n_folds, shuffle=shuffle, random_state=seed if shuffle else None)
    return tuple(
        tuple(sorted(set(groups.iloc[test_idx])))
        for _train_idx, test_idx in splitter.split(groups, strata, groups=groups)
    )


def report_data_and_clusters():
    print("=" * 78)
    print("2.1  data tables and cluster structure")
    print("=" * 78)
    for setname in ["set05", "set07", "set08"]:
        df = load_set(setname)
        summary = df.groupby("db").agg(
            n_res=("interface", "size"), n_pos=("interface", "sum"), n_prot=("acc_db", "nunique")
        )
        summary["frac_pos"] = (summary["n_pos"] / summary["n_res"]).round(3)
        print(
            f"\n{setname}: {len(df)} residues, {df['acc_db'].nunique()} proteins, "
            f"{int(df['interface'].sum())} positives ({df['interface'].mean():.3f})"
        )
        print(summary.to_string())
        cluster_ids = pd.Series(homology_groups(setname))
        sizes = cluster_ids.value_counts().value_counts().sort_index()
        print(f"  clusters: {cluster_ids.nunique()} groups, size histogram {dict(sizes)}")

    print("\nclustering-threshold sweep on set05 (identity cutoff -> n groups, largest cluster):")
    for cutoff in [20, 17, 15, 14, 11, 10]:
        counts = pd.Series(homology_groups_at_cutoff("set05", cutoff)).value_counts()
        print(f"  {cutoff:>2}%  ->  {len(counts):>2} groups, largest {counts.max()}")


def report_blind_test():
    print("\n" + "=" * 78)
    print("2.3  the blind test set's interval")
    print("=" * 78)
    train = pd.read_csv(
        RESULTS_DIR / "set08" / "train_data" / "03_train_data_after_first_feature_seln.csv", index_col=0
    )
    test = load_set("set07")
    X_test = test[selected_features()].dropna()
    y_test = test.loc[X_test.index, "interface"].astype(int)
    proteins = test.loc[X_test.index, "acc_db"]

    fitted = classifier().fit(train.drop(columns=["interface"]), train["interface"].astype(int))
    predictions = fitted.predict_proba(X_test)[:, 1]

    print(
        f"n_res={len(y_test)} n_pos={int(y_test.sum())} n_prot={proteins.nunique()}  "
        f"AUC={roc_auc_score(y_test, predictions):.3f}"
    )
    low, high = bootstrap_auc_ci(y_test.values, predictions, np.arange(len(y_test)))
    print(f"  resampling residues: {low:.3f}-{high:.3f}  width {high - low:.3f}")
    low, high = bootstrap_auc_ci(y_test.values, predictions, proteins.values)
    print(f"  resampling proteins: {low:.3f}-{high:.3f}  width {high - low:.3f}")

    per_protein = [
        roc_auc_score(y_test[proteins == acc], predictions[(proteins == acc).values])
        for acc in proteins.unique()
        if y_test[proteins == acc].nunique() == 2
    ]
    print(
        f"  per-protein AUC {np.round(sorted(per_protein), 2).tolist()}  "
        f"mean {np.mean(per_protein):.3f} sd {np.std(per_protein, ddof=1):.3f}"
    )


def report_splitting_schemes(df, X, y, groups, strata_db_label, strata_db):
    print("\n" + "=" * 78)
    print(f"2.4  splitting-scheme comparison, k={N_FOLDS}, shuffle=True, mean over seeds 0-{N_SEEDS - 1}")
    print("=" * 78)
    schemes = {
        "StratifiedKFold, no grouping": (
            lambda seed: StratifiedKFold(N_FOLDS, shuffle=True, random_state=seed),
            lambda: (X, y),
        ),
        "GroupKFold by protein": (
            lambda seed: GroupKFold(N_FOLDS, shuffle=True, random_state=seed),
            lambda: (X, y, df["acc_db"]),
        ),
        "SGKF by protein, strata=db x label": (
            lambda seed: StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=seed),
            lambda: (X, strata_db_label, df["acc_db"]),
        ),
        "SGKF by cluster, strata=db x label": (
            lambda seed: StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=seed),
            lambda: (X, strata_db_label, groups),
        ),
        "SGKF by cluster, strata=db only": (
            lambda seed: StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=seed),
            lambda: (X, strata_db, groups),
        ),
    }
    for name, (make_splitter, make_args) in schemes.items():
        aucs, average_precisions = [], []
        for seed in range(N_SEEDS):
            predictions = out_of_fold_predictions(X, y, make_splitter(seed), make_args())
            aucs.append(roc_auc_score(y, predictions))
            average_precisions.append(average_precision_score(y, predictions))
        print(
            f"  {name:38s} AUC {np.mean(aucs):.3f} ({np.std(aucs):.3f})  "
            f"AP {np.mean(average_precisions):.3f} ({np.std(average_precisions):.3f})"
        )


def report_proposed_design(X, y, groups, strata_db):
    print("\n" + "=" * 78)
    print("2.5  the proposed design's own interval, and the paired instrument")
    print("=" * 78)
    dropped_feature = "rate4site4mean"
    X_reduced = X.drop(columns=[dropped_feature])

    predictions_a = out_of_fold_predictions(
        X, y, StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=42), (X, strata_db, groups)
    )
    predictions_b = out_of_fold_predictions(
        X_reduced, y, StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=42), (X_reduced, strata_db, groups)
    )
    auc_a = roc_auc_score(y, predictions_a)
    auc_b = roc_auc_score(y, predictions_b)

    low, high = bootstrap_auc_ci(y.values, predictions_a, groups.values)
    print(f"arm A ({X.shape[1]} feat) AUC {auc_a:.4f}  marginal 95% CI {low:.3f}-{high:.3f}  width {high - low:.3f}")
    print(f"arm B ({X_reduced.shape[1]} feat, minus {dropped_feature}) AUC {auc_b:.4f}")
    low, high = bootstrap_auc_difference_ci(y.values, predictions_a, predictions_b, groups.values)
    print(f"paired difference {auc_a - auc_b:+.4f}  95% CI {low:+.4f} to {high:+.4f}  width {high - low:.4f}")


def report_fold_composition():
    print("\n" + "=" * 78)
    print("3.4  fold composition and distinct groups per stratum")
    print("=" * 78)
    for setname in ["set05", "set08"]:
        df = load_set(setname)
        groups = df["acc_db"].map(homology_groups(setname))
        strata = df["db"] + "_" + df["interface"].astype(str)
        print(f"\n{setname}: {groups.nunique()} groups")
        for label, stratum in [("database x label", strata), ("database", df["db"])]:
            per_stratum = pd.DataFrame({"stratum": stratum, "group": groups}).groupby("stratum")["group"].nunique()
            print(f"  groups per stratum ({label}): {dict(per_stratum)}  -> max feasible k = {per_stratum.min()}")
        for n_folds in [5, 10]:
            rows = []
            splitter = StratifiedGroupKFold(n_folds, shuffle=True, random_state=0)
            for _train_idx, test_idx in splitter.split(df, strata, groups=groups):
                fold = df.iloc[test_idx]
                rows.append(
                    {
                        "n_prot": fold["acc_db"].nunique(),
                        "n_pos": int(fold["interface"].sum()),
                        "NMR": int((fold["db"] == "NMR").sum()),
                    }
                )
            summary = pd.DataFrame(rows)
            print(
                f"  k={n_folds}: proteins/fold {summary['n_prot'].tolist()}  "
                f"zero-NMR folds {int((summary['NMR'] == 0).sum())}  "
                f"positives/fold {summary['n_pos'].min()}-{summary['n_pos'].max()}"
            )


def report_partition_diversity(groups, strata_db_label, strata_db):
    print("\n" + "=" * 78)
    print(f"3.5  shuffle behaviour and partition diversity (set05, k={N_FOLDS})")
    print("=" * 78)
    unshuffled = {partition_signature(groups, strata_db_label, seed, shuffle=False) for seed in range(3)}
    print(f"  shuffle=False, 3 seeds            : {len(unshuffled)} distinct partition(s)  <- random_state ignored")
    for label, strata in [("db x label", strata_db_label), ("db only", strata_db)]:
        distinct = {partition_signature(groups, strata, seed) for seed in range(50)}
        print(f"  shuffle=True, {label:10s} 50 seeds : {len(distinct)} distinct partitions")
    distinct_k10 = {partition_signature(groups, strata_db_label, seed, n_folds=10) for seed in range(50)}
    print(f"  shuffle=True, db x label, k=10, 50 seeds : {len(distinct_k10)} distinct partitions")


def report_timing(X_all, y):
    print("\n" + "=" * 78)
    print("6  per-fit timing")
    print("=" * 78)
    for n_estimators in [50, 200]:
        for n_jobs in [1, -1]:
            classifier(n_estimators=n_estimators, n_jobs=n_jobs).fit(X_all.iloc[:100], y.iloc[:100])
            start = time.perf_counter()
            for _ in range(5):
                classifier(n_estimators=n_estimators, n_jobs=n_jobs).fit(X_all, y)
            elapsed = (time.perf_counter() - start) / 5
            print(f"  n_estimators={n_estimators:<4} n_jobs={n_jobs:<3} : {elapsed * 1000:.0f} ms per fit")
    print(f"  in-fold grid = 5 folds x 3 repeats x 576 combinations x 3 inner folds = {5 * 3 * 576 * 3:,} fits")


def report_column_filter(X, X_all, y, groups, strata_db):
    print("\n" + "=" * 78)
    print("4  the column filter is not optional")
    print("=" * 78)
    features = pd.read_csv(MODEL_FEATURES_CSV, index_col=0)
    eligible = set(features.index[features["include"].astype(str).str.upper().isin(["1", "TRUE", "WAHR"])])
    ineligible = [column for column in X_all.columns if column not in eligible]
    print(f"  columns in the naive 'all numeric' matrix but NOT include=TRUE: {ineligible}")

    df = load_set("set05")
    correlation = df[["interface", "interface_score_norm"]].corr().iloc[0, 1]
    print(f"  corr(interface, interface_score_norm) = {correlation:.3f}")

    splitter = StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=0)
    eligible_auc = roc_auc_score(y, out_of_fold_predictions(X, y, splitter, (X, strata_db, groups)))
    splitter = StratifiedGroupKFold(N_FOLDS, shuffle=True, random_state=0)
    naive_auc = roc_auc_score(y, out_of_fold_predictions(X_all, y, splitter, (X_all, strata_db, groups)))
    print(f"  pooled OOF AUC, {X.shape[1]} eligible features : {eligible_auc:.3f}")
    print(f"  pooled OOF AUC, naive all-numeric    : {naive_auc:.3f}")


def main():
    if not RESULTS_DIR.is_dir():
        raise FileNotFoundError(
            f"The results directory is not available at {RESULTS_DIR}. Run 'dvc pull' from the "
            "repository root to fetch the data directory."
        )

    report_data_and_clusters()
    report_blind_test()

    df = load_set("set05")
    y = df["interface"].astype(int)
    X = df[selected_features()]
    groups = df["acc_db"].map(homology_groups("set05"))
    strata_db_label = (df["db"] + "_" + y.astype(str)).values
    strata_db = df["db"].values

    # The "naive" matrix a fold pipeline produces if it drops only the obvious label columns and
    # skips drop_cols_not_used_in_ML. Used to show what that omission costs.
    obvious_label_columns = ["acc_db", "interface", "interface_score", "residue_num", "residue_name"]
    X_all = (
        df.drop(columns=[c for c in obvious_label_columns if c in df.columns]).select_dtypes("number").dropna(axis=1)
    )

    report_splitting_schemes(df, X, y, groups, strata_db_label, strata_db)
    report_proposed_design(X, y, groups, strata_db)
    report_fold_composition()
    report_partition_diversity(groups, strata_db_label, strata_db)
    report_timing(X_all, y)
    report_column_filter(X, X_all, y, groups, strata_db)


if __name__ == "__main__":
    main()
