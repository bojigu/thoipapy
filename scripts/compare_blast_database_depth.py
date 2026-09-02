"""What does THOIPA lose by searching UniRef90 instead of NCBI nr?

    THOIPA_UNIREF_DATA=/path/to/uniref90/data python scripts/compare_blast_database_depth.py

The published model, and every number in the paper, rests on alignments from NCBI ``nr``. The
standalone predictor and the webserver search a local UniRef90 database instead, because a remote
nr query has been measured queueing for hours and a webserver job cannot wait. UniRef90 was assumed to be
shallower by construction: it clusters UniProt at 90% identity, so ~121 million clusters against
nr's billion-plus records.

Nobody had measured what that costs, and the assumption turns out to be wrong on depth. The
comparison is unavoidably confounded with time -- the nr archive dates from May 2020, a local
UniRef90 search runs against today's release -- and six years of sequence growth outweighs the
clustering. This script does, on three axes:

1. **Alignment depth.** The number of unique TMD homologues that survive filtering, per protein.
   This is the raw input every conservation and coevolution feature is computed from.
2. **Predictive accuracy.** Leave-one-protein-out over set08, and the blind set07 test set,
   scored by identical code on both feature sets.
3. **Which features degrade.** Per-feature correlation between the two datasets, so the loss can
   be attributed rather than just totalled.

Both feature sets are read, never written. The UniRef90 rebuild lives outside the repository and
is deliberately not DVC-tracked: it is the weaker dataset, and the published nr artefacts must
stay the reference.
"""

import os
import sys
from ast import literal_eval
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT))

NR_DATA = REPO_ROOT / "data"
_UNIREF_DATA_ENV_VAR = "THOIPA_UNIREF_DATA"


def _comparison_data_dir() -> Path:
    """The data directory to compare against, from the environment.

    No default. The directory being compared is whatever the caller rebuilt, and its location is a
    property of their machine, not of this repository.
    """
    value = os.environ.get(_UNIREF_DATA_ENV_VAR, "").strip()
    if not value:
        raise SystemExit(
            f"{_UNIREF_DATA_ENV_VAR} is not set. Point it at a data directory built from a "
            "different homologue source, for example the output of a pipeline run with "
            "THOIPA_LOCAL_BLAST_DB set to a local UniRef90 database."
        )
    return Path(value)


TRAIN_SET, TEST_SET = "set08", "set07"
BIND_COLUMN = "interface"
N_ESTIMATORS = 1000
FOREST_SEEDS = (0, 1, 2, 3, 4)


def training_table(data_dir, setname):
    return pd.read_csv(data_dir / f"results/{setname}/train_data/01_train_data_orig.csv", index_col=0)


def selected_features(data_dir, setname):
    df = pd.read_csv(data_dir / f"results/{setname}/train_data/03_train_data_after_first_feature_seln.csv", index_col=0)
    return [c for c in df.columns if c != BIND_COLUMN]


def homologue_clusters(data_dir, setname):
    clusters = pd.read_excel(
        data_dir / f"results/{setname}/clusters/{setname}_sim_matrix.xlsx",
        sheet_name="reduced_clusters",
        index_col=0,
    )
    members = (
        clusters["reduced_clusters"]
        .apply(literal_eval)
        .apply(lambda names: ["-".join(n.split("-")[1:]) for n in list(names)])
    )
    return {acc: set(c) for c in members for acc in c}, list(members)


def forest(tuned_csv, seed):
    params = pd.read_csv(tuned_csv, index_col=0)["GridSearchSlowMethod"]
    return ExtraTreesClassifier(
        n_estimators=N_ESTIMATORS,
        criterion=params["criterion"],
        min_samples_leaf=int(params["min_samples_leaf"]),
        max_depth=int(params["max_depth"]),
        max_features="sqrt",
        bootstrap=False,
        n_jobs=-1,
        random_state=seed,
    )


def leave_one_out(data_dir, setname):
    """Per-protein ROC and PR AUC, leaving out each protein and its whole CD-HIT cluster."""
    df = training_table(data_dir, setname)
    features = selected_features(data_dir, setname)
    tuned = data_dir / f"results/{setname}/train_data/04_tuned_ensemble_parameters.csv"
    y, groups = df[BIND_COLUMN], df["acc_db"]
    homologues, clusters = homologue_clusters(data_dir, setname)

    roc, pr = {}, {}
    for acc_db in groups.unique():
        train = (~groups.isin(homologues.get(acc_db, {acc_db}))).to_numpy()
        test = (groups == acc_db).to_numpy()
        if test.sum() == 0 or y[test].nunique() < 2:
            continue
        X = df[features]
        r, p = [], []
        for seed in FOREST_SEEDS:
            fitted = forest(tuned, seed).fit(X[train], y[train])
            prediction = fitted.predict_proba(X[test])[:, 1]
            r.append(roc_auc_score(y[test], prediction))
            p.append(average_precision_score(y[test], prediction))
        roc[acc_db], pr[acc_db] = float(np.mean(r)), float(np.mean(p))
    return pd.Series(roc), pd.Series(pr), clusters


def train_on_trainset_test_on_testset(data_dir):
    """Train on set08, score every protein of the blind set07. Truly held out: disjoint sets."""
    train_df = training_table(data_dir, TRAIN_SET)
    test_df = training_table(data_dir, TEST_SET)
    features = selected_features(data_dir, TRAIN_SET)
    missing = [f for f in features if f not in test_df.columns]
    if missing:
        raise ValueError(f"{TEST_SET} feature table is missing {missing}")
    tuned = data_dir / f"results/{TRAIN_SET}/train_data/04_tuned_ensemble_parameters.csv"

    # Fit once per seed, not once per seed per protein: the training set is the same every time.
    fitted_forests = [forest(tuned, seed).fit(train_df[features], train_df[BIND_COLUMN]) for seed in FOREST_SEEDS]

    roc, pr = {}, {}
    for acc_db, rows in test_df.groupby("acc_db"):
        if rows[BIND_COLUMN].nunique() < 2:
            continue
        r, p = [], []
        for fitted in fitted_forests:
            prediction = fitted.predict_proba(rows[features])[:, 1]
            r.append(roc_auc_score(rows[BIND_COLUMN], prediction))
            p.append(average_precision_score(rows[BIND_COLUMN], prediction))
        roc[acc_db], pr[acc_db] = float(np.mean(r)), float(np.mean(p))
    return pd.Series(roc), pd.Series(pr)


def cluster_means(differences, clusters):
    values = [
        differences[[a for a in cluster if a in differences.index]].mean()
        for cluster in clusters
        if any(a in differences.index for a in cluster)
    ]
    return np.array([v for v in values if not np.isnan(v)])


def bootstrap_ci(by_cluster, n_draws=10000, seed=0):
    rng = np.random.default_rng(seed)
    draws = rng.choice(by_cluster, size=(n_draws, len(by_cluster)), replace=True).mean(axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def alignment_depth(data_dir, setname):
    """Unique TMD homologues per protein, recovered from the cube root stored as n_homologues."""
    df = training_table(data_dir, setname)
    per_protein = df.groupby("acc_db")["n_homologues"].first()
    return (per_protein**3).round().astype(int)


def report_depth(UNIREF_DATA):
    print("=" * 100)
    print("1. ALIGNMENT DEPTH: unique TMD homologues surviving the identity and gap filters")
    print("=" * 100)
    rows = []
    for setname in (TRAIN_SET, TEST_SET):
        nr = alignment_depth(NR_DATA, setname)
        uniref = alignment_depth(UNIREF_DATA, setname)
        common = nr.index.intersection(uniref.index)
        nr, uniref = nr[common], uniref[common]
        ratio = uniref / nr.replace(0, np.nan)
        rows.append(
            {
                "set": setname,
                "n": len(common),
                "nr median": int(nr.median()),
                "UniRef90 median": int(uniref.median()),
                "median ratio": f"{ratio.median():.2f}",
                "deeper on UniRef90": int((uniref > nr).sum()),
                "shallower": int((uniref < nr).sum()),
            }
        )
        if setname == TRAIN_SET:
            worst = ratio.nsmallest(5)
            best = ratio.nlargest(3)
            depth_detail = (nr, uniref, worst, best)
    print(pd.DataFrame(rows).to_string(index=False))
    nr, uniref, worst, best = depth_detail
    print(
        f"\n  {TRAIN_SET} totals: nr {int(nr.sum()):,} homologues, UniRef90 {int(uniref.sum()):,} "
        f"({uniref.sum() / nr.sum():.0%} of nr)"
    )
    print("  worst-hit proteins (UniRef90/nr):  " + ", ".join(f"{a} {v:.2f}" for a, v in worst.items()))
    print("  best-served proteins:              " + ", ".join(f"{a} {v:.2f}" for a, v in best.items()))
    return nr, uniref


def report_accuracy(UNIREF_DATA):
    print("\n" + "=" * 100)
    print("2. PREDICTIVE ACCURACY: identical code, identical folds, only the alignment source differs")
    print("=" * 100)
    results = {}
    for label, data_dir in [("nr (published)", NR_DATA), ("UniRef90 (server)", UNIREF_DATA)]:
        roc, pr, clusters = leave_one_out(data_dir, TRAIN_SET)
        results[label] = (roc, pr, clusters)
        print(
            f"  {label:<20} set08 leave-one-out   ROC {roc.mean():.3f}   PR {pr.mean():.3f}   " f"({len(roc)} proteins)"
        )

    test = {}
    for label, data_dir in [("nr (published)", NR_DATA), ("UniRef90 (server)", UNIREF_DATA)]:
        try:
            roc, pr = train_on_trainset_test_on_testset(data_dir)
            test[label] = (roc, pr)
            print(
                f"  {label:<20} set07 blind test      ROC {roc.mean():.3f}   PR {pr.mean():.3f}   "
                f"({len(roc)} proteins)"
            )
        except (FileNotFoundError, ValueError) as e:
            print(f"  {label:<20} set07 blind test      unavailable: {e}")

    nr_roc, _, clusters = results["nr (published)"]
    uni_roc, _, _ = results["UniRef90 (server)"]
    common = nr_roc.index.intersection(uni_roc.index)
    differences = uni_roc[common] - nr_roc[common]
    by_cluster = cluster_means(differences, clusters)
    low, high = bootstrap_ci(by_cluster)
    p = wilcoxon(by_cluster).pvalue
    print(
        f"\n  set08 leave-one-out, UniRef90 minus nr: {by_cluster.mean():+.3f} AUC "
        f"[95% CI {low:+.3f}, {high:+.3f}], Wilcoxon p = {p:.3f}"
    )
    print(f"  ({len(by_cluster)} CD-HIT clusters, not {len(common)} proteins: several are homologues)")
    if len(test) == 2:
        nr_t, uni_t = test["nr (published)"][0], test["UniRef90 (server)"][0]
        shared = nr_t.index.intersection(uni_t.index)
        print(
            f"  set07 blind test,  UniRef90 minus nr: {(uni_t[shared] - nr_t[shared]).mean():+.3f} AUC "
            f"({len(shared)} proteins, descriptive only)"
        )
    return results


def report_feature_agreement(UNIREF_DATA):
    print("\n" + "=" * 100)
    print("3. WHICH FEATURES DEGRADE: Spearman correlation between the two datasets, per feature")
    print("=" * 100)
    nr = training_table(NR_DATA, TRAIN_SET)
    uni = training_table(UNIREF_DATA, TRAIN_SET)
    rows = nr.index.intersection(uni.index)
    features = [f for f in selected_features(NR_DATA, TRAIN_SET) if f in uni.columns]
    correlations = {
        f: spearmanr(nr.loc[rows, f], uni.loc[rows, f]).statistic
        for f in features
        if nr.loc[rows, f].nunique() > 1 and uni.loc[rows, f].nunique() > 1
    }
    ranked = pd.Series(correlations).sort_values()
    print(f"  {len(rows)} shared residues, {len(ranked)} features present in both\n")
    print("  least reproducible across databases:")
    for name, rho in ranked.head(8).items():
        print(f"    {name:<28} rho {rho:+.3f}")
    print("\n  most reproducible:")
    for name, rho in ranked.tail(5).items():
        print(f"    {name:<28} rho {rho:+.3f}")
    return ranked


if __name__ == "__main__":
    UNIREF_DATA = _comparison_data_dir()
    if not UNIREF_DATA.is_dir():
        raise SystemExit(f"comparison data directory not found: {UNIREF_DATA}")
    report_depth(UNIREF_DATA)
    report_accuracy(UNIREF_DATA)
    report_feature_agreement(UNIREF_DATA)
    print(
        "\nBoth datasets were read, neither written. The UniRef90 rebuild is intentionally not "
        "DVC-tracked:\nit is the weaker dataset, and the published nr artefacts remain the reference."
    )
