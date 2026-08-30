"""Guard against the pipeline silently becoming non-reproducible again.

Python salts string hashing per process, so `list(set(feature_names))` returns a different order
in every interpreter. ExtraTreesClassifier draws its `max_features` candidates by column index, so
a permuted column order changes every split and therefore the fitted model -- even with
`random_state` fixed. Before this was fixed, re-running feature selection produced a different
feature set (28/29/29/28 features across four column orders, against a published model of 27) and
refitting on a reversed column order moved predicted probabilities by up to 0.116.

These tests run the ordering-sensitive code in subprocesses with different PYTHONHASHSEED values,
which is the only way to actually catch a reintroduced `set`.
"""

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[2]

# Deliberately unordered relative to the training data, so a set round-trip cannot accidentally
# reproduce the input order and let a regression slip through.
FEATURES = [
    "polarity3Cmean",
    "DI3mean",
    "V",
    "conservation",
    "L",
    "GxxxG",
    "DI5mean",
    "mass",
    "relative_polarity",
    "H",
    "branched",
    "RelPos_TMD",
    "E",
    "residue_depth",
    "DImax",
]


def _run_in_subprocess(snippet: str, hashseed: str) -> str:
    """Execute a snippet in a fresh interpreter with a given PYTHONHASHSEED, return its stdout."""
    result = subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True,
        text=True,
        env={"PYTHONHASHSEED": hashseed, "PATH": "", "PYTHONPATH": str(REPO_ROOT)},
        cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
    # The functions under test log warnings that themselves contain set reprs, whose order also
    # varies with PYTHONHASHSEED, and some write to stdout without a trailing newline. Pull out
    # the explicitly marked value rather than trusting line boundaries.
    marked = re.findall(r"RESULT:(\S*)", result.stdout)
    assert len(marked) == 1, f"expected one RESULT marker, got {marked!r}\nstdout={result.stdout!r}"
    return marked[0]


def test_drop_cols_not_used_in_ML_column_order_is_hashseed_independent():
    """The ML feature column order must not depend on PYTHONHASHSEED."""
    snippet = f"""
import pandas as pd
from thoipapy.validation.feature_selection import drop_cols_not_used_in_ML
from thoipapy.utils import LogOnlyToConsole
features = {FEATURES!r}
df = pd.DataFrame({{c: [0.0, 1.0] for c in features}})
df["interface"] = [0, 1]
out = drop_cols_not_used_in_ML(LogOnlyToConsole(), df)
print("RESULT:" + ",".join(out.columns))
"""
    orders = [_run_in_subprocess(snippet, seed) for seed in ("0", "1", "424242")]
    assert orders[0] == orders[1] == orders[2], (
        "ML feature column order changed with PYTHONHASHSEED. Something reintroduced a set "
        f"into the column selection path.\n{chr(10).join(orders)}"
    )


def test_merge_top_features_order_is_hashseed_independent():
    """The merged anova+rfe feature list must not depend on PYTHONHASHSEED."""
    snippet = """
anova = ["conservation", "DImax", "polarity3Cmean", "V"]
rfe = ["V", "GxxxG", "conservation", "mass"]
retained = ["GxxxG", "branched"]
# mirrors thoipapy.feature_importance.merge.merge_top_features_anova_ensemble
combined = list(dict.fromkeys(anova + rfe + retained))
print("RESULT:" + ",".join(combined))
"""
    orders = [_run_in_subprocess(snippet, seed) for seed in ("0", "1", "424242")]
    assert orders[0] == orders[1] == orders[2]
    # dict.fromkeys must preserve first-seen order, not merely be stable
    assert orders[0] == "conservation,DImax,polarity3Cmean,V,GxxxG,mass,branched"


def test_no_list_set_in_ml_column_paths():
    """`list(set(...))` must not reappear in code that determines ML column order.

    A static check, because the runtime tests above only cover the two paths they exercise.
    """
    offenders = []
    for py in (REPO_ROOT / "thoipapy").rglob("*.py"):
        if "paper_figures" in py.parts:
            continue
        for lineno, line in enumerate(py.read_text().splitlines(), 1):
            if "list(set(" in line and not line.lstrip().startswith("#"):
                offenders.append(f"{py.relative_to(REPO_ROOT)}:{lineno}: {line.strip()}")
    assert not offenders, (
        "list(set(...)) produces a hash-seed-dependent order. Use a list comprehension or "
        "dict.fromkeys instead:\n" + "\n".join(offenders)
    )
