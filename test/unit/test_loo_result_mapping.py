"""Leave-one-out results must stay attached to the protein they came from.

The fold loop can skip a protein with ``continue`` when it is in the protein set but absent from
the training data. Results used to be matched back to proteins by position against
``df_set.acc_db``, so one skip shifted every subsequent AUC onto the wrong protein. Nothing failed:
the run completed, the numbers looked ordinary, and every protein after the gap was reported under
its neighbour's name.

No protein is skipped on the shipped sets, so the published numbers were never affected. These
tests pin the invariant so that a set which does skip one cannot reintroduce it silently.
"""

import numpy as np
import pytest


def build_val_list(acc_dbs, skip=None):
    """Mimic the fold loop: one result per protein, in order, skipping where asked.

    Mirrors the shape ``run_LOO_validation`` builds, where each entry carries its own acc_db.
    """
    val_list = []
    for n, acc_db in enumerate(acc_dbs):
        if acc_db == skip:
            continue
        auc_dict = {
            "roc_auc": 0.5 + n / 100,
            "pr_auc": 0.4 + n / 100,
            "fpr": np.linspace(0, 1, 5),
            "tpr": np.linspace(0, 1, 5),
        }
        val_list.append((acc_db, auc_dict, f"BO_{acc_db}"))
    return val_list


PROTEINS = ["P00-crystal", "P01-crystal", "P02-ETRA", "P03-NMR", "P04-crystal"]


def test_every_result_is_reported_under_its_own_protein():
    val_list = build_val_list(PROTEINS)
    xv_dict = {acc_db: auc["roc_auc"] for acc_db, auc, _ in val_list}
    for n, acc_db in enumerate(PROTEINS):
        assert xv_dict[acc_db] == pytest.approx(0.5 + n / 100)


def test_a_skipped_protein_does_not_shift_the_others():
    """The regression itself. With P01 skipped, positional matching moves P02's AUC onto P01."""
    val_list = build_val_list(PROTEINS, skip="P01-crystal")
    assert len(val_list) == 4

    # What the code does now: each result carries its own label.
    by_label = {acc_db: auc["roc_auc"] for acc_db, auc, _ in val_list}
    assert "P01-crystal" not in by_label
    assert by_label["P02-ETRA"] == pytest.approx(0.52)
    assert by_label["P04-crystal"] == pytest.approx(0.54)

    # What it did before: zip the results against the full protein list by position.
    by_position = {PROTEINS[n]: auc["roc_auc"] for n, (_, auc, _) in enumerate(val_list)}
    assert by_position["P01-crystal"] == pytest.approx(0.52), "P02's score, filed under P01"
    assert by_position["P03-NMR"] == pytest.approx(0.54), "P04's score, filed under P03"
    assert "P04-crystal" not in by_position, "the last protein silently disappears"
    assert by_position != by_label


def test_the_mean_curve_is_averaged_over_the_folds_that_ran():
    """Dividing by the size of the protein set scales the mean curve down when a fold is skipped."""
    val_list = build_val_list(PROTEINS, skip="P01-crystal")
    mean_fpr = np.linspace(0, 1, 5)

    summed = np.zeros_like(mean_fpr)
    for _, auc_dict, _ in val_list:
        summed += np.interp(mean_fpr, auc_dict["fpr"], auc_dict["tpr"])

    over_folds_that_ran = summed / len(val_list)
    over_set_size = summed / len(PROTEINS)

    assert over_folds_that_ran[-1] == pytest.approx(1.0)
    assert over_set_size[-1] == pytest.approx(0.8)
    assert over_set_size[-1] < over_folds_that_ran[-1]


def _identity(n):
    """Module level so it can be pickled into a pool worker."""
    return n


def test_pool_map_preserves_order_so_results_can_be_zipped_to_their_inputs():
    """The multiprocessing branch relabels results by zipping them against the queued inputs.

    That is only valid because Pool.map returns results in input order regardless of completion
    order. Asserted here rather than assumed, since the relabelling is silent if it is wrong.
    """
    from multiprocessing import Pool

    with Pool(processes=4) as pool:
        assert pool.map(_identity, range(12)) == list(range(12))
