# Which homologue source should THOIPA use?

## The question

Every conservation and coevolution feature THOIPA uses is computed from an alignment of TMD homologues, so the homologue search bounds what the model can learn. The published model used a single blastp pass against NCBI `nr` in May 2020. That source is now unavailable: NCBI does not archive old `nr` releases, so the snapshot behind the published numbers cannot be recovered, and a live remote `nr` query has been measured queueing for hours (a probe of a three-residue query returned an NCBI estimate of 3.1 hours). `nr` is the reference, not a candidate.

Two candidates remain, and they differ mainly in what they cost:

- **UniRef90 local** — a single blastp pass against a local database. 96 GB of server disk, rebuilt every UniProt release, depends on nothing outside the machine.
- **ColabFold MSA server** — MMseqs2 with three profile iterations against UniRef30 plus an environmental database (BFD, MGnify, MetaEuk, SMAG), with matched clusters expanded back to their members. No disk, no maintenance, but a dependency on a free academic service.

The specific hypothesis tested was that the greater evolutionary depth available in 2026 would improve prediction accuracy over the 2020 archive.

## Method

Both candidates were rebuilt into parallel data directories: every alignment-derived feature recomputed, feature selection and hyperparameter tuning rerun, nothing else changed. Everything alignment-independent was copied from the tracked directory, including the CD-HIT clusters that define the cross-validation folds, so the folds are identical across sources. Each rebuilt training table has the same shape as the published one (850 residues for set08, 194 for set07, 115 columns) with identical interface labels, so only the feature values differ.

**Every set08 accuracy figure below is weighted per CD-HIT cluster, not per protein.** Several proteins in set08 are homologues of each other, and treating them as independent counts the same evidence more than once. 40 proteins reduce to 34 clusters. Leave-one-out excludes the held-out protein *and its whole cluster* from training. The unit matters: on per-protein means ColabFold sits above UniRef90, on cluster means it sits below, and only the latter is consistent with the paired differences.

Reproduce with `scripts/compare_colabfold_alignment_depth.py` (downloads alignments, measures depth), `scripts/build_colabfold_dataset.py` (rebuilds features and retrains), `scripts/compare_homologue_sources.py` (scores the three sources), and `scripts/compare_alignment_depth_control.py` (the depth control below).

## Alignment depth

Unique TMD homologues surviving the identity and gap filters.

| source | set08 median | set08 total | set07 median | set07 total |
|---|---:|---:|---:|---:|
| nr 2020 | 235 | 26,672 | 734 | 8,640 |
| UniRef90 local | 348 | 37,996 | 517 | 8,580 |
| ColabFold | 703 | 233,836 | 2,104 | 69,727 |

Ratios, stated both ways because they differ by up to 3x and are easy to conflate:

| comparison | median of per-protein ratios | ratio of set totals | deeper on |
|---|---:|---:|---|
| UniRef90 / nr | 1.44x | 1.42x | 30 of 40 |
| ColabFold / nr | 2.63x | 8.77x | 36 of 40 |
| ColabFold / UniRef90 | 2.17x | 6.15x | 34 of 40 |

The total ratio is dominated by a handful of proteins whose alignments grew enormously; the median is what a typical protein gains. UniRef90 is *shallower* than the 2020 `nr` archive on set07 — the local database was adopted to escape the NCBI queue, not because it searches better.

## Predictive accuracy

| source | set08 LOO ROC | set08 LOO PR | set07 blind ROC | set07 blind PR |
|---|---:|---:|---:|---:|
| nr 2020 | 0.656 | 0.545 | 0.654 | 0.616 |
| UniRef90 local | 0.643 | 0.534 | 0.635 | 0.608 |
| ColabFold | 0.639 | 0.544 | 0.685 | 0.628 |

The set08 columns are cluster-weighted. The set07 columns are per-protein means over 10 proteins with no confidence interval, and are descriptive only.

Paired differences over the 34 clusters, set08 leave-one-out ROC AUC:

| comparison | difference | 95% CI | Wilcoxon p |
|---|---:|---|---:|
| ColabFold minus UniRef90 | -0.005 | [-0.026, +0.016] | 0.78 |
| UniRef90 minus nr 2020 | -0.012 | [-0.034, +0.007] | 0.30 |
| ColabFold minus nr 2020 | -0.017 | [-0.049, +0.012] | 0.57 |

Every interval contains zero. Rather than reading that as "no difference", which a failure to reject does not establish, the useful statement is what the intervals exclude: **any ColabFold advantage over UniRef90 larger than +0.016 AUC, or over nr larger than +0.012 AUC, is ruled out at 95%.** An improvement smaller than that would not be detectable here — the minimum detectable effect at 80% power is 0.031 AUC against UniRef90 and 0.044 against nr.

On set07 ColabFold leads UniRef90 by +0.050 ROC, but that is 10 proteins, has no confidence interval, and disagrees in sign with the far better powered leave-one-out. It should not carry weight.

## Does depth alone change anything?

Every comparison above is confounded. `nr` differs from a 2026 search by six years of sequence growth as well as by depth; UniRef90 differs from ColabFold by search algorithm as well as by depth. Neither separates "more sequences" from "different sequences".

`scripts/compare_alignment_depth_control.py` does. Both arms are the same server, the same three-iteration MMseqs2 search, the same query and the same downstream filters. The only difference is `colabfold_mode`: with the diversity filter on, MMseqs2 thins near-identical hits with `--qsc 0.8 --max-seq-id 0.95`; with it off they are kept. The feature set and the hyperparameters are pinned to the published `nr` selection for both arms, so not even feature selection varies.

| | median depth | total depth |
|---|---:|---:|
| ColabFold `env` | 472 | 90,964 |
| ColabFold `env-nofilter` | 703 | 233,836 |
| manipulation | 1.37x | 2.57x |

| metric | deeper minus shallower | 95% CI | Wilcoxon p |
|---|---:|---|---:|
| ROC AUC | **+0.0001** | [-0.0156, +0.0157] | 0.95 |
| PR AUC | +0.0071 | [-0.0108, +0.0266] | 0.55 |

A 2.57x increase in total alignment depth, with everything else held constant, moves ROC AUC by one ten-thousandth. This is the cleanest evidence in the analysis, and unlike the source comparisons it carries no time or algorithm confound.

## What the alignment source does change

Spearman correlation of each selected feature against its published `nr` value, same residues:

| source | median rho | coevolution features median rho | worst |
|---|---:|---:|---|
| UniRef90 local | +0.974 | +0.525 | `DImax` +0.452 |
| ColabFold | +0.898 | +0.303 | `DImax` +0.256 |

Sequence-derived features (`GxxxG`, `SmxxxSm`, `mass`, `branched`) reproduce exactly, as they must. The FreeContact coevolution features do not: `DImax`, one of the features listed in `features_to_be_retained_during_selection`, correlates at rho 0.26 between the `nr` and ColabFold datasets.

So the homologue source substantially rewrites the coevolution signal while leaving accuracy unchanged. That is a question about the coevolution block rather than about the alignment source, and it is probably where the next accuracy work belongs. Two caveats weaken the section, and both should be read before acting on it:

- **One of the "coevolution" features carries no alignment information at all.** `thoipapy/features/freecontact.py` initialises `highest_XI_face_value = 0` and never updates it inside the face loop, so the best-face test is always `> 0` and the winner is whichever face was evaluated last. `MI_highest_face` is therefore a fixed heptad mask, byte-identical across all four datasets, and it is among the features selected for the ColabFold model. Pre-existing, not introduced by the homologue work, and worth fixing before anyone draws conclusions about coevolution stability.
- **The feature-selection churn is not attributed.** 26 features selected on `nr`, 29 on ColabFold, 19 in common, with `DI3mean`, `DI5mean`, `MIall_mean` and `rate4site4mean` dropping out, and also `LIPS_polarity`, `LIPS_surface_ranked_norm` and `N`. With `n_top_features_to_keep = 20` for both ANOVA and RFE, the union of two top-20 lists is intrinsically unstable, and no control establishes how much it moves under a perturbation that is *not* a source change. The churn is real; blaming it on the source is not established.

## Exploratory: is very deep better or worse?

Post-hoc, computed per protein rather than per cluster, and run after seeing the per-protein differences. Hypothesis-generating only.

Splitting set08 by the depth of the new alignment, the deepest quartile (median around 15,000 homologues) is the only one with a mean accuracy loss, and the rank correlation between depth gained and accuracy gained is negative. Nothing reaches p < 0.05, the correlations weaken further when computed per cluster, and proteins that scored badly on `nr` tend to improve on any other source, which is regression to the mean rather than a source effect. By structure type the pattern is not monotone: crystal proteins gained the most depth and lost the most accuracy, but NMR proteins gained intermediate depth and gained the most accuracy.

The mechanism this would suggest — a very deep alignment diluting family-specific covariation with distant environmental homologues — is plausible, given that the downstream filter is `frac_ident_TMD > 0.2` over a roughly 20-residue TMD, which is five matching residues. It is not established.

## Cost

| | disk | maintenance | external dependency |
|---|---|---|---|
| UniRef90 local | 96 GB | rebuild every UniProt release | none |
| ColabFold, public server | none | none | free academic service |
| ColabFold, self-hosted | 103 GB + 118 GB as downloads, nearer 1 TB built with indexes | database updates | none |

A ColabFold search typically takes 10 to 70 seconds per protein, against hours for a remote `nr` query; over the 50 archives downloaded here the median was 32 s and the slowest 455 s. Rebuilding all features and retraining took 39 minutes for set08 and 10 minutes for set07.

The public server processes a few thousand alignments a day and asks for serial queries from a single IP. A 40-protein set is three orders of magnitude inside that, and the webserver has never carried enough traffic to come near it, so the public server is a legitimate default for both. `THOIPA_COLABFOLD_HOST` points the same protocol at a self-hosted MMseqs2 server if that changes — but self-hosting costs roughly an order of magnitude more disk than the UniRef90 database it would replace, so it buys independence, not space.

The archive returned by the server includes `msa.sh`, the exact MMseqs2 command line and database versions used. The BLAST path had no equivalent, so an alignment's provenance was previously unrecoverable.

## Licensing

Nothing in this stack restricts commercial use, and none of it is academic-only:

| component | licence | commercial use |
|---|---|---|
| ColabFold | MIT | permitted |
| MMseqs2 | MIT | permitted |
| ColabFold databases (UniRef30, ColabFoldDB, BFD/MGnify as distributed) | CC BY 4.0 | permitted with attribution |
| UniProt / UniRef upstream | CC BY 4.0 | permitted with attribution |

The attribution condition is the only obligation, and it is met by citing ColabFold, MMseqs2 and
the databases in anything published from THOIPA output.

The public server is a different matter, and it is not a licensing one. `api.colabfold.com`
publishes no terms of service — `/`, `/terms` and `/license` all return 404, and neither the
ColabFold README nor its notebooks carry an academic-only clause. What exists is a fair-use
expectation: the compute is donated by KOBIC and the Söding Lab, the operators describe it as a
limited shared resource of a few thousand alignments a day, ask for serial queries from a single
IP, and reserve the right to limit access case by case. There is correspondingly no service
guarantee, and the server was overwhelmed once already, in August 2025, by a single user
submitting a large batch.

So the exposure is not legal. It is that a courtesy can be withdrawn, without notice and without
recourse, from a service with no contract behind it. The mitigations are the ones already in the
client: serial submission, a user agent carrying a contact address so the operators can ask before
they block, and `THOIPA_COLABFOLD_HOST` to move to a self-hosted MMseqs2 server if the volume ever
justifies it.

## Conclusions

1. ColabFold and UniRef90 are indistinguishable on accuracy (-0.005 ROC, with any ColabFold advantage above +0.016 excluded). The 96 GB buys nothing.
2. Neither beats the 2020 `nr` archive on the cluster-weighted leave-one-out, and both differences sit inside their intervals.
3. **Depth alone does not improve THOIPA.** A 2.57x depth manipulation of a single search, with everything else held constant, moves accuracy by +0.0001 AUC.
4. This does **not** show that no better alignment exists. Outside the control, depth and alignment composition were varied together; inside it, the manipulation is a diversity filter rather than a different search. A deeper alignment built by a different search, or filtered more carefully for the TM context, has not been tested. The defensible claim is that these homologue sources are interchangeable, and that adding sequences by relaxing a diversity filter does not help.
5. For the training pipeline, ColabFold via the public server is the better default: comparable accuracy, deeper alignments, no disk, no maintenance, recorded provenance, and no NCBI queue.
6. If accuracy is the goal, the evidence points at the coevolution features rather than at the homologue search — starting with the `MI_highest_face` bug above.

## Limits of this analysis

**Power.** 34 clusters detects roughly a 0.03 to 0.05 AUC effect. A smaller real improvement is invisible here.

**These are not clean generalisation estimates.** Feature selection (`03_train_data_after_first_feature_seln.csv`) and hyperparameters (`04_tuned_ensemble_parameters.csv`) are fit on all 40 proteins and then reused inside every leave-one-out fold. The absolute AUCs are optimistic, equally so for every arm. The paired *differences* are the trustworthy part. The depth control avoids this for its own comparison by pinning both arms to one model.

**Retuning per source is deliberate** in the three-source comparison: each arm gets its own feature selection and hyperparameters, which is the realistic deployment comparison rather than an isolation of the alignment. `scripts/compare_blast_database_depth.py` previously hardcoded `max_features="sqrt"` while reading every other parameter from the tuned file, so UniRef90 — which tunes to `log2` — was scored with a parameter it had not been tuned for. Fixed; it moves UniRef90 by +0.002 AUC and the paired differences by at most 0.003.

**The arms were not all produced by the same code version.** `residue_depth` is geometric and alignment-independent, yet it differs between the published table and both rebuilds on 23 of 850 rows by ±0.1, while the three rebuilds agree with each other exactly. The cause is a float round-trip through CSV in `thoipapy/features/relative_position.py` that flips a rounding tie. Immaterial in size, but "identical code, only the feature values differ" is not strictly true of the `nr` arm.

**The e-value cutoff means different things across sources.** `e_value_cutoff = 100` is numerically identical everywhere, but an e-value from a three-iteration MMseqs2 profile search is not the same quantity as one from a single blastp pass. The filter is nominally constant and effectively is not.

**The a3m path keeps hits that BLAST dropped.** BLAST returned a partial query span per hit, so a homologue that did not cover the whole TMD failed the TMD regex and was discarded. An a3m is query-anchored, so the same homologue arrives padded to full length and survives if its TMD slice has few enough gaps. On 1xioA4 this is 6 of 390 retained rows. It is a filter-semantics difference between sources, not a conversion error, and it applies equally to both ColabFold arms.

**The comparison against nr is confounded with time**, and because NCBI does not archive old releases, the unconfounded version is not available to anyone.

**set07 is not fully independent of set08.** `5nkqA3` (set08) and `5nkqA1` (set07) are different chains of the same PDB entry. Accessions are disjoint, structures are not, and set07 received no CD-HIT treatment. This changes no conclusion, since set07 is already treated as descriptive.

**Both rebuilds are read-only with respect to the repository.** They live outside it and are not DVC-tracked; the published `nr` artefacts remain the reference.
