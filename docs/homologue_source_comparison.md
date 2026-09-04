# Which homologue source should THOIPA use?

## The question

Every conservation and coevolution feature THOIPA uses is computed from an alignment of TMD homologues, so the homologue search bounds what the model can learn. The published model used a single blastp pass against NCBI `nr` in May 2020. That source is now unavailable: NCBI does not archive old `nr` releases, so the snapshot behind the published numbers cannot be recovered, and a live remote `nr` query has been measured queueing for hours (a probe of a three-residue query returned an NCBI estimate of 3.1 hours). `nr` is the reference, not a candidate.

Two candidates remain, and they differ mainly in what they cost:

- **UniRef90 local** — a single blastp pass against a local database. 96 GB of server disk, rebuilt every UniProt release, depends on nothing outside the machine.
- **ColabFold MSA server** — MMseqs2 with three profile iterations against UniRef30 plus an environmental database (BFD, MGnify, MetaEuk, SMAG), with matched clusters expanded back to their members. No disk, no maintenance, but a dependency on a free academic service.

The specific hypothesis tested was that the greater evolutionary depth available in 2026 would improve prediction accuracy over the 2020 archive.

## Method

Both candidates were rebuilt into parallel data directories: every alignment-derived feature recomputed, feature selection and hyperparameter tuning rerun, nothing else changed. Everything alignment-independent was copied from the tracked directory, including the CD-HIT clusters that define the cross-validation folds, so the folds are identical across sources. Each rebuilt training table has the same shape as the published one (850 residues for set08, 194 for set07, 115 columns) with identical interface labels, so only the feature values differ.

Accuracy is reported per CD-HIT cluster rather than per protein. Several proteins in set08 are homologues of each other, and treating them as independent would count the same evidence twice and make the confidence interval too narrow. 40 proteins reduce to 34 clusters.

Reproduce with `scripts/compare_colabfold_alignment_depth.py` (downloads alignments, measures depth), then `scripts/build_colabfold_dataset.py` (rebuilds features and retrains), then `scripts/compare_homologue_sources.py` (scores all three).

## Alignment depth

Unique TMD homologues surviving the identity and gap filters.

| source | set08 median | set08 total | set07 median | set07 total |
|---|---:|---:|---:|---:|
| nr 2020 | 235 | 26,672 | 734 | 8,640 |
| UniRef90 local | 348 | 37,996 | 517 | 8,580 |
| ColabFold | 703 | 233,836 | 2,104 | 69,727 |

ColabFold is 2.17x deeper than UniRef90 on the median set08 protein and deeper on 34 of 40. Note that UniRef90 is *shallower* than the 2020 `nr` archive on set07: the local database was adopted to escape the NCBI queue, not because it searched better.

## Predictive accuracy

| source | set08 LOO ROC | set08 LOO PR | set07 blind ROC | set07 blind PR |
|---|---:|---:|---:|---:|
| nr 2020 | 0.641 | 0.534 | 0.654 | 0.616 |
| UniRef90 local | 0.630 | 0.513 | 0.639 | 0.607 |
| ColabFold | 0.636 | 0.540 | 0.685 | 0.628 |

Paired differences over the 34 CD-HIT clusters, set08 leave-one-out ROC AUC:

| comparison | difference | 95% CI | Wilcoxon p |
|---|---:|---|---:|
| ColabFold minus UniRef90 | -0.002 | [-0.024, +0.020] | 0.80 |
| UniRef90 minus nr 2020 | -0.015 | [-0.034, +0.002] | 0.24 |
| ColabFold minus nr 2020 | -0.017 | [-0.049, +0.012] | 0.57 |

On set07 ColabFold beats UniRef90 by +0.046 ROC, but that is 10 proteins with no confidence interval and it disagrees in sign with the far better powered leave-one-out. It is descriptive only.

The same comparison with the MMseqs2 diversity filter left on (`env` rather than `env-nofilter`, 2.3x deeper than `nr` instead of 8.8x) gives -0.004 ROC [-0.029, +0.021], p = 0.79. Two independent runs at different depths, both flat.

**Deeper alignments did not improve accuracy.** Roughly nine times more homologues in total changed nothing measurable.

## More depth is, if anything, mildly harmful

| set08 quartile by final depth | median depth | mean change in ROC AUC |
|---|---:|---:|
| Q1 shallowest | 119 | +0.022 |
| Q2 | 393 | -0.002 |
| Q3 | 1,691 | +0.032 |
| Q4 deepest | 14,931 | -0.072 |

The loss concentrates entirely in the deepest quartile, and correspondingly by structure type: crystal proteins gained the most depth (7.7x median) and lost the most accuracy (-0.048), while ETRA proteins gained least (1.8x) and were unchanged (+0.007). The rank correlation between depth gained and accuracy gained is -0.25 (p = 0.12).

This is consistent with a very deep alignment diluting a family-specific covariation signal with distant environmental homologues, but it is not established: no correlation here reaches p < 0.05, and the analysis was run after seeing the per-protein differences.

## The coevolution features are the unstable part

Spearman correlation of each selected feature against its published `nr` value, same residues:

| source | median rho | coevolution features median rho |
|---|---:|---:|
| UniRef90 local | +0.974 | +0.525 |
| ColabFold | +0.898 | +0.303 |

The sequence-derived features (`GxxxG`, `SmxxxSm`, `mass`, `branched`) reproduce exactly, as they must. The FreeContact coevolution features do not: `DImax` correlates at rho +0.26 between the `nr` and ColabFold datasets, `DI5mean` at +0.30, `DI3mean` at +0.30. `DImax` is one of the features listed in `features_to_be_retained_during_selection`.

Feature selection destabilised with it: 26 features selected on `nr`, 29 on ColabFold, 19 in common, with `DI3mean`, `DI5mean`, `MIall_mean` and `rate4site4mean` dropping out entirely.

So changing the homologue source substantially rewrites the coevolution signal while leaving accuracy unchanged. The most economical reading is that those features carry less stable information than their prominence in the model suggests. That is a question about the coevolution block rather than about the alignment source, and it is where the next accuracy work should probably go.

## Cost

| | disk | maintenance | external dependency |
|---|---|---|---|
| UniRef90 local | 96 GB | rebuild every UniProt release | none |
| ColabFold, public server | none | none | free academic service |
| ColabFold, self-hosted | ~100 GB and up | database updates | none |

A ColabFold search takes 10 to 70 seconds per protein depending on queue depth, against hours for a remote `nr` query. Rebuilding all features and retraining took 39 minutes for set08 and 10 minutes for set07.

The public server processes a few thousand alignments a day and asks for serial queries from a single IP. A 40-protein set is three orders of magnitude inside that, and the webserver has never carried enough traffic to come near it, so the public server is a legitimate default for both. `THOIPA_COLABFOLD_HOST` points the same protocol at a self-hosted MMseqs2 server if that changes, but note that self-hosting trades the disk cost back rather than removing it.

The archive returned by the server includes `msa.sh`, the exact MMseqs2 command line and database versions used. The BLAST path had no equivalent, so an alignment's provenance was previously unrecoverable.

## Conclusions

1. ColabFold and UniRef90 are statistically indistinguishable on accuracy (-0.002 ROC, p = 0.80). The 96 GB buys nothing.
2. Neither beats the 2020 `nr` archive, and UniRef90 is slightly worse than it.
3. Greater evolutionary depth does not improve THOIPA. This was the hypothesis under test and it is not supported at two different depths.
4. For the training pipeline, ColabFold via the public server is the better default: deeper alignments, no disk, no maintenance, recorded provenance, and no NCBI queue.
5. If accuracy is the goal, the evidence points at the coevolution features rather than at the homologue search.

## Limits of this analysis

34 clusters gives power to detect roughly a 0.03 to 0.05 AUC effect; a smaller real improvement would be invisible here. Proteins scoring poorly on `nr` tend to improve on any other source (rho -0.30 between baseline score and change, p = 0.06), which is regression to the mean and inflates some per-protein gains. Feature selection and hyperparameters were retuned per source, which is the realistic deployment comparison rather than an isolation of the alignment alone. The comparison against `nr` is unavoidably confounded with six years of sequence growth, and because NCBI does not archive old releases, the unconfounded version is not available to anyone.
