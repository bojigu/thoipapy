==========================================================
Protein language model embeddings, and how old nr data is
==========================================================

:Date: September 2026
:Training set: set08, 850 residues, 40 proteins, 34 CD-HIT clusters
:Blind test set: set07, 194 residues, 10 proteins
:Model: ExtraTreesClassifier, 1000 trees, results averaged over 5 seeds

Two questions were asked of the same dataset. Neither produced the answer that was expected, and
both are worth recording because the negative result is the useful part.

1. Do ESM-C protein language model embeddings improve THOIPA's interface prediction? **No.**
2. What does THOIPA lose by searching UniRef90 instead of NCBI nr? **About 0.015 AUC, and the
   alignments are deeper rather than shallower.**

Reproducing this needs ``scripts/compare_plm_embedding_configurations.py`` on the
``feature/plm-embeddings`` branch (question 1) and ``scripts/compare_blast_database_depth.py``
(question 2). The embedding work was deliberately not merged; the branch is kept as an archive.


Question 1: embeddings add nothing
===================================

Method
------

Per-residue embeddings were taken from **ESM-C** (``EvolutionaryScale/esmc-600m-2024-12``, 1152
dimensions per residue), a published model. The **full protein sequence** was submitted rather
than the transmembrane domain alone, so that each residue's representation is conditioned on the
rest of the chain, and the transmembrane rows were sliced out afterwards. In this study the model
was reached through an internal HTTP service, but nothing about the result depends on that: the
same vectors come from loading the public model directly.

1152 extra columns against 850 training residues is a 1.4:1 feature-to-row ratio, so the block was
reduced before the classifier saw it. An elastic-net-penalised logistic regression ranked the
dimensions by absolute coefficient, with the penalty strength chosen by a ``GroupKFold``
cross-validation grouped by protein, and only the top 10 or 20 were passed on.

**The ranking was refitted inside every leave-one-out fold.** Selecting 20 dimensions from 1152
once over the whole training set and then validating leave-one-out would let each fold's held-out
protein help choose the features it is then scored on. With 1152 candidates and 850 rows that bias
is not small. Each fold excludes the held-out protein and its entire CD-HIT cluster.

Results
-------

Mean per-protein ROC and PR AUC. The difference, interval and p-value are computed over the 34
CD-HIT clusters rather than the 40 proteins, because five receptor tyrosine kinases share a cluster
and are not five independent observations.

===============================================  =====  =====  ======  ================  =====  ====
configuration                                    ROC    PR     dROC    95% CI            p      Holm
===============================================  =====  =====  ======  ================  =====  ====
baseline, 26 alignment features                  0.641  0.534
+ top 20 embedding dimensions                    0.631  0.509  -0.015  [-0.041, +0.011]  0.369  1.000
+ top 10 embedding dimensions                    0.627  0.510  -0.020  [-0.040, -0.002]  0.085  0.376
embeddings alone, nothing else                   0.564  0.475  -0.086  [-0.167, +0.005]  0.030  0.180
embeddings + non-alignment features              0.613  0.495  -0.034  [-0.085, +0.021]  0.075  0.376
+ per-protein mean-centred embeddings            0.638  0.524  -0.006  [-0.030, +0.017]  0.510  1.000
+ amino-acid-mean embeddings, context removed    0.641  0.522  -0.005  [-0.033, +0.026]  0.469  1.000
===============================================  =====  =====  ======  ================  =====  ====

Every arm lands at or below the baseline, and nothing survives the Holm correction. Non-inferiority
against a 0.03 margin is not established for either alignment-free arm: the best of them has a
confidence interval reaching down to -0.085.

For scale, the mean AUC moves by only **0.003** across five forest seeds at 1000 trees, so these
differences are not reseeds. At the tuned 50 trees the seed range is 0.024, which would have
swamped every effect in the table; that is why the comparison was run at 1000.

The decisive row
----------------

The last row replaces each residue's contextual embedding with the **mean embedding for its
amino-acid type** and scores 0.641, matching the baseline and beating the real embeddings.

This is not an extraction error. On the GxxxG motif of ``1orqC4``, the two motif glycines have
cosine similarity 0.678 to each other and 0.490 to the domain mean, and all four glycines in that
helix are distinct vectors. The context is measurably present in the data. The predictor simply
does no better with it than with a flat substitution table.

Why the interface is not predictable from a language model
-----------------------------------------------------------

* **The embeddings encode the alignment, and THOIPA already has the alignment.** Masked language
  modelling learns the evolutionary distribution at each position, which is conservation and the
  PSSM implicitly. The embeddings are not a new view of the problem; they are a lossier
  re-derivation of one the model already has explicitly. That is consistent with the
  embeddings-alone arm scoring well *below* the alignment features rather than near them.

* **The interface is a quaternary property and the model sees one chain.** Which helix face
  self-associates is a fact about a complex. THOIPA's most informative features are coevolution
  scores, which work because they capture pairwise couplings across an alignment; a per-residue
  embedding vector exposes no pairwise information at all.

* **The scoring metric discards most of the embedding.** AUC ranks the ~25 residues within one
  transmembrane domain, so any component of a representation that is constant across that domain
  contributes nothing, and much of a language model embedding's variance is protein-level. The
  per-protein mean-centred arm changes almost nothing, which says the same thing from the other
  direction.

What this does not rule out
----------------------------

850 residues across 34 independent clusters is a small dataset and the intervals are wide: the
first arm's spans -0.041 to +0.011. This rules out a large gain, not a small one, and it says
nothing about fine-tuning rather than frozen embeddings, about larger models, or about attention
maps and pair representations, which is where a language model actually stores the pairwise
information this task needs.


Question 2: what UniRef90 costs
================================

The published model rests on NCBI nr. The standalone predictor and the webserver search a local
UniRef90 database instead, because a remote nr query has been measured queueing for hours and a
webserver job cannot wait. The cost of that substitution had never been measured.

The entire feature pipeline was rebuilt from UniRef90 -- BLAST, alignments, PSSM, entropy,
rate4site, coevolution, LIPS, motifs, its own feature selection and its own tuned hyperparameters
-- and scored by identical code on identical folds.

The premise was wrong: UniRef90 is deeper
------------------------------------------

UniRef90 returns **142% of nr's homologues** on set08 and is deeper for **30 of the 40** training
proteins, despite clustering UniProt at 90% identity.

The reason is a confound that cannot be removed. The nr archive was downloaded in **May 2020**; a
local UniRef90 search runs against **today's** release. Six years of sequence growth outweighs the
clustering. This is a fair description of the choice the server actually faces, but it is not a
comparison of two databases at one moment in time.

Accuracy
--------

=========================================  =======  =============  ======  ================  =====
evaluation                                 2020 nr  2026 UniRef90  change  95% CI            p
=========================================  =======  =============  ======  ================  =====
set08 leave-one-out, 40 proteins           0.641    0.630          -0.015  [-0.034, +0.002]  0.235
set07 blind test, 10 proteins              0.654    0.639          -0.015  descriptive
=========================================  =======  =============  ======  ================  =====

Two independent evaluations agree to three decimal places. Neither is significant on its own, but
the direction and magnitude replicate and 0.015 is five times the forest's own seed noise. Treat
it as a real but small penalty.

More homologues yet slightly worse accuracy means depth is not the binding constraint.

Which features degrade
-----------------------

Spearman correlation between the two feature tables over the 850 shared residues:

===============================  ================  ===========
feature                          kind              rho
===============================  ================  ===========
DImax                            coevolution       0.452
DI5mean                          coevolution       0.465
DI3mean                          coevolution       0.585
MIall_mean                       coevolution       0.756
E, H, N, G                       PSSM frequency    0.727-0.872
GxxxG, SmxxxSm, branched, mass   sequence only     1.000
===============================  ================  ===========

**The coevolution features are the least reproducible across alignment sources, and they are among
the features THOIPA relies on most.** Changing the database perturbs precisely the signal the
predictor leans on hardest. The sequence-only features are identical, as they must be, because they
never touch the alignment.


Notes on method
================

* Leave-one-out folds exclude the held-out protein and its whole CD-HIT cluster.
* Statistics are computed over 34 clusters rather than 40 proteins.
* Hyperparameters are frozen across all arms, so the pre-existing tuning optimism in
  ``04_tuned_ensemble_parameters.csv`` is a constant of every comparison rather than something that
  varies between arms. It is not removed: ``tune.py`` fits an ungrouped ``KFold`` over residues, so
  residues of one domain sit on both sides of every tuning split. That predates this work.
* PR AUC in these scripts is ``average_precision_score``. Note that
  ``validation/leave_one_out.py`` still reports it as ``auc(recall, precision)``, a trapezoidal
  interpolation that is biased at 11 to 29 points per protein. The two are not comparable.
* The UniRef90 rebuild lives outside the repository and is deliberately not DVC-tracked. The
  published nr artefacts remain the reference for the paper's numbers.
