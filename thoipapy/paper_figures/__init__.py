"""Figures for the 2020 CSBJ publication. NOT MAINTAINED.

This code was written to produce the figures in

    Xiao, Zeng, Berner, Frishman, Langosch and Teese (2020),
    Computational and Structural Biotechnology Journal.

It was last exercised in 2020 and is not covered by any test. Every pipeline stage that calls it
is switched off in both shipped settings files, and the functional test disables the validation
stage explicitly ("currently failing due to divide by zero error").

It is kept because regenerating a published figure is occasionally useful, and deleted code is
harder to recover than dormant code. It is deliberately excluded from the refactoring work applied
to the rest of the package: do not assume anything here reflects current conventions, and if you
need a figure regenerated, expect to fix that module first.

The maintained surface is thoipapy/features, homologues, feature_importance, ML_model, clustering,
common, utils, predict and run.
"""
