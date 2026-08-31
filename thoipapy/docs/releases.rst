=================
THOIPApy releases
=================

2.1.0
-----
* the example heatmap in the readme is served from the ``develop`` branch rather than from
  ``master``. The two branches hold the same commit, and ``master`` is being retired, which would
  have left the image unresolvable. No code changed.

2.0.0
-----
The first release since 1.2.0. It was developed as two internal milestones, neither of which was
published, so their notes are combined here.

* **the shipped model is retrained on 26 features, and the ``n_TMDs`` feature is removed.
  Predictions differ from 1.2.0 and stored output should be regenerated.** The 2020 scikit-learn
  0.23 pickle cannot be loaded by any current version, so retraining was unavoidable. Computing
  ``n_TMDs`` required Phobius, which is licensed academic software that was never actually
  installed, and the feature was silently constant on every prediction ever made. Removing it costs
  no accuracy. See ``docs/n_TMDs_dropped.md``.
* BLAST can search a local database instead of querying NCBI. Set ``THOIPA_LOCAL_BLAST_DB``; unset,
  the NCBI path is unchanged. ``THOIPA_NCBI_CONTACT_EMAIL`` is now required for NCBI queries, as
  their usage policy demands. ``scripts/build_local_blast_db.py`` builds the database.
* every dependency is pinned to an exact version, in ``pyproject.toml``, ``requirements.txt`` and
  ``environment.yml``. A version range resolves differently on every install, which for scientific
  software means a silently different result rather than an error.
* the training data, models and test data are tracked with DVC and served from a public store, so
  ``dvc pull`` reproduces them without credentials.
* the protein sets, ``protein_names`` and the ETRA scanning-mutagenesis data are CSV rather than
  Excel. They are single flat tables, so the workbook format made them undiffable for no gain.
* ``utils.Command`` records the subprocess return code and whether it timed out, and callers now
  check it. Previously an external tool that was missing or that failed was indistinguishable from
  one that succeeded.
* python 3.12-3.13, numpy 2, pandas 3, scikit-learn 1.9. Drops python 3.8 and Django.
* bugfix: feature selection was not reproducible. Column order came from a python set, so it
  varied with the interpreter hash seed and the selected feature set differed between runs.
  Model training and hyperparameter tuning are now seeded.
* bugfix: the standalone predictor numbered residues from 0 while the training pipeline numbered
  from 1. ``res_num_full_seq`` in ``THOIPA_full_out.csv`` was one lower than the true UniProt
  residue number. Now 1-based throughout. Stored output from earlier versions will differ.
* bugfix: CD-HIT redundancy reduction never removed a protein.
* bugfix: a stage set to FALSE in the settings spreadsheet could still run, because a non-empty
  string is truthy in python.
* bugfix: five sites hardcoded ``surr20.gaps5`` while the writers built the suffix from the
  settings, so changing either setting made the pipeline read and write different filenames.
* 10-fold cross-validation now groups by protein rather than splitting residues.
* pipeline functions take explicit parameters instead of a settings dictionary; paths come from
  ``ArtefactPaths`` and stage flags from ``RunSettings``.
* removes ~3,300 lines of unreachable code. ``import thoipapy`` now loads 21 of the package's
  modules rather than 70, and no longer pulls in scikit-learn or statsmodels.
* publication figure code moved to ``thoipapy/paper_figures``, unmaintained.
* packaging moved to pyproject.toml; fixes a bug where the settings file the predictor loads on
  every run was missing from the built wheel.
* settings moved from a multi-sheet Excel workbook to CSV, so they can be diffed in git and
  edited without Excel. Values are typed on load, so a setting of 0 can no longer be read as a
  non-empty string and behave as if it were 1. ~20 settings that no code read were removed.
* see the readme for how to reproduce the published 2020 results.

1.2.0
-----
* bugfix: prevent error when (highly random) sequence in rate4site output did not match initial TMD
* feature: add test for feature extraction and machine-learning model creation for a protein list
* feature: rename several inputs in settings file, including base_dir & data_dir. This affects the creation of new ML models, but not standalone predictions.

1.1.3
-----
* fix bug where cd-hit not found when installed using apt-get
* improve functional tests by using pre-downloaded homologue

1.1.2
-----
* fix requirements.txt

1.1.1
-----
* fix bug where CREATE_VIEW indicates malformed xml from ncbi blast servers

1.1.0
-----
* update settings files

1.0.1
-----
* fix missing psutil dependency
* fix biopython syntax to remove warning message

1.0.0
-----
* added new prediction features including conservation calculated with rate4site
* removed blind test TMDs from the training dataset
* added a feature selection pipeline, including the removal of duplicate features, and selection of best predictive features
* added automatic tuning of machine-learning predictor
* excluded putative distant homologues from training in each iteration of leave-one-out validation
* added scripts for conducting bootstrapped t-test, comparing interface and non-interface residues
* extended the output for a single file/experiment to include the mean EC50 values for replicates with identical sample names (issue #8)

0.0.7
-----
* initial public release
