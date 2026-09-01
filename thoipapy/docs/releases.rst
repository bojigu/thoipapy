=================
THOIPApy releases
=================

3.0.0
-----
* **security: external programs are no longer run through a shell.** ``utils.Command`` takes a
  list of arguments and runs it directly. Every external tool (rate4site, cd-hit, freecontact) was
  invoked as a single string with ``shell=True``, and rate4site's and cd-hit's filenames are built
  from the protein name, which for a standalone run or a webserver submission is the user's FASTA
  header. Every character of it reached ``/bin/sh``. An ordinary UniProt header contains pipes, so
  it was split into shell commands and the prediction failed; anything else could be injected the
  same way. The protein name is now also reduced to ``[A-Za-z0-9._-]`` and capped at 60 characters
  before it becomes part of a filename.
* **bugfix: rate4site scored the wrong sequence whenever the protein name contained a hyphen, and
  the conservation feature was then attributed to the wrong residues. Predictions change for any
  such protein and stored output should be regenerated.** The alignment handed to cd-hit was built
  by stripping gap characters from every line of the FASTA including the deflines, so the query --
  the only record whose defline is not a bare integer -- came back under a different name, was
  dropped from the alignment, and rate4site scored whichever homologue was left at the top.
  rate4site is now given the query as record zero and told to use it by name (``-a``), the residue
  identities come from the query rather than from rate4site's echo of its reference, and a
  disagreement raises instead of being merged away. **The trained model is unchanged**; this is a
  fix to one input feature, not a new predictor.
* **bugfix: a protein whose alignment was still too large at 40% identity aborted the prediction.**
  The redundancy-reduction loop walked the cd-hit threshold down to 0.20, and cd-hit rejects
  anything below 0.40 with a fatal error. It now stops at 0.40 and truncates the alignment, which
  is what the unreachable branch below it always intended. Deep alignments became normal when
  local UniRef90 searches replaced NCBI ``nr``. The threshold is stepped in whole percentage
  points rather than by repeated subtraction of 0.01 from 1.0, which drifted far enough to lose a
  round and to pick a cd-hit word size one step too small at each of its boundaries.
* **bugfix: thoipapy silenced the logging of any application that embedded it.**
  ``setup_error_logging`` called ``dictConfig`` without ``disable_existing_loggers``, which
  defaults to ``True``, so every logger the caller had already created was disabled for the rest
  of the run. A failed prediction could leave no trace anywhere, including in the handler written
  specifically to record it.
* ``THOIPA_out.csv`` and ``THOIPA_out.xlsx`` gain a ``residue number in full sequence`` column.
  The value was already computed and written to ``THOIPA_full_out.csv``, and dropped from both
  files a user actually receives, so residue numbers ran 1..N within the TMD and had to be
  converted by hand for any use against the real protein.
* ``THOIPA_out.csv`` is a genuine comma-separated CSV. It was named ``.csv``, written with tabs,
  and padded with spaces to line up in a terminal, so Excel and mail-client previews showed the
  whole prediction as a single column. Code that parsed it with ``sep="\t"`` or ``sep=r"\s+"``
  must be updated.
* ``run_THOIPA_prediction`` raises ``ValueError`` when the TMD sequence is absent from the full
  sequence, or occurs in it more than once. It used to log a warning and return normally, having
  written no output, which callers could not distinguish from success.
* ``utils.slugify`` is removed; use ``utils.safe_filename_component``. slugify lowercased and
  dropped dots and underscores, so ``P21860_ERBB3`` became ``p21860_erbb3``. The CLI's
  ``input.csv`` records the cleaned name under ``safe_name`` rather than ``slugified_name``.
* the TMD sequence is no longer used as a regular expression when locating it in the full
  sequence, so a sequence containing a regex metacharacter reports an input error rather than
  raising ``re.error`` from inside the pipeline.
* rate4site runs in its own output directory. It writes ``r4s.res``, ``r4sOrig.res`` and
  ``TheTree.txt`` into the working directory whatever ``-o`` says, so two predictions running at
  once overwrote each other's temporary files.
* the error raised when the combined features do not describe the submitted TMD now names which
  feature file disagrees and at which position, instead of printing eight sequences to compare by
  eye.
* bugfix: ``SurroundingSequence`` returned float offsets when the training pipeline built it from
  a dataframe column, so any code slicing a sequence with one raised ``TypeError``.
* bugfix: running the test suite deleted committed test fixtures. A cleanup step globbed every
  ``.xml`` and ``.txt`` under ``test/test_inputs`` and unlinked them, and because it ran as a test
  the deletion was reported as a pass.

**Exceptions that callers may be catching.** ``combine_all_features`` raises ``ValueError`` where
it raised ``IndexError``; nothing is being indexed, the features disagree about the sequence they
describe. ``rate4site_calculation`` raises ``ValueError`` in four new situations, all of which
previously produced a wrong answer or a misleading error two stages later: the first record of the
homologue alignment is not the query, that record contains a gap, the TMD is not at the expected
offset within it, or rate4site produced no score for some TMD residue. A training-set alignment
built under a different accession from the one passed to ``rate4site_calculation`` now fails here
rather than silently scoring the wrong sequence.

**Re-running over an existing output tree.** Because the accession is sanitised before it becomes
part of a filename, a protein whose name contains a character outside ``[A-Za-z0-9._-]`` gets
different intermediate filenames under ``datafiles/`` than it did in 2.x. The "skip if the output
already exists" checks will not find the old files, so the intermediates are recomputed under the
new names beside the old ones. Delete the old output directory rather than reusing it.

**Other API changes.** ``utils.Command`` gains ``stdin_bytes``, ``stdout_path`` and ``cwd``
parameters, and a ``command_string`` property for log messages; on timeout it now kills the child
rather than sending SIGTERM and then waiting for it indefinitely. The length of rate4site's output
banner is measured rather than assumed to be 13 lines. ``get_word_size`` raises instead of
returning the string ``"error"``, which would previously have been passed to cd-hit as an
argument.

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
