.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/develop/thoipapy/docs/THOIPA_banner.png

THOIPApy
========

The Transmembrane HOmodimer Interface Prediction Algorithm (THOIPA) is a machine learning method for the analysis of protein-protein-interactions.

THOIPA predicts transmembrane homodimer interface residues from evolutionary sequence information.

THOIPA helps predict potential homotypic transmembrane interface residues, which can then be verified experimentally.
THOIPA also aids in the energy-based modelling of transmembrane homodimers.

Important links:

* `THOIPA webserver <http://www.thoipa.org>`_
* `THOIPA FAQ <https://github.com/bojigu/thoipapy/wiki/What-is-THOIPA%3F>`_
* `THOIPA wiki main page <https://github.com/bojigu/thoipapy/wiki/THOIPA-wiki-main-page>`_


How does thoipapy work?
-----------------------

* downloads protein homologues with BLAST
* extracts residue properties (e.g. residue conservation and polarity)
* trains a machine learning classifier
* validates the prediction performance
* creates heatmaps of residue properties and THOIPA prediction


Installation
------------

THOIPA is on PyPI, but read this before installing it.

.. code:: bash

    pip install "thoipapy>=2.1"

**Ask for 2.1 or newer.** Releases up to 1.2.0 were published between 2018 and 2021 and no longer
work: the model they ship is a scikit-learn 0.23 pickle, which no current scikit-learn can load,
so they install without error and then fail on the first prediction.

**Install it into a new, empty environment.** Every dependency is pinned to an exact version
(``numpy==2.5.2``, ``pandas==3.0.5``, and so on), so pip will almost certainly refuse to install
it alongside packages you already have. That is intended, not a defect: THOIPA is not actively
maintained, and exact pins are what stop an unattended install from silently resolving to a
different numerical result years from now.

The same pinning means that on a future python, where those exact versions have no wheel, the
install will fail outright rather than succeed and misbehave. If that happens, use the conda
environment below, which also pins the interpreter.

.. code:: bash

    conda env create -f environment.yml
    conda activate thoipapy
    pip install -e .

This is the supported path, and the only one that pins python itself along with the external
command-line tools.

THOIPA has only been tested on Linux, because of its reliance on external programs such as
FreeContact, CD-HIT and rate4site.


Dependencies
------------

THOIPApy requires python 3.12 or newer and is tested on python 3.13.

**Every dependency is pinned to an exact version.** This is deliberate. THOIPA is scientific
software whose output has to stay reproducible, and a version range resolves to something
different every time it is installed. The failure that causes is not a helpful error at install
time; it is a prediction that quietly differs from the published one. Pinning trades the ability
to install alongside arbitrary other packages for the ability to still work in five years. Use a
dedicated environment.

Pins live in three places, which are kept in step:

=========================  ==================================================================
file                       pins
=========================  ==================================================================
``pyproject.toml``         python packages; what ``pip install .`` reads
``requirements.txt``       the same packages plus their transitive dependencies
``environment.yml``        the python interpreter and the external command-line tools
=========================  ==================================================================

The external tools are pinned for the same reason as the python packages: a different
``rate4site`` or ``cd-hit`` changes the conservation scores and the redundancy reduction, and
therefore the features and the prediction.

You are of course free to relax the pins and test a newer combination. Nothing guarantees that
it resolves, or that predictions still match. If you do, run the full test suite, including the
golden-file regression tests, before trusting the result.

Several pipeline steps call external command-line programs rather than python libraries.
``blast``, ``cd-hit`` and ``rate4site`` are installed by ``environment.yml``. ``freecontact`` is
not packaged for conda on any channel, so it is installed separately, either with root:

.. code:: bash

    sudo apt-get install freecontact

or without, into the active conda environment:

.. code:: bash

    ./scripts/install_freecontact.sh

Phobius is no longer required. It was needed only for the ``n_TMDs`` feature, which was removed
in 2.1.0. See ``docs/n_TMDs_dropped.md``.

Getting the data (training sets, models, test data)
---------------------------------------------------

**If you only want to run predictions, you do not need any of this.** The trained model ships
inside the package, so ``pip install thoipapy`` gives you a working predictor with no data
download, no clone, and no DVC. Everything below is for developers who want to retrain the model,
re-run the validation, or run the test suite.

The code lives in git; the data does not. ``data/`` is about 128 MB of training sets, homologue
alignments, extracted features and prediction outputs, which is too large and too binary to keep
in a git repository. It is stored separately and fetched with a tool called **DVC** (Data Version
Control).

**If you have not used DVC before**, the idea is simple. DVC leaves a small text file in git for
each data directory, for example ``data/features.dvc``. That file records a checksum, not the
data. When you run ``dvc pull``, DVC reads those checksums and downloads the matching files from
a data store into ``data/``. It is, in effect, ``git clone`` for large files. You do not need an
account, a login, or any credentials to fetch this project's data.

Full documentation: https://dvc.org/doc . A good short introduction is
https://dvc.org/doc/start/data-management/data-versioning .

Step 1: install DVC
~~~~~~~~~~~~~~~~~~~

If you created the conda environment from ``environment.yml`` you already have it, and can skip to
step 2. Otherwise, pick whichever you prefer::

    # with conda or mamba (recommended, matches environment.yml)
    conda install -c conda-forge dvc dvc-http

    # or with pip
    pip install "dvc[http]"

``dvc-http`` matters: this project's public data store is served over HTTPS, and plain ``dvc``
cannot read an HTTPS remote without it. If you see ``URL 'https://...' is supported but requires
these missing dependencies``, this is what is missing.

Check it worked::

    dvc --version

Step 2: fetch the data
~~~~~~~~~~~~~~~~~~~~~~

From the root of your clone::

    dvc pull -j 4

Expect roughly 3800 files, 128 MB, and a minute or two on a normal connection.

**Use** ``-j 4``. The ``-j`` flag sets how many files DVC downloads at once. Its default is four
times your CPU count, which is enough parallel requests that the data store starts refusing them,
and the pull fails with ``429 Too Many Requests``. ``-j 4`` avoids this. If you still see 429
errors, lower it further with ``-j 2``.

To confirm everything arrived::

    dvc status

``Data and pipelines are up to date.`` means the contents of ``data/`` match the checksums
recorded in git, so you have exactly the files the authors had.

Fetching only part of the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The full set is rarely needed. To fetch one directory, name its ``.dvc`` file::

    dvc pull -j 4 test/regression_data.dvc  # run the test suite                    28 KB
    dvc pull -j 4 data/features.dvc         # retrain the model                      18 MB
    dvc pull -j 4 data/results.dvc          # trained models and validation output   25 MB
    dvc pull -j 4 data/homologues.dvc       # re-extract features from alignments    60 MB

The remaining directories are small: ``data/Proteins`` (420 KB), ``data/Input_data`` (272 KB) and
``data/ETRA_data`` (200 KB). ``data/Predictions`` (26 MB) holds the published predictions and is
needed only to compare against them.

Common problems
~~~~~~~~~~~~~~~

``429 Too Many Requests``
    Too many parallel downloads. Use ``dvc pull -j 4``, or ``-j 2`` on a fast connection.

``URL 'https://...' is supported but requires these missing dependencies``
    ``dvc-http`` is not installed. See step 1.

``No DVC remote is specified``
    You are not in the repository root, or ``.dvc/config`` is missing. Run ``dvc remote list``;
    it should print a remote named ``public``.

``dvc: command not found``
    DVC installed into a different environment than the one you have active. Run
    ``conda activate thoipapy`` first.

The pull is slow or stalls
    The data store is rate-limited. This is expected and ``-j 4`` handles it. DVC resumes where it
    stopped, so re-running ``dvc pull -j 4`` after an interruption is safe and will not re-download
    what you already have.

Where the data comes from
~~~~~~~~~~~~~~~~~~~~~~~~~

``dvc pull`` fetches from a publicly readable Cloudflare R2 bucket, declared in ``.dvc/config``.
It is read-only and anonymous by design: anyone can fetch the data, nobody can overwrite it.
Maintainers push with a separate credentialed remote configured in ``.dvc/config.local``, which is
not tracked by git.

The data originates from the dataset published with the paper, archived at the Open Science
Framework: https://osf.io/txjev/. The OSF archive is the citable record of what the paper used and
does not change. The DVC store is the working copy the pipeline reads, and it tracks the code.
Where the two differ, OSF is authoritative for the published results.

The protein sets
~~~~~~~~~~~~~~~~

``data/sets/`` holds the three protein sets published with the paper, matching
``data/protein_lists/`` in the `OSF repository <https://osf.io/txjev/>`_:

========================================  ===  =========================================
set                                       n    contents
========================================  ===  =========================================
``set05_ETRA_NMR_crystal_nr.csv``         50   all proteins, redundancy-reduced
``set07_test.csv``                        10   blind test set
``set08_train.csv``                       40   training set
========================================  ===  =========================================

set05 is the union of set07 and set08. The copies here carry a few extra annotation columns
beyond the OSF versions; the shared columns are identical.

These were Excel workbooks and are now CSV. They are single flat tables, so the workbook format
bought nothing and made them undiffable in git and awkward to read without Excel. The same applies
to ``data/protein_names.csv`` and the ETRA scanning-mutagenesis data in
``data/ETRA_data/Average_with_interface/``. Pipeline outputs under ``data/results/`` remain Excel,
because those genuinely hold several tables per file.


Running BLAST locally instead of at NCBI
-----------------------------------------

**Recommended: run THOIPA against a local UniRef90 database rather than querying NCBI.**

Why this changed
~~~~~~~~~~~~~~~~

THOIPA was built in 2017-2020, when a BLAST query to NCBI returned in a few minutes. The sequence
databases have grown enormously since: ``nr`` now holds over a billion sequences, and NCBI queues
automated clients accordingly. A probe of a *three-residue* query in August 2026 returned
``RTOE=11165`` -- NCBI's own estimate of **3.1 hours** before results would be ready. A prediction
that took minutes when the paper was published can now take hours, or be refused outright if the
same client queries repeatedly.

That is not a fault in THOIPA and it is not something THOIPA can fix from the client side. The
only reliable answer is to stop depending on a shared, queued, remote service.

What you gain and what you lose
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**This is a breaking change and predictions will not match the published results.** A different
database yields a different homologue set, a different alignment, and therefore different
conservation, coevolution and polarity features. The scores move.

In exchange, a full prediction goes **from hours to a few minutes** -- seconds, on a small
database -- and it runs offline and deterministically, with no queue and no third-party service in
the path.

Setting up UniRef90
~~~~~~~~~~~~~~~~~~~

UniRef90 clusters UniProt at 90% identity: 121 million clusters against ``nr``'s billion-plus
sequences, with little practical loss for this pipeline, which collapses exact-duplicate TMD
sequences and runs CD-HIT anyway.

Requirements: about **90 GB of free disk** (30 GB download, 60-80 GB built) and BLAST+.

::

    mamba install -c bioconda blast

    # Downloads UniRef90 and runs makeblastdb. Idempotent: an existing download is not refetched.
    # Expect several hours for the download and 1-2 hours to build.
    python scripts/build_local_blast_db.py

    export THOIPA_LOCAL_BLAST_DB=/mnt/shared/blastdb/uniref90

Or by hand::

    wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    gunzip uniref90.fasta.gz
    makeblastdb -in uniref90.fasta -dbtype prot -out uniref90 -title uniref90 -parse_seqids

With ``THOIPA_LOCAL_BLAST_DB`` set, predictions search that database. Unset, they query
NCBI as before, and ``THOIPA_NCBI_CONTACT_EMAIL`` must be set because NCBI requires
automated clients to identify themselves.

A smaller database for tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set ``DATABASE_TO_BUILD = "swissprot"`` in the script for a 575k-sequence, ~340 MB database that
builds in ten seconds and answers in under one. Use it for tests, CI and development.

**Do not use Swiss-Prot for production.** It yields roughly a quarter of ``nr``'s alignment depth:
26 hits against 76 unique TMD homologues for ``1xioA4``, 18 against 101 for ``4hksA1``, and only 10
survived filtering for glycophorin A. Conservation, rate4site and coevolution are the 3rd, 4th and
7th most important model features and all derive from that alignment.

A note on database choice
~~~~~~~~~~~~~~~~~~~~~~~~~

A database filtered to membrane proteins does **not** help. THOIPA needs homologues of the query
protein itself -- glycophorin A's useful hits are glycophorins in other species -- so what matters
is taxonomic depth within the query's own family, not enrichment for membrane proteins. Filtering
on a transmembrane annotation would also silently drop genuine homologues that merely lack the
annotation, and alignment depth is precisely what must not be lost.

NCBI does not host UniRef90, so the remote path cannot use it. NCBI's BLAST service offers ``nr``,
``refseq_protein``, ``swissprot``, ``pdb``, ``env_nr``, ``tsa_nr``, ``landmark`` and ``pataa``; of
those only ``nr`` has comparable depth, which is why the remote default is unchanged.


Reproducing the published results
---------------------------------

Versions from 2.0.0 onwards include corrections and improvements made after publication, so their
output differs slightly from the numbers in the 2020 paper. The differences are small and do not
change the conclusions, but they are real.

**To reproduce the results as published**, use the state of the repository at the time of
publication rather than the current release:

* thoipapy **v1.2.0**, git commit ``aca9692``
* the dataset published via the `Open Science Foundation <https://osf.io/txjev/>`_

Note that v1.2.0 requires python 3.8 and the dependency versions pinned in its ``requirements.txt``;
it will not run on a current scientific python stack.

**For any new work**, use the current release. The main post-publication corrections are that
feature selection is now deterministic (it previously depended on the interpreter's hash seed, so
the selected feature set varied between runs), model training and hyperparameter tuning are seeded,
and several ineffective filters and dead assertions were repaired.


Usage as a standalone predictor
-------------------------------

* for local predictions on linux, first install NCBI_BLAST, biopython, freecontact, CD-HIT and rate4site
* see ``test/functional/test_standalone_prediction.py`` for the current run syntax, typically

.. code:: python

    from thoipapy import run_THOIPA_prediction
    from thoipapy.predict import get_md5_checksum
    from thoipapy.utils import make_sure_path_exists

    protein_name = "ERBB3"
    TMD_seq = "MALTVIAGLVVIFMMLGGTFL"
    full_seq = "MVQNECRPCHENCTQGCKGPELQDCLGQTLVLIGKTHLTMALTVIAGLVVIFMMLGGTFLYWRGRRIQNKRAMRRYLERGESIEPLDPSEKANKVLA"
    out_dir = "/path/to/your/desired/output/folder"
    make_sure_path_exists(out_dir)
    md5 = get_md5_checksum(TMD_seq, full_seq)
    run_THOIPA_prediction(protein_name, md5, TMD_seq, full_seq, out_dir)


**Example Output**

* the output includes a csv showing the THOIPA prediction for each residue, as well as a heatmap figure as a summary
* below is a heatmap showing the THOIPA prediction, and underlying conservation, relative polarity, and coevolution

.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/master/thoipapy/docs/standalone_heatmap_example.png


Create your own machine learning predictor
------------------------------------------

* THOIPA can be retrained to any dataset of your choice
* the original set of training sequences and other resources are published via the `Open Science Foundation <https://osf.io/txjev/>`_, and are also what ``dvc pull`` fetches into ``data/``
* the THOIPA feature extraction, feature selection, and training pipeline is fully automated
* contact us for an introduction to the THOIPA software pipeline and settings

The data, sets and base directories are resolved from the repository itself (see
``thoipapy/paths.py``); they are no longer read from the settings spreadsheet. Protein sets are
defined in ``data/sets/``. Because the training data is not shipped inside the installed package,
the training pipeline can only be run from a checkout of this repository, not from
``pip install thoipapy``.


Where things are
----------------

============================================  ==========================================
path                                          contents
============================================  ==========================================
``thoipapy/predict.py``                       standalone prediction for a single sequence
``thoipapy/run.py``                           the training and validation pipeline
``thoipapy/features/``                        one module per feature
``thoipapy/homologues/``                      BLAST search and alignment parsing
``thoipapy/ML_model/``                        training, tuning, and the shipped model
``thoipapy/validation/``                      cross-validation and performance metrics
``thoipapy/paths.py``                         resolves ``data/``, ``data/sets/`` and settings
``thoipapy/artefacts.py``                     ``ArtefactPaths``: the path of every pipeline file
``thoipapy/run_settings.py``                  ``RunSettings``: which pipeline stages run
``thoipapy/paper_figures/``                   figure code for the 2020 paper, unmaintained
============================================  ==========================================

Settings are CSV. ``thoipapy/setting/run_settings_example.csv`` drives a training run,
``standalone_run_settings.csv`` a single prediction, and ``model_features.csv`` lists the features
eligible for the model. Each row carries a ``type`` column and is converted on load.

Pipeline functions take explicit named parameters rather than a settings dictionary, so a
signature states what the function reads::

    def lipo_from_pssm_mult_prot(paths, df_set, lipophilicity_scale, logging):

Artefact filenames encode the settings that produced them (``surr20.gaps5``). ``ArtefactPaths``
computes them, so writers and readers cannot disagree when a setting changes.

Run the training pipeline with::

    python -m thoipapy.run

The settings path is hardcoded in ``run.py``. To use a different one, import ``run()`` and pass
your own settings dict.


License
-------

THOIPApy is free software distributed under the permissive MIT License.


Contribute
-------------

* Contributors are welcome.
* For feedback or troubleshooting, please email us directly and initiate an issue in Github.


Contact
-------

* Mark Teese, `22DataCatalysis GmbH <https://www.datacatalysis.com>`_, formerly of the `Langosch Lab <http://cbp.wzw.tum.de/index.php?id=10>`_ at the `Technical University of Munich <https://www.tum.de/en/>`_
* `Bo Zeng <http://frishman.wzw.tum.de/index.php?id=50>`_, `Chinese Academy of Sciences, Beijing <http://english.cas.cn/>`_ formerly of the `Frishman Lab <http://frishman.wzw.tum.de/index.php?id=2>`_ at the `Technical University of Munich <https://www.tum.de/en/>`_

.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/develop/thoipapy/docs/signac_seine_bei_samois_mt.png
   :height: 150px
   :width: 250px

.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/develop/thoipapy/docs/signac_notredame_bz.png
   :height: 120px
   :width: 250px


Citation
--------

`Yao Xiao, Bo Zeng, Nicola Berner, Dmitrij Frishman, Dieter Langosch, and Mark George Teese (2020)
Experimental determination and data-driven prediction of homotypic transmembrane domain interfaces,
Computational and Structural Biotechnology Journal <https://doi.org/10.1016/j.csbj.2020.09.035>`_
