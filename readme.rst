.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/master/thoipapy/docs/THOIPApy_logo.png

THOIPApy
========

The Transmembrane HOmodimer Interface Prediction Algorithm (THOIPA) is a machine learning method for the analysis of protein-protein-interactions.

THOIPA predicts TM homodimer interface residues from evolutionary sequence information alone.

THOIPA was designed to complement experimental approaches, and also energy-based modelling of TM homodimers.

See the `FAQ`__ in the THOIPA wiki for more information.

.. _FAQ: https://github.com/bojigu/thoipapy/wiki/What-is-THOIPA%3F
__ FAQ_

What does thoipapy do?
----------------------

* download protein homologues with BLAST
* extract residue properties (e.g. residue conservation and polarity)
* trains a machine learning classifier
* validates the prediction performance
* creates heatmaps of residue properties and THOIPA prediction


Installation
------------
.. code:: python

    pip install thoipapy

THOIPA has external dependencies such as FreeContact and Phobius that are only working on Linux.

Use the "Docker" implementation detailed in the `Wiki`__ to:
* run THOIPA standalone on platforms such as Windows or MacOS
* run THOIPA standalone on Linux without installing the dependencies on your system

.. _Wiki: https://github.com/bojigu/thoipapy/wiki
__ Wiki_


Dependencies
------------

We recommend the `Anaconda python distribution`__, which contains all the required python modules (numpy, scipy, pandas,biopython and matplotlib). THOIPApy is currently tested for python 3.6.

Pip should automatically install the pytoxr package of Mark Teese.

.. _AnacondaLink: https://www.continuum.io/downloads
__ AnacondaLink_

THOIPApy depends on the command-line programs phobius and freecontact.
Both of these are only available for Linux. THOIPApy itself has been tested on several different systems running Windows and Linux.

Development status
------------------

The code has been extensively updated and annotated for public release. However is released "as is" with some known issues, limitations and legacy code.
The THOIPA standalone predictor is currently available to use. The settings file and databases used for THOIPA training are not yet released.

Usage as a standalone predictor
-------------------------------

For TMD interface residue predictions of a protein of interest, we recommend running THOIPA as a standalone program via Docker, as described in the `Wiki`__ .

.. _Wiki: https://github.com/bojigu/thoipapy/wiki
__ Wiki_

THOIPA can also be installed in Linux and used as a standalone predictor:
* The operating system needs to have freecontact, phobius, and NCBI_BLAST installed.
* The biopython wrapper for NCBIblast should be installed.

.. code:: python

    from thoipapy.thoipa import get_md5_checksum, run_THOIPA_prediction
    from thoipapy.utils import make_sure_path_exists

    protein_name = "ERBB3"
    TMD_seq = "MALTVIAGLVVIFMMLGGTFL"
    full_seq = "MVQNECRPCHENCTQGCKGPELQDCLGQTLVLIGKTHLTMALTVIAGLVVIFMMLGGTFLYWRGRRIQNKRAMRRYLERGESIEPLDPSEKANKVLA"
    out_dir = "/path/to/your/desired/output/folder"
    make_sure_path_exists(out_dir)
    md5 = get_md5_checksum(TMD_seq, full_seq)
    run_THOIPA_prediction(protein_name, md5, TMD_seq, full_seq, out_dir)

**Example Output**

.. image:: https://raw.githubusercontent.com/bojigu/thoipapy/master/thoipapy/docs/standalone_heatmap_example.png

Create your own machine learning predictor
------------------------------------------

Details on how to train THOIPA on your own datasets will be released after publication.

.. code:: python

    import THOIPApy
    settings = r"D:\data\THOIPApy_settings.xlsx"
    THOIPApy.run(settings)

License
-------

THOIPApy is free software distributed under the permissive MIT License.


Contribute
-------------

THOIPApy is not yet officially published. However, feedback regarding the installation and usage of the standalone version is appreciated. Simply email us directly, or initiate an issue in Github.


Contact
-------

For contact details, see the relevant TU-Munich websites:

Author: `Bo Zeng`__  of the `Frishman lab`__, TU-Munich, Weihenstephan Campus

Further coding and supervision: `Mark Teese`__ of the `Langosch lab`__, TU-Munich, Weihenstephan Campus

.. _BoWebsite: http://frishman.wzw.tum.de/index.php?id=50
.. _FrishmanWebsite: http://frishman.wzw.tum.de/index.php?id=2
.. _MarkWebsite: http://cbp.wzw.tum.de/index.php?id=49&L=1
.. _LangoschWebsite: http://cbp.wzw.tum.de/index.php?id=10
__ BoWebsite_
__ FrishmanWebsite_
__ MarkWebsite_
__ LangoschWebsite_


Citation
--------

Citation to be added.
Full Credits: Bo Zeng, Yao Xiao, Dmitrij Frishman, Dieter Langosch, Mark Teese

