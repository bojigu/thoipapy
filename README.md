![THOIPApy logo](docs/logo/THOIPApy_logo.png)

# THOIPApy
THOIPApy is a program of Transmembrane HOmodimer Interface Prediction Algorithm (THOIPA) in python.

THOIPA predicts TM homodimer interface residues from evolutionary sequence information alone. THOIPA is unique in that it
was trained exclusively on TM homodimers within crystal and NMR structures, and yet has the ability to predict
the dynamic interfaces identified with ToxR-based assays in live cells.


# Features
## Robust data analysis
 * fully automated:
  * transmembrane protein homologous download with BLAST against NCBI nr.
  * evolution and co-evolutionary feature calculation.
  * interface prediction with random forest classifier
  * prediction output be fitted with sine curve
  * prediciton output analysis
 * send result through email

## Designed for humans
 * easy-to-use excel files:
  * excel settings file
  * excel input file with setting data
  * excel output file with interface prediction output data
 * simple graphical output:
  * sine curve with predicted interface residue mark

## Customisable
 - simple python syntax
 - open-source software
 - built on powerful numpy, scipy, and pandas packages

# Development status

THOIPApy has been used extensively for the analysis of ToxR datasets from the lab of Dieter Langosch at the Technical University of Munich in Germany.

The code has been extensively updated and annotated for public release.

However the module is still under development and is released "as is" with some known issues, limitations and legacy code. As a newly released module, bugfixing related to diverse operating systems, python versions, data formats, and experimental data types should be expected.

# Installation

THOIPApy requires python 3.x (currently written for 3.5). We recommend the Anaconda python distribution, which contains all the required python modules (numpy, scipy, pandas and matplotlib).
https://www.continuum.io/downloads

To install THOIPApy:
 * download and unpack the module from Github
 * open the command console. Navigate to the THOIPApy folder that contains setup.py
 * run the following command:
   `python setup.py install`

# Usage
Using THOIPApy requires only the following:

## 1. Prepare your data
 * use the excel or microplate templates in the THOIPApy/templates folder
 * for the generic excel format, simply open the template and paste in your dose and response data.

## 2. Update an excel settings file
 * copy the THOIPApy_settings_template.xlsx from THOIPApy/templates
 * open the excel file, input the name and location of your datafiles, and the desired location for your output files
 * write "TRUE" next to the files you want to examine
![01_run_curvefit_settings](docs/images/01_run_curvefit_settings.png)

## 3. tell THOIPApy to "run"
 * run the ipython/jupyter notebook, which opens a python interpreter in your web browser
 * paste in the following three lines. Replace the location of your settings file.
 * hit Ctrl-Enter to run
 * based on your output, adjust the quality thresholds in the settings file to suit your data
```
import THOIPApy
settings = r"D:\data\THOIPApy_settings.xlsx"
THOIPApy.run_curvefit(settings)

# Contribute#
If you encounter a bug or THOIPApy doesn't work for any reason, please send an email to zeng /at/ wzw.tum.de or initiate an issue in Github.

Non-programmers can contribute by:
 - testing THOIPApy with your particular datasets
 - suggesting features
 - improving the readme and documentation

Programmer contributions are very welcome:
 - adapting THOIPApy for more diverse input files and datatypes. Currently accepted are only excel.
 - adding your own desired features
 - improving code, or fixing known issues.

# License#
THOIPApy is free software distributed under the MIT License.

# Citation#
If you use THOIPApy in your research, please cite as follows:
"Determination and prediction of interface residues from transmembrane helix dimers. Yao Xiao?, Bo Zeng?, Dmitrij Frishman, Dieter Langosch, Mark George Teese*
."
