![THOIPApy logo](docs/logo/THOIPApy_logo.png)

# THOIPApy#
THOIPApy is a program of Transmembrane HOmodimer Interface Prediction Algorithm (THOIPA) in python.

THOIPA predicts TM homodimer interface residues from evolutionary sequence information alone. THOIPA is unique in that it
was trained exclusively on TM homodimers within crystal and NMR structures, and yet has the ability to predict
the dynamic interfaces identified with ToxR-based assays in live cells.


# Features#
<h3>Robust data analysis.</h3>
 - fully automated:
  - fitting of sigmoidal curves to dose-response data
  - calculation of EC50, LD50 or IC50 values.
  - high-throughput analysis
  - comparison of EC50 values from different experiments
  - calculation of EC25 or EC75 values
 - accepts REAL biological data
  - pre-filtering excludes nonsense data
  - judgefit module identifies low-quality EC50 values

<h3>Designed for humans.</h3>
 - easy-to-use excel files:
  - excel settings file
  - excel input file with dose and response data
  - excel output file with EC50 values
 - simple graphical output:
  - sigmoidal curves with EC50 shown on graph
  - daily summary barcharts and curves

<h3>Customisable.</h3>
 - simple python syntax
 - open-source software
 - built on powerful numpy, scipy, and pandas packages

# Development status#

THOIPApy has been used extensively for the analysis of LD50 assays (Z-shaped curves) by bachelor, master and PhD students within the lab of Dieter Langosch at the Technical University of Munich in Germany.

The code has been extensively updated and annotated for public release.

However the module is still under development and is released "as is" with some known issues, limitations and legacy code. As a newly released module, bugfixing related to diverse operating systems, python versions, data formats, and experimental data types should be expected.

# Installation#

THOIPApy requires python 3.x (currently written for 3.5). We recommend the Anaconda python distribution, which contains all the required python modules (numpy, scipy, pandas and matplotlib).
https://www.continuum.io/downloads

To install THOIPApy:
 - download and unpack the module from Github
 - open the command console. Navigate to the THOIPApy folder that contains setup.py
 - run the following command:
   `python setup.py install`

# Usage#
Using THOIPApy requires only the following:
1) Prepare your data, 2) update an excel settings file, 3) tell THOIPApy to "run".
<br />
<h3>1) Prepare your data.</h3>
 - use the excel or microplate templates in the THOIPApy/templates folder
 - for the generic excel format, simply open the template and paste in your dose and response data.

<h3>2) Update an excel settings file</h3>
 - copy the THOIPApy_settings_template.xlsx from THOIPApy/templates
 - open the excel file, input the name and location of your datafiles, and the desired location for your output files
 - write "TRUE" next to the files you want to examine
![01_run_curvefit_settings](docs/images/01_run_curvefit_settings.png)

<h3>3) tell THOIPApy to "run".</h3>
 - run the ipython/jupyter notebook, which opens a python interpreter in your web browser
 - paste in the following three lines. Replace the location of your settings file.
 - hit Ctrl-Enter to run
 - based on your output, adjust the quality thresholds in the settings file to suit your data
```
import THOIPApy
settings = r"D:\data\THOIPApy_settings.xlsx"
THOIPApy.run_curvefit(settings)
THOIPApy.run_gatherer(settings)
```

# Test#
 - try the example excel files in the THOIPApy/examples folder before switching to your own data.

# THOIPApy output#

<h3>run_curvefit program</h3>
 - individual dose-response curves
 - automatic judging of data quality
 - daily summary curves, barchart and more!
![curve_fit_output_sample3](docs/images/curve_fit_output_sample3.png)
![generated_data_0EC50_analysis_fig](docs/images/generated_data_0EC50_analysis_fig.png)

<h3>run_gatherer program</h3>
 - combines data from multiple experiments
 - excludes EC50 values that are not of sufficient quality, according to user-defined thresholds
 - bar charts with mean and SEM over all selected experiments
 - scatter plots showing individual datapoints for each day/experiment, and more!
![analysis_output_figures](docs/images/analysis_output_figures.png)

<h3>compare_rawdata program</h3>
 - collects raw data and dose-response curves from multiple experiments
 - compares datapoints and fitted curves between the selected samples
![20160527_0_compare_raw](docs/images/20160527_0_compare_raw.png)

# Contribute#
If you encounter a bug or THOIPApy doesn't work for any reason, please send an email to mark.teese /at/ tum.de or initiate an issue in Github.

Non-programmers can contribute by:
 - testing THOIPApy with your particular datasets
 - suggesting features
 - improving the readme and documentation

Programmer contributions are very welcome:
 - adapting THOIPApy for more diverse input files and datatypes. Currently accepted are A) excel, B) BMG FluoStar, and C) Molecular Devices SoftMax(VersaMax) files.
 - adding your own desired features
 - improving code, or fixing known issues.

# License#
THOIPApy is free software distributed under the GNU General Public License version 3.

# Citation#
Currently there in no scientific article associated with this module. If you use THOIPApy in your research, please cite as follows:
"EC50 values were calculated using the open-source THOIPApy module in python (Mark Teese, Technical University of Munich)."
