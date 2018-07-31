from setuptools import setup, find_packages
from os import path

# grab the long_description from the readme file
# I'm not using unicode due to a possible distutils incompatibility. Readme cannot contain non-ASCII text.
here = path.abspath(path.dirname(__file__))
with open(path.join(here, "readme.rst")) as f:
    long_description = f.read()

setup(
    name="thoipapy",
    author="Bo Zeng",
    author_email="zeng@checkmytumhomepage.de",
    description="Machine-learning prediction of residues driving homotypic transmembrane interactions.",
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="https://github.com/bojigu/thoipapy",
    download_url='https://github.com/bojigu/thoipapy/archive/0.0.1.tar.gz',
    license='MIT',
    classifiers=
    [
    'Programming Language :: Python :: 3.6',
    'License :: OSI Approved :: MIT License',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=["pandas", "numpy", "scipy", "biopython", "matplotlib", "seaborn", "django", "pytoxr","korbinian","statsmodels"],
    project_urls=
    {
    'BoZeng': 'http://frishman.wzw.tum.de/index.php?id=50',
    'FrishmanLab': 'http://frishman.wzw.tum.de/index.php?id=2',
    'LangoschLab': 'http://cbp.wzw.tum.de/index.php?id=10',
    "TU_Munich": "https://www.tum.de"
    },
    keywords="bioinformatics protein transmembrane residue conservation coevolution covariance evolutionary "
             "couplings polarity hydrophobicity randomforest machinelearning interface LIPS evolution",
    packages=find_packages(),
    version = "0.0.1",
    )