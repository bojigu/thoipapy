import os
import platform
import re

from thoipapy import utils as utils


def return_num_tmd(s, acc,  full_seq, full_seq_fasta_file,phobius_outfile,logging):
    """Calculate the number of TMDs for the protein using the phobius prediction algorithm.

    Important when mixing crystal dataset (multipass) with single-pass protein datasets.
    Gives the extra feature n_TMDs (number of TMDs) in the protein.
    Helps the machine learning technique.

    Requires Phobius (http://phobius.sbc.su.se/data.html). Follow instructions in readme,
    including link creation so "phobius" is recognised in the console.

    Parameters
    ----------
    s : dict
        Settings dictionary
    acc : str
        Protein accession number. (e.g. UniProt acc, PDB_acc + chain letter)
    database : str
        Database name, e.g. "crystal", "NMR" or "ETRA".
    full_seq : str
        Full protein sequence
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    utils.make_sure_path_exists(full_seq_fasta_file, isfile=True)
    with open(full_seq_fasta_file, 'w') as f:
        f.write(">{}\n{}".format(acc, full_seq))
    f.close()
    if "Windows" in platform.system():
        logging.warning("phobius currently not run for Windows! Skipping phobius prediction.")
    else:
        #perl_dir = s["perl_dir"]
        #phobius_dir = s["phobius_dir"]
        #exect_str = "{} {} {}> {}".format(perl_dir, phobius_dir, full_seq_fasta_file, phobius_outfile)
        # use sudo ln -s /path/to/phobius.pl /usr/local/bin/phobius to create a link,
        # so the perl and phobius directory are not necessary
        exect_str = "phobius {}> {}".format(full_seq_fasta_file, phobius_outfile)
        command = utils.Command(exect_str)
        command.run(timeout=400, log_stderr=False)

    if os.path.exists(phobius_outfile):
        tm_num = 1
        with open(phobius_outfile) as file:
            for line in file:
                if re.search('TRANSMEM', line):
                    tm_num = tm_num + 1

        # Options are only 0 (no TMD predicted by phobius), 1, 2, 3, or 4 (4 or more TMDs predicted by phobius)
        if tm_num > 4:
            tm_num = 4
        return tm_num
    else:
        #sys.stdout.write("no phobius output file found, try to check the reason")
        #return None
        raise FileNotFoundError("{} Phobius output not found ({})".format(acc, phobius_outfile))