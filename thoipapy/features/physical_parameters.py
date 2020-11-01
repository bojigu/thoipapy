import csv
import os
import re
from pathlib import Path
from shutil import copyfile
import pandas as pd
import thoipapy


def add_physical_parameters_to_features_mult_prot(s, df_set, logging):
    """Run add_physical_parameters_to_features for multiple proteins.

    Parameters
    ----------
    s : dict
        Settings dictionary
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    logging.info('adding physical parameters into traindata')
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        database = df_set.loc[i, "database"]
        feature_combined_file = os.path.join(s["thoipapy_data_folder"], "features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        # feature_combined_file_incl_phys_param = os.path.join(s["thoipapy_data_folder"], "features", "combined", database,
        #                                                     "{}.surr{}.gaps{}.combined_features_incl_phys_param.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
        add_physical_parameters_to_features(acc, feature_combined_file, logging)


def add_physical_parameters_to_features(acc, feature_combined_file, logging):
    """Add physical parameters (e.g. PI of amino acid) to the features.

    Parameters
    ----------
	acc : str
        Protein accession (e.g. UniProt, PDB)
    feature_combined_file : str
        Path to csv with all features combined
    logging : logging.Logger
        Python object with settings for logging to console and file.
    """
    feature_combined_file_incl_phys_param: Path = Path(str(feature_combined_file)[:-4] + "_incl_phys_param.csv")
    thoipapy_module_path = os.path.dirname(os.path.abspath(thoipapy.__file__))
    physical_parameter_file = os.path.join(thoipapy_module_path, "setting", "Physical_Property_csv.txt")

    dictionary = {}

    with open(physical_parameter_file, "r") as physical_parameter_file_handle:

        if os.path.isfile(feature_combined_file):
            with open(feature_combined_file, "r") as train_data_file_handle:
                # train_data_add_physical_parameter_file = os.path.join("/scratch/zeng/thoipapy/features/coevolution/zfullfreecontact","%s.physipara.traindata1.csv") %acc
                train_data_add_physical_parameter_file_handle = open(feature_combined_file_incl_phys_param, "w")
                # train_data_physical_parameter_file_handle = open(r"/scratch/zeng/thoipapy/features/5hej_A2.mem.2gap.physipara.traindata.csv", "w")
                writer = csv.writer(train_data_add_physical_parameter_file_handle, delimiter=',', quotechar='"',
                                    lineterminator='\n',
                                    quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                for row in physical_parameter_file_handle:
                    if re.search("^Hydro", row):
                        continue
                    else:
                        array = row.split()
                        dictionary[array[0]] = array[1:16]
                for row1 in train_data_file_handle:
                    if re.search(r"residue_num", row1):
                        array1 = row1.rstrip().split(",")
                        array1[42:14] = ["Hydrophobicity_sAA", "Charge_sAA", "PI_sAA", "LIPSI_sAA", "LIPSM_sAA", "Hydrophobic_sAA", "Aliphatic_sAA", "Aromatic_sAA", "Polar_sAA", "Negative_sAA", "Positive_sAA", "Small_sAA", "branched",
                                         "mass", "Volume_sAA"]  # Mass_sAA
                        # array2 = array1[0:31]
                        # array2.extend(["Hydrophobicity", "Charge", "PI", "LIPS", "LIPSM", "Hydrophobic", "Aliphatic", "Aromatic", "Polar","Negative", "Positive", "Small", "Cbbranched", "Mass", "Volumn", array1[30].rstrip()])
                        writer.writerow(array1)
                    else:
                        array3 = row1.rstrip().split(",")
                        array3[42:14] = dictionary[array3[2]]
                        writer.writerow(array3)
                train_data_add_physical_parameter_file_handle.close()

        else:
            logging.warning("{} does not exist".format(feature_combined_file))

    df = pd.read_csv(feature_combined_file_incl_phys_param, index_col=0)
    cols_to_delete = ["Charge_sAA", "LIPSI_sAA", "LIPSM_sAA", "Hydrophobic_sAA", "Aliphatic_sAA", "Negative_sAA", "Positive_sAA", "Volume_sAA"]
    df.drop(cols_to_delete, axis=1, inplace=True)
    df.to_csv(feature_combined_file_incl_phys_param)

    # overwrite the original feature_combined_file, and delete feature_combined_file_incl_phys_param
    copyfile(feature_combined_file_incl_phys_param, feature_combined_file)
    try:
        os.remove(feature_combined_file_incl_phys_param)
    except:
        logging.warning("{} could not be removed".format(feature_combined_file_incl_phys_param))

    logging.info("{} add_physical_parameters_to_features_mult_prot finished. (updated {})".format(acc, feature_combined_file))
