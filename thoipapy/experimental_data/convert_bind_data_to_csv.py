import csv
import os
import sys


def convert_bind_data_to_csv(s, df_set, logging):
    """Convert bind data (interface-residue or non-interface-residue) to a csv, for all proteins in a list.

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
    for i in df_set.index:
        acc = df_set.loc[i, "acc"]
        bind_file = os.path.join(s["data_dir"], "features", "structure", "%s.4.0closedist") % acc
        csv_output_file = os.path.join(s["data_dir"], "features", "structure", "%s.4.0closedist.csv") % acc
        if os.path.isfile(bind_file):
            try:
                with open(bind_file, "r") as bind_file_handle:
                    # csv_output_file=os.path.join(s["data_dir"], "features", "structure","NoRedundPro/%s.csv") %acc
                    with open(csv_output_file, "w") as csv_output_file_handle:
                        writer = csv.writer(csv_output_file_handle, delimiter=',', lineterminator='\n')
                        writer.writerow(["residue_num", "residue_name", "bind", "closedist"])
                        i = 1
                        for row in bind_file_handle:
                            array = row.split()
                            if len(array) == 6:
                                csv_header_for_bind_file = [i, array[3].strip('"'), array[5], array[4]]
                                writer.writerow(csv_header_for_bind_file)
                                i = i + 1
            except:
                sys.stdout.write("bind file parsing occures errors")
        else:
            sys.stdout.write("{} convert_bind_data_to_csv failed. {} not found".format(acc, bind_file))
