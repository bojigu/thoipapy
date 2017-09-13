import pandas as pd
import numpy as np
def parse_Predicted_Output(thoipapyset_,out_put_file_loc,output_parse_file,logging):
    logging.info('starting to parsing prediction output file')
    orig_csv = out_put_file_loc
    test_out_csv = output_parse_file
    # read the csv. Skip the final number at the bottom (8 for GpA).
    #df = pd.read_csv(orig_csv, index_col=0, skiprows = [24])
    df = pd.read_csv(orig_csv, index_col=0)
    # set the TM RESIDUE as the index
    df.set_index("Res", inplace=True)
    # add the position in the full sequence
    # for GpA, I've just added 91. I guess you'll use REGEX to find the position in the original sequence, and give an error.
    #df["Pos"] = df.index + 91
    df["Pos"] = df.index + int(set_["tm_start"])-1
    # reindex so that "Pos" is first
    df = df.reindex(columns = ["Pos", "AA", "Score"])
    # get the index of the top 7 predicted residues
    top7 = df.Score.sort_values(ascending=False).index[0:7]
    # label top 7 with "Y"
    df.loc[top7,"Int"] = "Y"
    # Put the Sine results into the dataframe (HERE I JUST USED RANDOM SINE NUMBERS)
    df["Sin"] = df.Pos.apply(lambda x : np.sin(x))
    
    # get the index of the top 10 predicted residues in the sine wave
    # count how many interface residues to show (here, a third of the total number of residues)
    number_of_sin_interface_res = int(df.shape[0]/ 3)
    # get the index of the sine interface residues
    sin_top_third = df.Sin.sort_values(ascending=False).index[0:number_of_sin_interface_res]
    # label Sine interface residues
    df.loc[sin_top_third, "SinInt"] = "Y"

    # convert any floats to 3 significant figures (e.g. 0.200)
    columns_to_convert = ["Score", "Sin"]
    for col in columns_to_convert:   
        df[col] = df[col].round(3)
        df[col] = df[col].apply(lambda x : "{:0.03f}".format(x))
    # save as a tab-separated CSV, which is easier for people to use with excel
    # end in .csv.txt, as most computers don't open csv files simply with a double-click
    # do NOT add the extra quotation marks, which prevent copying from txt to excel
    df.to_csv(test_out_csv, sep="\t")
    df