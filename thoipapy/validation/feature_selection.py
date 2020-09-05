import pandas as pd

from thoipapy.utils import convert_truelike_to_bool, convert_falselike_to_bool


def drop_cols_not_used_in_ML(logging, df_data, excel_file_with_settings, i=0):
    """Remove columns not used in machine learning training or testing.

    This includes
        1) text and name columns (e.g. residue_name)
        2) bind columns (interface, interface_score).
        3) columns of non-normalised features, or features that we don't want to include anymore.

    Parameters
    ----------
    logging : logging.Logger
        Python object with settings for logging to console and file.
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns = 'acc_db', "residue_num", "residue_name", etc
        rows = range of number of residues
    excel_file_with_settings : str
        Path to excel file with all settings.
        Necessary for getting the list of features included in the THOIPA algorithm.


    Returns
    -------
    df_data : pd.DataFrame
        Dataframe with either the training or test dataset
        columns :
    """
    # read the features tab of the excel settings file
    features_df = pd.read_excel(excel_file_with_settings, sheet_name="features", index_col=0)
    # convert "WAHR" etc to true and false
    features_df["include"] = features_df["include"].apply(convert_truelike_to_bool, convert_nontrue=False)
    features_df["include"] = features_df["include"].apply(convert_falselike_to_bool)

    # print any features that are not shared between the two lists
    unused_features_in_excel = set(features_df.index) - set(df_data.columns)
    features_missing_from_excel = set(df_data.columns) - set(features_df.index)

    # only print this stuff for the first TMD analysed, if in a list
    if i == 0:
        if len(features_missing_from_excel) > 1:
            logging.info("\nfeatures_missing_from_excel, {}".format(features_missing_from_excel))
            if "Unnamed: 0" in features_missing_from_excel:
                raise IndexError("Unnamed column is in dataframe. Try opening csv with index_col=0.")
        if len(unused_features_in_excel) > 1:
            logging.info("\nunused_features_in_excel, {}".format(unused_features_in_excel))
            if len(unused_features_in_excel) > 2:
                logging.warning("There are more than two unused features in the excel settings file. "
                                "This is acceptable for standalone predictor, but otherwise the settings file "
                                "might need to be checked.")
    # drop any features that are not labeled TRUE for inclusion
    features_df = features_df.loc[features_df.include == True]
    # filter df_data to only keep the desired feature columns
    feature_list = features_df.index.tolist()
    feature_list_shared_cols = list(set(feature_list).intersection(set(df_data.columns)))
    df_data = df_data.loc[:, feature_list_shared_cols]
    return df_data