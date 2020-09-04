from weighslide import calculate_weighted_windows

from thoipapy.utils import normalise_0_1
import numpy as np


def normalise_features(df_features_single_protein):
    """Normalise selected features.

    Number of homologues -> normalised between 1 and 8
    conservation -> as minus entropy
    LIPS L*E -> calculated here
    LIPS surface ranked -> calculated based on LIPS_surface and LIPS L*E
    CS, DE, etc -> %C+%S, %D+%E, etc

    Parameters
    ----------
    df_features_single_protein : pd.DataFrame
        Dataframe with all features for a protein
        Index : range index
        Columns : "residue_num", "residue_name", "entropy", etc
    """
    # calculate LIPS L*E for later validation. LIPS_L*E = LIPS_polarity * log(LIPS_entropy)
    df_features_single_protein["LIPS_L*E"] = df_features_single_protein.LIPS_polarity * np.log(df_features_single_protein.LIPS_entropy)
    # rank the LIPS score by adding a fraction of the L*E to the predicted interface (0 or 1)
    df_features_single_protein["LIPS_surface_ranked"] = df_features_single_protein["LIPS_surface"] - (df_features_single_protein["LIPS_L*E"] / 20)
    df_features_single_protein["LIPS_surface_ranked_norm"] = normalise_0_1(df_features_single_protein["LIPS_surface_ranked"])[0]

    # create our own conservation + polarity and conservation*polarity
    df_features_single_protein["cons+polarity"] = df_features_single_protein["conservation"] + df_features_single_protein["polarity"]
    df_features_single_protein["cons*polarity"] = df_features_single_protein["conservation"] * df_features_single_protein["polarity"]

    # add the mean polarity or conservation of positions i, i+4 and i-4
    window = [1,0,0,0,1,0,0,0,1]
    df_features_single_protein["cons4mean"] = calculate_weighted_windows(df_features_single_protein["conservation"], window, statistic="mean", full_output=False)
    df_features_single_protein["polarity4mean"] = calculate_weighted_windows(df_features_single_protein["polarity"], window, statistic="mean", full_output=False)

    df_features_single_protein["CS"] = df_features_single_protein["C"] + df_features_single_protein["S"]
    df_features_single_protein["DE"] = df_features_single_protein["D"] + df_features_single_protein["E"]
    df_features_single_protein["KR"] = df_features_single_protein["K"] + df_features_single_protein["R"]
    df_features_single_protein["QN"] = df_features_single_protein["Q"] + df_features_single_protein["N"]
    df_features_single_protein["LIV"] = df_features_single_protein["L"] + df_features_single_protein["I"] + df_features_single_protein["V"]

    # round the relative positions so there are 10 options (to 10%)
    df_features_single_protein["RelPos_TMD"] = df_features_single_protein["RelPos_TMD"].round(1)
    df_features_single_protein["RelPos_fullseq"] = df_features_single_protein["RelPos_fullseq"].round(1)

    return df_features_single_protein


def normalise_number_of_homologues_to_categorical_variable(x):
    """ DEPRECATED. Replaced with cubed root of the n_homologues.
    Convert non-linear number of homologues to an integer value.

    Parameters
    ----------
    x : int
        Number of homologues in fasta alignment of TMD region including gaps.
    """
    if x <= 75:
        return 1
    elif x <= 100:
        return 2
    elif x <= 200:
        return 3
    elif x <= 400:
        return 4
    elif x <= 800:
        return 5
    elif x <= 1600:
        return 6
    elif x <= 3200:
        return 7
    else:
        return 8