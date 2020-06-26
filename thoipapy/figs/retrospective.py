import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pytoxr.mathfunctions import residuals, sine_perfect_helix
from scipy.optimize import leastsq
from scipy.stats import ttest_ind
from thoipapy.utils import create_colour_lists

def get_pivot_table_coev_data(s, i, XI, df_set):
    acc = df_set.loc[i, "acc"]
    database = df_set.loc[i, "database"]
    TMD_start = int(df_set.loc[i, "TMD_start"])
    TMD_end = int(df_set.loc[i, "TMD_end"])
    freecontact_file = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

    df = pd.read_csv(freecontact_file, sep=" ", header=None)
    df.columns = ["n1", "res1", "n2", "res2", "MI", "DI"]
    # due to uniprot indexing, residue_num should equal res + TMD_start - 1
    df["n1"] = df.n1 + TMD_start - 1
    df["n2"] = df.n2 + TMD_start - 1
    """df is a csv like this:
           n1 res1  n2 res2        MI        DI
    0  92    I  93    T  0.243618  0.454792
    1  92    I  94    L  0.404760  0.445580
    2  92    I  95    I  0.017704 -1.066260
    3  92    I  96    I  0.106223 -0.731704
    4  92    I  97    F  0.244482 -0.252246
    """

    dfp = df.pivot_table(index="n1", columns="n2", values=XI)

    """ asymmetrical pivoted data

         235       236       237       238       239       240 ...        252       253       254       255       256       
    n1                                                     ...                                                              
    235   0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523 
    236        NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577 
    237        NaN       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414 
    238        NaN       NaN       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822 
    """
    # get full list of residues
    position_list = range(TMD_start, TMD_end + 1)
    dfp = dfp.reindex(index=position_list, columns=position_list)
    # put data on both sides of the table for easy indexing
    for col in dfp.columns:
        start = col + 1
        dfp.loc[start:, col] = dfp.loc[col, start:]

    return dfp


def calc_coev_vs_res_dist(s, df_set, logging):
    """Calculate mean coevolution scores for each residue distance

    Plots the distance between two residues on the x-axis, against the mean coevolution scores.
    Calculated for both DI and MI.
    Does not take interface or non-interface definitions into account.

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

    Returns
    -------

    """

    logging.info('calc_coev_vs_res_dist starting')
    coev_vs_res_dist_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_coev_vs_res_dist.xlsx".format(s["setname"]))
    writer = pd.ExcelWriter(coev_vs_res_dist_xlsx)

    for XI in ["MI", "DI"]:
        nested_coev_dist_dict = {}
        nested_Cterm_dist_dict = {}

        for i in df_set.index:
            sys.stdout.write(".")
            sys.stdout.flush()
            acc_db = df_set.loc[i, "acc_db"]
            TMD_start = int(df_set.loc[i, "TMD_start"])
            TMD_end = int(df_set.loc[i, "TMD_end"])

            dfp = get_pivot_table_coev_data(s, i, XI, df_set)

            """dfp pivot table has a symmetric coevolution values between all residues
            
                    307       308       309       310       311       312       313       314       315       316    ...          326       327       328       329       330       331       332       333       334       335
            n1                                                                                                         ...                                                                                                       
            307       NaN  0.233388  0.251910  0.257193  0.365270  0.468933  0.313458  0.253943  0.278989  0.297606    ...     0.206634  0.153770  0.271118  0.364004  0.185186  0.286575  0.166321  0.313355  0.269962  0.307910
            308  0.233388       NaN  0.194896  0.290690  0.230827  0.300403  0.371423  0.149162  0.250657  0.283342    ...     0.130392  0.135035  0.317070  0.266557  0.134991  0.244770  0.164932  0.211624  0.185211  0.155684
            309  0.251910  0.194896       NaN  0.376993  0.401137  0.480202  0.313931  0.298846  0.336291  0.317149    ...     0.253053  0.229490  0.359081  0.366537  0.203667  0.264654  0.221240  0.373255  0.300027  0.240920
            310  0.257193  0.290690  0.376993       NaN  0.563525  0.651108  0.411248  0.353645  0.366177  0.455358    ...     0.264179  0.280039  0.398423  0.492760  0.226311  0.401377  0.242655  0.401358  0.316718  0.305772
            311  0.365270  0.230827  0.401137  0.563525       NaN  0.686111  0.367423  0.446554  0.457673  0.488534    ...     0.327568  0.286183  0.489552  0.468001  0.259249  0.447306  0.250619  0.455895  0.444784  0.345875

            """

            # set the "buffer", which controls the distance between the residues
            # will also be used to buffer the dataframe
            buffer = 20
            # extend axes by adding buffer
            rng = range(TMD_start - buffer, TMD_end + buffer + 1)
            dfp = dfp.reindex(index=rng, columns=rng)

            coev_dist_dict = {}
            N_coev_dist_dict = {}

            # iterate through each distance, b, in the buffer
            for b in range(1, buffer + 1):

                ###############################################
                #         Method 1 : mean of i+1 and i-1      #
                ###############################################
                single_dist_list = []
                for i in range(TMD_start, TMD_end + 1):
                    df_seln = dfp.loc[i, [i - b, i + b]]
                    """ The selection contains the two values, whose mean is taken.
                    Often, as b increases, one of the values is a nan.
                    [0.582511, nan]
                    [nan, 0.612737]
                    [nan, 0.5906779999999999]
                    """
                    mean_ = df_seln.mean()
                    assert isinstance(mean_, float)

                    if not np.isnan(mean_):
                       single_dist_list.append(mean_)
                    #single_dist_list.append(mean_)

                    # if np.isnan(mean_):
                    #     sys.stdout.write("({}.{}.{}.{})".format(b,i,seln, mean_))
                    #     sys.stdout.flush()

                ###############################################
                #     Method 2 : collect only i+1, i+2, etc   #
                ###############################################
                i_plus_b_coev_list = []
                for i in range(TMD_start, TMD_end - b + 1):
                    #i_plus_b is always a single coevolution value, between residue i, and residue i+b
                    i_plus_b = dfp.loc[i, i + b]
                    if not np.isnan(i_plus_b):
                        i_plus_b_coev_list.append(i_plus_b)
                    else:
                        # there should be no nans here, but check anyway.
                        sys.stdout.write("*")

                # get mean for both methods, for this TMD
                mean_for_this_dist = np.mean(single_dist_list)
                mean_Cterm = np.mean(i_plus_b_coev_list)

                # add mean for this distance(b) to summary dictionary
                if not np.isnan(mean_for_this_dist):
                    coev_dist_dict[b] = float("{:.03f}".format(mean_for_this_dist))
                if not np.isnan(mean_Cterm):
                    N_coev_dist_dict[b] = float("{:.03f}".format(mean_Cterm))

            # add summary values for this TMD to the output dicts
            nested_coev_dist_dict[acc_db] = coev_dist_dict
            nested_Cterm_dist_dict[acc_db] = N_coev_dist_dict

        # save output dicts with all TMD values to excel
        pd.DataFrame(nested_coev_dist_dict).to_excel(writer, sheet_name="coev_{}".format(XI))
        pd.DataFrame(nested_Cterm_dist_dict).to_excel(writer, sheet_name="C_{}".format(XI))

    writer.close()
    logging.info('calc_coev_vs_res_dist finished')


def plot_coev_vs_res_dist(s, logging):
    """Plot figure comparing coevolution values against residue distance.

    Uses excel output file from calc_coev_vs_res_dist.

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

    Returns
    -------

    """
    logging.info('plot_coev_vs_res_dist starting')
    plt.rcParams["font.family"] = "Verdana"
    plt.rcParams["font.family"] = "Verdana"
    colour_dict = create_colour_lists()
    blue1 = colour_dict["TUM_colours"]['TUM1']
    blue5 = colour_dict["TUM_colours"]['TUM5']
    TUMblue = colour_dict["TUM_colours"]['TUMBlue']

    fontsize = 9
    coev_vs_res_dist_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "{}_coev_vs_res_dist.xlsx".format(s["setname"]))

    fig, ax = plt.subplots(figsize=(4.5, 3.42))


    #######################################################################################################
    #                                                                                                     #
    #                                    add alpha helical sine wave                                      #
    #                                                                                                     #
    #######################################################################################################
    # [  0. ,   3.6,   7.2,  10.8,  14.4,  18. ,  21.6] for desired length
    x_helical = np.ogrid[0:23.4:14j]
    # [1,0, etc to train fitted sine as desired
    y_helical = [1,0]*7
    # fit to a perfect helix starting at 0 using leastsq
    sine_constants_guess_perfhelix = [1.575, 0.5]
    sine_constants_perfhelix1, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = leastsq(residuals,
                                                                        sine_constants_guess_perfhelix,
                                                                        args=(sine_perfect_helix,x_helical,y_helical),
                                                                        full_output=1)

    # create smooth sine curve x-points
    x_rng = np.linspace(0, 20, 200)
    logging.info("sine_constants_perfhelix1 : {}".format(sine_constants_perfhelix1))
    # adjust the height as desired
    sine_constants_fixed_centre = (sine_constants_perfhelix1[0], 0)
    # fit and plot
    yvalues_fitted_perfhelix1 = sine_perfect_helix(sine_constants_fixed_centre, x_rng)
    ax.plot(x_rng, yvalues_fitted_perfhelix1, color="0.8", linestyle="--", label=r"$\alpha$-helical periodicity")

    # only to show that heptad periodicity is not desirable
    plot_heptad_periodicity = False
    if plot_heptad_periodicity:
        # test of heptad motif
        x = range(0,21)
        y = np.array([1,0,0,1,1,0,1] * 3)/4
        ax.plot(x, y, color="r", alpha=0.5, linestyle=":", label="heptad periodicity")

    #######################################################################################################
    #                                                                                                     #
    #                                     add mean MI and DI values                                       #
    #                                                                                                     #
    #######################################################################################################

    # tab is "coev_" or "C_"
    # this corresponds to Method 1 or Method 2 above
    # they both give similar values, but i+b method is less likely to count values twice and is preferred
    excel_tab = "C_"
    #excel_tab = "coev_"

    # plot MI on secondary axis
    ax2 = ax.twinx()
    df = pd.read_excel(coev_vs_res_dist_xlsx, sheet_name="{}MI".format(excel_tab))
    mean_ser = df.mean(axis=1)
    mean_ser.plot(ax=ax2, label="mutual information (MI)", fontsize=fontsize, color=colour_dict["TUM_accents"]['orange'])

    df = pd.read_excel(coev_vs_res_dist_xlsx, sheet_name="{}DI".format(excel_tab))
    mean_ser = df.mean(axis=1)
    mean_ser.plot(ax=ax, label="direct information (DI)", fontsize=fontsize, color=TUMblue)

    # axis colours
    ax2.tick_params("y", colors=colour_dict["TUM_accents"]['orange'])
    ax.tick_params("y", colors=TUMblue)
    # axis labels
    ax.set_ylabel("mean DI coevolution score", labelpad=-3, fontsize=fontsize, color = TUMblue)
    ax2.set_ylabel("mean MI coevolution score", labelpad=1, fontsize=fontsize, color=colour_dict["TUM_accents"]['orange'])

    ax.set_xlabel("residue distance", fontsize=fontsize)

    ax.set_xticks(range(0, df.index.max()))
    # data gets messy due to low numbers after 15 residues.
    ax.set_xlim(0, 15)
    # add a bit of height, so the legend does not overlap data
    ax.set_ylim(mean_ser.min(), mean_ser.max() + 0.1)
    figpath = coev_vs_res_dist_xlsx[:-5] + "_coev" + ".png"
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    ax.legend(handles, labels, ncol=1, loc=1, fontsize=fontsize, frameon=True, facecolor='white')
    #fig.legend(fontsize=fontsize, loc="upper right", bbox_to_anchor=[0.85, 0.95])
    fig.tight_layout()
    fig.savefig(figpath, dpi=240)
    fig.savefig(figpath[:-4] + ".pdf")


def calc_retrospective_coev_from_list_interf_res(s, df_set, logging):

    logging.info('calc_retrospective_coev_from_list_interf_res starting')

    for randomise_int_res in [False, True]:

        if randomise_int_res == True:
            retrospective_coev_from_list_interf_res_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective", "{}_retrospective_coev_from_list_interf_res_random.xlsx".format(s["setname"]))
        else:
            retrospective_coev_from_list_interf_res_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective", "{}_retrospective_coev_from_list_interf_res.xlsx".format(s["setname"]))

        if not os.path.exists(os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective")):
            os.makedirs(os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective"))
        writer = pd.ExcelWriter(retrospective_coev_from_list_interf_res_xlsx)

        #randomise_int_res = False
        remove_residues_outside_interface_region = False
        logging.info("randomise_int_res = {}, remove_residues_outside_interface_region = {}".format(randomise_int_res, remove_residues_outside_interface_region))
        InterResList_of_last_TMD = None
        NoninterResList_of_last_TMD = None
        TMD_start_of_last_TMD = None

        for XI in ["MI", "DI"]:
            sub_dict = {}

            for i in df_set.index:
                sys.stdout.write(".")
                sys.stdout.flush()
                if i == 0:
                    is_first_TMD = True
                else:
                    is_first_TMD = False
                sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD = calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set,
                                                                                                                                                                  randomise_int_res, InterResList_of_last_TMD,
                                                                                                                                                                  NoninterResList_of_last_TMD, TMD_start_of_last_TMD,
                                                                                                                                                                  remove_residues_outside_interface_region)

            if randomise_int_res == True:
                # need to add the data for the first TMD, which was skipped above
                i = 0
                is_first_TMD = False
                sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD = calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set,
                                                                                                                                                                  randomise_int_res, InterResList_of_last_TMD,
                                                                                                                                                                  NoninterResList_of_last_TMD, TMD_start_of_last_TMD,
                                                                                                                                                                  remove_residues_outside_interface_region)
            # save each MI and DI separately
            df_retro = pd.DataFrame(sub_dict).T
            df_retro.to_excel(writer, sheet_name = XI)

            create_quick_plot = True
            if create_quick_plot:
                retrospective_coev_plot = retrospective_coev_from_list_interf_res_xlsx[:-5] + XI + ".png"
                df_retro["inter_larger"] = df_retro.AverageInter > df_retro.AverageNoninter
                fig, ax = plt.subplots()
                df_retro[["AverageInter", "AverageNoninter"]].plot(kind="bar", ax=ax)
                fig.tight_layout()
                fig.savefig(retrospective_coev_plot)

            # drop any (rare) rows without data, where the interface region was outside the length of the TMD?
            df_retro.dropna(how="any", inplace=True)

            vc = df_retro["inter_larger"].value_counts()
            if True in vc.index.tolist():
                n_TMDs_with_higher_int = vc[True]
            else:
                n_TMDs_with_higher_int = 0

            perc_higher_int = n_TMDs_with_higher_int / df_retro.shape[0]
            logging.info("\n{:.2f} % ({}/{}) of TMDs have higher {} of interface than non-interface".format(perc_higher_int*100, n_TMDs_with_higher_int, df_retro.shape[0], XI))

            logging.info("\n\nmean values\n{}\n".format(df_retro.mean()))

            t_value, p_value = ttest_ind(df_retro.AverageInter, df_retro.AverageNoninter)

            #logging.info("remove_residues_outside_interface_region = {}".format(remove_residues_outside_interface_region))
            logging.info("\np-value for average {} coevolution of interface vs non-interface = {:.03f}".format(XI, p_value))

        writer.close()

    sys.stdout.write("\n")
    logging.info('calc_retrospective_coev_from_list_interf_res finished')


def calc_retrospective_coev_from_list_interf_res_single_prot(sub_dict, s, logging, i, XI, is_first_TMD, df_set, randomise_int_res, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region):
    """Calculate average fraction of DI for a single protein.

    PANDAS METHOD USED FOR ETRA DATASET.

    SPLIT INTO A SUB-FUNCTION FOR USE DURING RANDOMISATION.

    Parameters
    ----------
    sub_dict : dict
        dictionary for each TMD. Separate dicts are made for MI and DI.
    s : dict
        Settings dictionary
    logging : logging.Logger
        Python object with settings for logging to console and file.
    i : int
        TMD number
    XI : str
        "MI" or "DI"
    is_first_TMD : bool
        whether the TMD is the first one, without initial randomised data
    df_set : pd.DataFrame
        Dataframe containing the list of proteins to process, including their TMD sequences and full-length sequences
        index : range(0, ..)
        columns : ['acc', 'seqlen', 'TMD_start', 'TMD_end', 'tm_surr_left', 'tm_surr_right', 'database',  ....]
    randomise_int_res : bool
        whether the interface residues shold be randomised
    InterResList_of_last_TMD : list
        List of interface residues from the last TMD. To be used during randomisation.
    NoninterResList_of_last_TMD : list
        List of non-interface residues from the last TMD. To be used during randomisation.
    TMD_start_of_last_TMD : int
        TMD start used to convert residue positions to range index

    Returns
    -------
    sub_dict : dict
        dictionary for each TMD. Separate dicts are made for MI and DI.
    InterResList_of_last_TMD : list
        List of interface residues from the last TMD. To be used during randomisation.
    NoninterResList_of_last_TMD : list
        List of non-interface residues from the last TMD. To be used during randomisation.
    TMD_start_of_last_TMD : int
        TMD start used to convert residue positions to range index
    """
    acc = df_set.loc[i, "acc"]
    database = df_set.loc[i, "database"]
    TMD_start = int(df_set.loc[i, "TMD_start"])
    TMD_end = int(df_set.loc[i, "TMD_end"])
    TMD_len = TMD_end - TMD_start
    freecontact_file = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr{}.gaps{}.freecontact.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))
    feature_combined_file = os.path.join(s["thoipapy_data_folder"], "Features", "combined", database, "{}.surr{}.gaps{}.combined_features.csv".format(acc, s["num_of_sur_residues"], s["max_n_gaps_in_TMD_subject_seq"]))

    df = pd.read_csv(freecontact_file, sep=" ", header=None)
    df.columns = ["n1", "res1", "n2", "res2", "MI", "DI"]
    # due to uniprot indexing, residue_num should equal res + TMD_start - 1
    df["n1"] = df.n1 + TMD_start - 1
    df["n2"] = df.n2 + TMD_start - 1
    """df is a csv like this:
           n1 res1  n2 res2        MI        DI
    0  92    I  93    T  0.243618  0.454792
    1  92    I  94    L  0.404760  0.445580
    2  92    I  95    I  0.017704 -1.066260
    3  92    I  96    I  0.106223 -0.731704
    4  92    I  97    F  0.244482 -0.252246
    """

    dfp = df.pivot_table(index="n1", columns="n2", values=XI)

    """ asymmetrical pivoted data

         235       236       237       238       239       240 ...        252       253       254       255       256       
    n1                                                     ...                                                              
    235   0.243618  0.404760  0.017704  0.106223  0.244482 ...   0.132235  0.219876  0.198667  0.360217  0.320984  0.145523 
    236        NaN  0.332451  0.140595  0.000747  0.151737 ...   0.217048  0.403469  0.174750  0.286540  0.357700  0.044577 
    237        NaN       NaN  0.062405  0.173925  0.353367 ...   0.336857  0.657512  0.418125  0.521322  0.538269  0.229414 
    238        NaN       NaN       NaN  0.049759  0.044692 ...   0.119658  0.236728  0.080722  0.114663  0.064796  0.096822 
    """
    # get full list of residues
    position_list = range(TMD_start, TMD_end + 1)
    dfp = dfp.reindex(index=position_list, columns=position_list)
    # put data on both sides of the table for easy indexing
    for col in dfp.columns:
        start = col + 1
        dfp.loc[start:, col] = dfp.loc[col, start:]


    """dfp now contains the coevolution data as a symmetrical dataframe, from each residue to the other

              192       193       194       195       196       197       198       199       200       201    ...          210       211       212       213       214       215       216       217       218       219
    n1                                                                                                         ...                                                                                                       
    192       NaN  0.092237  0.026186  0.126701  0.108622  0.107383  0.075048  0.070287  0.084822  0.037957    ...     0.152848  0.074908  0.073767  0.159693  0.044335  0.092576  0.057039  0.176549  0.076715  0.066157
    193  0.092237       NaN  0.089528  0.137392  0.112203  0.153103  0.114659  0.173971  0.134006  0.091982    ...     0.237441  0.107704  0.097004  0.216488  0.146309  0.100271  0.101273  0.301949  0.105543  0.193257
    194  0.026186  0.089528       NaN  0.102470  0.078647  0.138274  0.141817  0.142261  0.133799  0.079009    ...     0.172375  0.111071  0.121039  0.171232  0.106160  0.095982  0.188747  0.230212  0.093526  0.217379
    195  0.126701  0.137392  0.102470       NaN  0.162021  0.124095  0.131162  0.248673  0.167416  0.094939    ...     0.179470  0.168239  0.139384  0.193543  0.102942  0.172607  0.153524  0.289339  0.113594  0.181711
    196  0.108622  0.112203  0.078647  0.162021       NaN  0.147395  0.106920  0.186598  0.170876  0.074893    ...     0.152920  0.130958  0.104620  0.165248  0.071461  0.117822  0.113831  0.243438  0.097208  0.153550
    197  0.107383  0.153103  0.138274  0.124095  0.147395       NaN  0.185372  0.300418  0.254464  0.135116    ...     0.294558  0.214323  0.237466  0.396039  0.111643  0.203568  0.221890  0.442481  0.167183  0.255704
    198  0.075048  0.114659  0.141817  0.131162  0.106920  0.185372       NaN  0.188028  0.174667  0.158833    ...     0.145839  0.134066  0.147938  0.256873  0.098789  0.146614  0.202526  0.266566  0.114003  0.211277
    199  
    """

    # open combined file with interface definitions
    dfc = pd.read_csv(feature_combined_file, index_col=0)

    # set the residue numbering as the index
    dfc.index = dfc.res_num_full_seq.astype(int)
    # get list of interface and noninterface residue positions
    InterResList = dfc.loc[dfc.interface == 1].index
    NoninterResList = list(dfc.loc[dfc.interface == 0].index)

    # get dataframe of distances between the residues


    if remove_residues_outside_interface_region:
        #logging.info("orig NoninterResList = {}".format(NoninterResList))
        lowest_interface_res = InterResList.min()
        highest_interface_res = InterResList.max()
        #logging.info("interface residue range = {} to {}".format(lowest_interface_res, highest_interface_res))
        NoninterResList = [x for x in NoninterResList if lowest_interface_res < x < highest_interface_res]
        #logging.info("final NoninterResList = {}".format(NoninterResList))

    if randomise_int_res == False:
        # calculate mean coevolution values for the desired selection of the dataframe.
        # Note that the values are symmetric and doubled ([235,236] and also [236,235])
        # but the mean will be unaffected
        mean_XI_interface = dfp.loc[InterResList, InterResList].mean().mean()
        mean_XI_noninterface = dfp.loc[NoninterResList, NoninterResList].mean().mean()

        sub_dict[acc] = {"AverageInter": mean_XI_interface, "AverageNoninter": mean_XI_noninterface}
    elif randomise_int_res == True and is_first_TMD == True:
        # can't used the interface residues from previous protein
        pass
    elif randomise_int_res == True and is_first_TMD != True:
        # FOR RANDOMISATION, THE ORIGINAL INDEXING BASED ON AA NUMBERS MUST BE REPLACED BY RANGE INDEXING

        # convert from amino acid numbering to a range index by subtracting the TMD_start
        InterResList_of_last_TMD = pd.Series(InterResList_of_last_TMD) - TMD_start_of_last_TMD
        NoninterResList_of_last_TMD = pd.Series(NoninterResList_of_last_TMD) - TMD_start_of_last_TMD
        # drop any residue positions that are longer than the TMD (i.e., where index TMD is longer than TMD supplying coevolution data)
        InterResList_of_last_TMD = InterResList_of_last_TMD[InterResList_of_last_TMD <= TMD_len].tolist()
        NoninterResList_of_last_TMD = NoninterResList_of_last_TMD[NoninterResList_of_last_TMD <= TMD_len].tolist()

        # convert pivoted data to range index
        dfp.index = range(dfp.shape[0])
        dfp.columns = range(dfp.shape[1])

        # slice out the interface and non-interface residues, as above
        mean_XI_interface = dfp.loc[InterResList_of_last_TMD, InterResList_of_last_TMD].mean().mean()
        mean_XI_noninterface = dfp.loc[NoninterResList_of_last_TMD, NoninterResList_of_last_TMD].mean().mean()

        sub_dict[acc] = {"AverageInter": mean_XI_interface, "AverageNoninter": mean_XI_noninterface}

    InterResList_of_last_TMD = InterResList
    NoninterResList_of_last_TMD = NoninterResList
    TMD_start_of_last_TMD = TMD_start

    return sub_dict, InterResList_of_last_TMD, NoninterResList_of_last_TMD, TMD_start_of_last_TMD


def get_array_dist_separating_res_in_list(pos_list):
    """Get array of the distance separating the residues in a list.

    Parameters
    ----------
    pos_list : list
        list of positions that are interacting, e.g [28,29,32,33,34,36]

    Returns
    -------
    dist_2D_arr np.ndarray
        numpy array showing the distances between all combinations in the list of residues.
        E.g.
        array([[ nan,   1.,   4.,   5.,   6.,   8.],
               [  1.,  nan,   3.,   4.,   5.,   7.],
               [  4.,   3.,  nan,   1.,   2.,   4.],
               [  5.,   4.,   1.,  nan,   1.,   3.],
               [  6.,   5.,   2.,   1.,  nan,   2.],
               [  8.,   7.,   4.,   3.,   2.,  nan]])

        the dataframe equivalent looks like this:
             28   29   32   33   34   36
        28  NaN  1.0  4.0  5.0  6.0  8.0
        29  1.0  NaN  3.0  4.0  5.0  7.0
        32  4.0  3.0  NaN  1.0  2.0  4.0
        33  5.0  4.0  1.0  NaN  1.0  3.0
        34  6.0  5.0  2.0  1.0  NaN  2.0
        36  8.0  7.0  4.0  3.0  2.0  NaN

    """
    # make empty array
    le = len(pos_list)
    dist_2D_arr = np.zeros(le*le).reshape(le, le)

    # calculate distances for each position separately. Add to array
    for i, pos in enumerate(pos_list):
        arr = np.array(pos_list)
        arr = abs(arr - pos)
        dist_2D_arr[i, :] = arr
    # replace 0 with nan, so that it doesn't count in the averages
    dist_2D_arr[dist_2D_arr == 0] = np.nan
    return dist_2D_arr


def calc_retrospective_coev_from_struct_contacts(s, dfset, logging):
    """Calculate retrospective coevolution from structural data with interpair contacts

    Calculates coevolution of interface (contacting) residues from homodimer structures
    using the method of Wang and Barth (2015).

    Takes csv with list of interacting residues as an input:
    acc	    inter1	inter2
    1orqC4	14	    20
    1orqC4	14	    23
    1orqC4	14	    24
    1orqC4	18	    24
    1xioA4	1	    2
    1xioA4	5	    6

    Converts this to a list of interacting pairs:
     - saved as InterPairList by the calc_retrospective_coev_from_struct_contacts_single_prot function
     - this is named as "last TMD" because the randomisation takes took? the positions of first TMD and applied them to the last

    Takes pairwise coevolution values from FreeContact output file.
    E.g. D:\Dropbox\tm_homodimer_dropbox\THOIPA_data\Features\cumulative_coevolution\ NMR\O15455.surr20.gaps5.freecontact.csv
    Data looks like this:
    1 F 2 F 0.195863 0.552187
    1 F 3 M 0.172853 -0.530669
    1 F 4 I 0.406909 0.445122
    1 F 5 N 0.245805 2.64414

    Finds coevolution scores for all interacting residue pairs.
     - calculates average

    Finds coevolution scores for all non-interacting pairs, separated by no longer than 8 residues.
     - creates list of non-interacting pairs (NonInterPairList)

    REPEATS THE ABOVE WITH RANDOM POSITIONS, TAKEN FROM A DIFFERENT TMD IN THE LIST OF PROTEINS/TMDS TO TEST
     - sorts the protein list by length of TMD (ascending, shorter TMDs first)
     - for each protein, saves contact positions from last TMD for use in calculating mean coevolution, etc
     - note that the output InterPairList is the real one, which is then applied to the next TMD

    FOR THE FIRST TMD
     - skipped in first round of iteration as no positions from prev TMD available
     - re-run separately afterwards
     - ANY CONTACT POSITIONS FROM LAST TMD THAT ARE LONGER THAN YIDNC TMD ARE EXCLUDED! (see InterPairList_rand creation)


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

    Returns
    -------

    """
    logging.info('calc_retrospective_coev_from_struct_contacts starting')

    # if randomise_int_res == True:
    #     retrospective_coev_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective", "{}_retrospective_coev_from_struct_contact_random.xlsx".format(s["setname"]))
    # else:
    #     retrospective_coev_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective", "{}_retrospective_coev_from_struct_contact.xlsx".format(s["setname"]))

    retrospective_coev_xlsx = os.path.join(s["thoipapy_data_folder"], "Results", s["setname"], "retrospective", "{}_retrospective_coev_from_struct_contact.xlsx".format(s["setname"]))
    writer = pd.ExcelWriter(retrospective_coev_xlsx)

    for randomise_int_res in [False, True]:
        # suffix for excel tabs e.g. "DI" and "DI_rand"
        suffix = "_rand" if randomise_int_res == True else ""

        remove_residues_outside_interface_region = False
        logging.info("randomise_int_res = {}, remove_residues_outside_interface_region = {}".format(randomise_int_res, remove_residues_outside_interface_region))
        InterPairList_of_last_TMD = None
        NonInterPairList_of_last_TMD = None
        TMD_start_of_last_TMD = None
        crystal_NMR_interpair_file = os.path.join(s["thoipapy_data_folder"], "Results", "Average_Fraction_DI", "Crystal_NMR_interpair.csv")
        pd_int = pd.read_csv(crystal_NMR_interpair_file, engine="python")
        """pd_int looks like this    
        acc	    inter1	inter2
        1orqC4	14	    20
        1orqC4	14	    23
        1orqC4	14	    24
        1orqC4	18	    24
        1xioA4	1	    2
        1xioA4	5	    6
        """

        # get list of proteins sort by TMD length
        # note that this will probably only work for set04
        dfset.index = dfset['TMD_seq'].str.len()
        #sys.stdout.write("test of random sorting")
        #dfset.index = np.random.random(dfset.shape[0])
        dfset = dfset.sort_index(ascending=True).reset_index(drop=True)

        sub_dict = {}
        NoninterPairList_dict = {}
        InterPairList_dict = {}

        for i in dfset.index:
            acc = dfset.at[i, "acc"]
            TMD_len = len(dfset.at[i, "TMD_seq"])
            sys.stdout.write(".")
            sys.stdout.flush()
            if i == 0:
                is_first_TMD = True
            else:
                is_first_TMD = False

            #if not all([randomise_int_res == True and i == 0]):

            #is_first_TMD = False
            sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD = calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, dfset, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region)
            # save all lists of interacting and non-interacting pairs
            InterPairList_dict[acc] = str(InterPairList)
            NoninterPairList_dict[acc] = str(NoninterPairList)

            InterPairList_of_last_TMD = InterPairList
            #sys.stdout.write("IN LOOP InterPairList_of_last_TMD", InterPairList_of_last_TMD)
            #sys.stdout.write("acc", acc, "IN LOOP TMD_len", len(dfset.at[i, "TMD_seq"]))
            #InterPairList_of_last_TMD_max = np.array(InterPairList_of_last_TMD).max()
            #NonInterPairList_of_last_TMD_max = np.array(NonInterPairList_of_last_TMD).max()
            NonInterPairList_of_last_TMD = NoninterPairList


        if randomise_int_res == True:
            # need to add the data for the first TMD, which was skipped above
            i = 0
            is_first_TMD = False
            logging.info("Randomisation is True. All other proteins finished. Starting first TMD, using the contacts from the last (longest) TMD above.")
            sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD = calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, dfset, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region)

        df_retro = pd.DataFrame(sub_dict).T
        df_retro.to_excel(writer, sheet_name="DI{}".format(suffix))

        """Save lists of interacting or non-interacting residue pairs:
        E.g. 
                InterPair	                                 NonInterPair
        1orqC4	[[14, 20], [14, 23], [14, 24], [18, 24]]	[[1, 2], [2, 1], [1, 3], [3, 1], [1, 4], [4, 1], [1, 5], [5, 1], ..........
        1xioA4	[[1, 2], [5, 6], [12, 12], [16, 16]]	    [[1, 3], [3, 1], [1, 4], [4, 1], [1, 5], [5, 1], [1, 6], [6, 1], .............
        """
        if not randomise_int_res:
            InterPair_ser = pd.Series(InterPairList_dict)
            NonInterPair_ser = pd.Series(NoninterPairList_dict)
            df_res_list = pd.DataFrame([InterPair_ser, NonInterPair_ser], index=["InterPair", "NonInterPair"]).T
            df_res_list.to_excel(writer, sheet_name="PairLists")

        create_quick_plot = True
        if create_quick_plot:
            retrospective_coev_plot = retrospective_coev_xlsx[:-5]  + "DI{}.png".format(suffix)
            df_retro["inter_larger"] = df_retro.AverageInter > df_retro.AverageNoninter
            fig, ax = plt.subplots()
            df_retro[["AverageInter", "AverageNoninter"]].plot(kind="bar", ax=ax)
            fig.tight_layout()
            fig.savefig(retrospective_coev_plot)

        # drop any (rare) rows without data, where the interface region was outside the length of the TMD?
        df_retro.dropna(how="any", inplace=True)

        vc = df_retro["inter_larger"].value_counts()
        if True in vc.index.tolist():
            n_TMDs_with_higher_int = vc[True]
        else:
            n_TMDs_with_higher_int = 0

        perc_higher_int = n_TMDs_with_higher_int / df_retro.shape[0]
        logging.info("\n{:.2f} % ({}/{}) of TMDs have higher DI of interface than non-interface".format(perc_higher_int * 100,
                                                                                               n_TMDs_with_higher_int,
                                                                                               df_retro.shape[0]))

        logging.info("\n\nmean values\n{}\n".format(df_retro.mean()))

        t_value, p_value = ttest_ind(df_retro.AverageInter, df_retro.AverageNoninter)

        # logging.info("remove_residues_outside_interface_region = {}".format(remove_residues_outside_interface_region))
        logging.info("\np-value for average DI coevolution of interface vs non-interface = {:.03f}".format( p_value))



    writer.close()
    sys.stdout.write("\n")
    logging.info('calc_retrospective_coev_from_struct_contacts finished')


def calc_retrospective_coev_from_struct_contacts_single_prot(sub_dict, s, pd_int, i, logging, is_first_TMD, dfset, randomise_int_res, InterPairList_of_last_TMD, NonInterPairList_of_last_TMD, TMD_start_of_last_TMD, remove_residues_outside_interface_region):
    """Calculate retrospective coevolution from structural data with interpair contacts, for a single protein.

    see docstring above for calc_retrospective_coev_from_struct_contacts

    Parameters
    ----------
    sub_dict
    s
    pd_int
    i
    logging
    is_first_TMD
    dfset
    randomise_int_res
    InterPairList
    NonInterPairList
    TMD_start_of_last_TMD
    remove_residues_outside_interface_region

    Returns
    -------

    """
    # IMPORTANT!
    # hard-coded max distance between non-contacting residue pairs
    # Wang and Barth used a value of 8. However this leads to a completely different distribution of contacting and non-contacting residues.
    max_dist_noncontact_res_pairs = 8

    acc = dfset.loc[i,'acc']
    database = dfset.loc[i,'database']
    InterPairList = pd_int[["inter1", "inter2"]][pd_int["acc"] == acc].values.tolist()
    interlist = []
    for x in InterPairList:
        interlist.extend(x)
    lowest_interface_res = min(interlist)
    highest_interface_res = max(interlist)
    NoninterPairList = []
    freecontact_file = os.path.join(s["thoipapy_data_folder"], "Features", "cumulative_coevolution", database, "{}.surr20.gaps5.freecontact.csv".format(acc))
    inter_within8_dict = {}
    DI_dict = {}


    if os.path.isfile(freecontact_file):
        with open(freecontact_file, 'r') as f:

            for line in f:
                # arr: < class 'list'>: ['16', 'V', '17', 'E', '0.289263', '5.95475']
                arr = line.strip().split()
                tmd_len = int(arr[2])
                # add forward and backward pair to DI_dict (e.g. 16_17, or 17_16)
                # DI_dict and inter_within8_dict looks like this
                # {'1_2': '3.12626', '2_1': '3.12626', '1_3': '0.589193', '3_1': '0.589193',
                DI_dict[arr[0] + '_' + arr[2]] = arr[5]
                DI_dict[arr[2] + '_' + arr[0]] = arr[5]

                # if separating distance is less than max_dist_noncontact_res_pairs, add pair scores forwards and backwards
                if int(arr[2]) - int(arr[0]) <= max_dist_noncontact_res_pairs:
                    #arr[5] is the DI score
                    inter_within8_dict[arr[0] + '_' + arr[2]] = arr[5]
                    inter_within8_dict[arr[2] + '_' + arr[0]] = arr[5]
                    # if the pair (e.g. 16_17, or 17_16) is not in the list of interface pairs
                    # add it to the list of non-interface residues
                    if [int(arr[0]), int(arr[2])] not in InterPairList and [int(arr[0]), int(arr[2])] not in InterPairList:
                        if remove_residues_outside_interface_region:
                            if lowest_interface_res < int(arr[0]) < highest_interface_res and lowest_interface_res < int(arr[2]) < highest_interface_res:
                                NoninterPairList.append([int(arr[0]), int(arr[2])])
                                NoninterPairList.append([int(arr[2]), int(arr[0])])
                        else:
                            NoninterPairList.append([int(arr[0]), int(arr[2])])
                            NoninterPairList.append([int(arr[2]), int(arr[0])])
        f.close()


    if randomise_int_res == False:
        logging.info("No randomisation. Real InterPairList for this TMD = {}".format(InterPairList))

        average_DI_inter, average_DI_non_inter = calc_average_DI_inter_and_average_DI_non_inter(InterPairList, NoninterPairList, DI_dict, inter_within8_dict)

        #sys.stdout.write(non_inter_DI)
        #sys.stdout.write(acc,average_DI_inter,average_DI_non_inter, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)
        sub_dict[acc] = {"AverageInter": average_DI_inter, "AverageNoninter": average_DI_non_inter}

    # Randomisation of first TMD not possible. Skip and return the true interface only.
    elif randomise_int_res == True and is_first_TMD == True:
        # First TMD is skipped now, and added separately later
        # No InterPairList_of_last_TMD available for interface randomisation!
        pass

    elif randomise_int_res == True and is_first_TMD != True:
        # TEST: add 2 residues to all interface residues
        #InterPairList_of_last_TMD =[list(np.array(x) + 5) for x in InterPairList_of_last_TMD]

        # delete any interpair positions from a previous TMD that are out of range for index of THIS TMD
        # should only apply to the first TMD, which is processed at the end
        InterPairList_rand =[x for x in InterPairList_of_last_TMD if x[0] <= tmd_len and x[1] <= tmd_len ]
        NonInterPairList_rand = [x for x in NonInterPairList_of_last_TMD if x[0] <= tmd_len and x[1] <= tmd_len ]

        logging.info("Randomisation of contacting residues for TMD # {} (acc = {})".format(i, acc))
        logging.info("Real InterPairList for this TMD = {}".format(InterPairList))
        logging.info("Random NonInterPairList_rand for TMD, taken from last TMD in list = {}".format(InterPairList_rand))

        orig_max = np.array(InterPairList_of_last_TMD).max()
        final_max_after_filtering_out_res_too_long = np.array(InterPairList_rand).max()

        if orig_max != final_max_after_filtering_out_res_too_long:
            sys.stdout.write("\nmust be the last TMD, if sorting is done correctly")
            sys.stdout.write("\norig_max {} final_max_after_filtering_out_res_too_long {}".format(orig_max, final_max_after_filtering_out_res_too_long))

        # #sys.stdout.write(acc,InterPairList_rand,NonInterPairList_rand, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)
        # inter_DI = []
        # non_inter_DI = []
        # for key, value in inter_within8_dict.items():
        #     inter_pair = [int(x) for x in key.split('_')]
        #     if inter_pair in InterPairList_rand:
        #         inter_DI.append(float(inter_within8_dict[key]))
        #     elif inter_pair in NonInterPairList_rand:
        #         non_inter_DI.append(float(inter_within8_dict[key]))
        # average_DI_inter = np.mean(inter_DI)
        # average_DI_non_inter = np.mean(non_inter_DI)
        # #sys.stdout.write(acc,average_DI_inter,average_DI_non_inter, "is first TMD", is_first_TMD, "randomise_int_res", randomise_int_res)

        average_DI_inter, average_DI_non_inter = calc_average_DI_inter_and_average_DI_non_inter(InterPairList_rand, NonInterPairList_rand, DI_dict, inter_within8_dict)

        sub_dict[acc] = {"AverageInter": average_DI_inter, "AverageNoninter": average_DI_non_inter}

    #InterPairList = InterPairList
    #NonInterPairList = NoninterPairList


    return sub_dict, InterPairList, NoninterPairList, TMD_start_of_last_TMD


def calc_average_DI_inter_and_average_DI_non_inter(InterPairList, NoninterPairList, DI_dict, inter_within8_dict):
    """Calculate average DI of interacting and non-interacting residues

    Parameters
    ----------
    InterPairList : list
        List of interacting residue pairs
    NoninterPairList : list
        List of non-interacting residue pairs
    DI_dict : dict
        Dict of all DI coevolution values for all residues
    inter_within8_dict : dict
        Dect of all DI coevolution values, excluding those separated by more than max_dist_noncontact_res_pairs

    Returns
    -------
    average_DI_inter : float
        Average DI of interacting residues
    average_DI_non_inter : float
        Average DI of non-interacting residues, separated by a max of max_dist_noncontact_res_pairs (typically 4 residues)
    """
    inter_DI = []
    non_inter_DI = []
    # DEPRECATED OLD CODE ALSO APPLIED 8 residue limit to the interacting residues!!!
    # for key, value in inter_within8_dict.items():
    #     inter_pair = [int(x) for x in key.split('_')]
    #     if inter_pair in InterPairList:
    #         inter_DI.append(float(inter_within8_dict[key]))
    #     elif inter_pair in NoninterPairList:
    #         non_inter_DI.append(float(inter_within8_dict[key]))

    # interface residues have no filter
    for key, value in DI_dict.items():
        inter_pair = [int(x) for x in key.split('_')]
        if inter_pair in InterPairList:
            inter_DI.append(float(DI_dict[key]))

    # non-interface residues filtered by max_dist_noncontact_res_pairs
    for key, value in inter_within8_dict.items():
        inter_pair = [int(x) for x in key.split('_')]
        if inter_pair in NoninterPairList:
            non_inter_DI.append(float(inter_within8_dict[key]))

    average_DI_inter = np.nanmean(inter_DI)
    average_DI_non_inter = np.nanmean(non_inter_DI)
    len_non_inter_DI = len(non_inter_DI)
    if len_non_inter_DI == 0 or np.isnan(average_DI_inter):
        sys.stdout.write(inter_DI)
        sys.stdout.write(len(non_inter_DI))
        sys.stdout.write(non_inter_DI)
        sys.stdout.write(NoninterPairList)
        sys.stdout.write(InterPairList)
        sys.stdout.write(pd.Series(non_inter_DI).isnull().value_counts())

    return average_DI_inter, average_DI_non_inter