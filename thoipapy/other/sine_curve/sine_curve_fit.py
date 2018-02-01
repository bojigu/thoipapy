import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
#from __future__ import unicode_literals
#from thoipapy.sine_curve.tlabassays.mathfunctions import sine, sine_perfect_helix, residuals
from tlabassays.mathfunctions import sine, sine_perfect_helix, residuals
import scipy.optimize
import os
#from thoipapy.sine_curve.tlabtools import tools
#import tlabtools as tools
#from korbinian.utils import create_colour_lists
from thoipapy.utils import create_colour_lists
colour_lists = create_colour_lists()
fontsize = 12
data_dir = 'TMD16/'   #set the data_dir where you save all the data
file_dir_extension = os.path.join(data_dir, '*.xlsx')


def make_interpolater(left_min, left_max, right_min, right_max):
    # Figure out how 'wide' each range is
    leftSpan = left_max - left_min
    rightSpan = right_max - right_min

    # Compute the scale factor between left and right values
    scaleFactor = float(rightSpan) / float(leftSpan)

    # create interpolation function using pre-calculated scaleFactor
    def interp_fn(value):
        return right_min + (value - left_min) * scaleFactor

    return interp_fn


def save_sine_vurve_result(s, output_file_loc, output_png_loc):
    uniprot_acc=s["tm_protein_name"]
    #Pred_disruption=os.path.join(s["Sine_Curve_loc"],"TMD16/%s_prediction.xls") %uniprot_acc
    #Pred_disruption=r"/home/students/zeng/workspace/test2/out/58f795953dc45/output.csv"
    #Pred_disruption = r"sine_curve/TMD16/%s_prediction.xls" %uniprot_acc
    #prd_disruption = pd.read_csv(Pred_disruption)
    prd_disruption = pd.read_csv(output_file_loc)
    prd_disruption.dropna(inplace = True)
    # define the x-axis for the sine calculation as the contiguous index integer
    #x_sin = prd_disruption['Num']
    x_sin = prd_disruption['Res']
    ## set residue number and name as xaxis ticks
    indexamino=[]
    #amino=prd_disruption['Amino_acid']
    amino=prd_disruption['AA']
    i=0
    while i < len(x_sin):
        string=str(x_sin[i])+str(amino[i])
        indexamino.append(string)
        i+=1
    #y_sin_disruption = prd_disruption['Disruption'].values
    #y_sin_pred = prd_disruption['Prediction'].values
    y_sin_pred = prd_disruption['Score'].values

    #scaler1 = make_interpolater(min(y_sin_disruption), max(y_sin_disruption), 0.2, 0.8)
    # now convert to scaled values using map 
    #y_sin_disruption_scaled_data = list(map(scaler1, y_sin_disruption))

    scaler2 = make_interpolater(min(y_sin_pred), max(y_sin_pred), 0.2, 0.8)

    y_sin_pred_scaled_data=list(map(scaler2,y_sin_pred))

    plt.close("all")
    fig, ax = plt.subplots()
    x_smooth = np.linspace(x_sin.min(), x_sin.max(), 500)
    sine_constants_guess = [1.0,1.7,0.2,0]

    # sine_constants, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,
    #                                                                 sine_constants_guess,
    #                                                                 args=(sine,x_sin,y_sin_disruption_scaled_data),
    #                                                                 full_output=1)
    sine_constants1, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,
                                                                    sine_constants_guess,
                                                                    args=(sine,x_sin,y_sin_pred_scaled_data),
                                                                    full_output=1)

    #periodicity_fitted_curve = 2 * np.pi / sine_constants[1]
    periodicity_fitted_curve1 = 2 * np.pi / sine_constants1[1]
    #yvalues_fitted = sine(sine_constants, x_smooth)
    yvalues_fitted1 = sine(sine_constants1, x_smooth)
    sine_constants_guess_perfhelix = [1.0,0.5]
    # sine_constants_perfhelix, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = scipy.optimize.leastsq(residuals,
    #                                                                     sine_constants_guess_perfhelix,
    #                                                                     args=(sine_perfect_helix,x_sin,y_sin_disruption_scaled_data),
    #                                                                     full_output=1)
    sine_constants_perfhelix1, cov_perfhelix, infodict_perfhelix, mesg_perfhelix, ier_perfhelix = scipy.optimize.leastsq(residuals,
                                                                        sine_constants_guess_perfhelix,
                                                                        args=(sine_perfect_helix,x_sin,y_sin_pred_scaled_data),
                                                                        full_output=1)
    #yvalues_fitted_perfhelix = sine_perfect_helix(sine_constants_perfhelix, x_smooth)
    yvalues_fitted_perfhelix1 = sine_perfect_helix(sine_constants_perfhelix1, x_smooth)
    yvalues_fitted_perfhelix_pred = sine_perfect_helix(sine_constants_perfhelix1, x_sin)
    #yvalues_fitted_perfhelix_exp = sine_perfect_helix(sine_constants_perfhelix, x_sin)
    av_sin_pred=np.mean(yvalues_fitted_perfhelix_pred)
    #av_sin_exp=np.mean(yvalues_fitted_perfhelix_exp)
    y_pred_aa_norm=yvalues_fitted_perfhelix_pred -av_sin_pred
    #y_exp_aa_norm=yvalues_fitted_perfhelix_exp -av_sin_exp
    f=lambda x: 1 if x>0 else 0
    sine_peak_pred=map(f,y_pred_aa_norm)
    #sine_peak_exp=map(f,y_exp_aa_norm)
    sine_peak_pred_index=[]
    #sine_peak_exp_index=[]
    index_pred=0

    for x in sine_peak_pred:
        if(x==1):
            sine_peak_pred_index.append(index_pred)
        index_pred+=1
    # index_exp=0
    # for x in sine_peak_exp:
    #     if(x==1):
    #         sine_peak_exp_index.append(index_exp)
    #     index_exp+=1

    #pred_top5_index=np.argsort(y_sin_pred_scaled_data)[(len(y_sin_pred_scaled_data)-5):(len(y_sin_pred_scaled_data))]
    #exp_top5_index=np.argsort(y_sin_disruption_scaled_data)[(len(y_sin_disruption_scaled_data)-5):(len(y_sin_disruption_scaled_data))]
    #comele_pred_sine=list(set(sine_peak_pred_index).intersection(pred_top5_index))
    #comele_exp_sine=list(set(sine_peak_exp_index).intersection(exp_top5_index))
    sys.stdout.write(uniprot_acc)
    # print(sine_peak_pred)
    # print(comele_pred_sine)
    # print(comele_exp_sine)
    #ax.plot(x_sin, y_sin_disruption_scaled_data, color = 'orange',linestyle=":" ,linewidth=2.5,label = 'experimental data')
    ax.plot(x_sin, y_sin_pred_scaled_data, color = 'green',linewidth=2.5,linestyle=":", label = "predicted data")
    #ax.plot(x_smooth, yvalues_fitted_perfhelix, color = 'red',linewidth=2.5, linestyle="-", label = "experimental fit to sine")
    ax.plot(x_smooth, yvalues_fitted_perfhelix1, color = 'blue',linewidth=2.5, linestyle="-", label = "fitted data")
    plt.xticks(x_sin,amino)
    plt.ylim(0,1)

    #title="Sine Curve Comparison of %s" %uniprot_acc
    #ax.set_title(title) 
    plt.legend(loc='upper right', frameon=False)
    #ax.annotate(s=uniprot_acc,xy = (0.04,0.9),color = colour_lists['TUM_accents']['green'],fontsize=16, xytext=None, xycoords='axes fraction', alpha=1.00)
    #dir="TMD16\SineCurveRawData%s.png" %uniprot_acc
    #plt.savefig(dir)
    plt.grid()
    ax.tick_params(axis='x',direction='out', length=6, width=2,right='off',labeltop=(x_sin,amino))
    pred_top5_index=np.argsort(y_sin_pred_scaled_data)[(len(y_sin_pred_scaled_data)-5):(len(y_sin_pred_scaled_data))]
    i=0
    for tickname in ax.get_xticklabels():
        if(i in pred_top5_index):
            tickname.set_color('red')
        i+=1
    #plt.show()
    #plt.savefig('Results\Q7L4S7.png')
    plt.savefig(output_png_loc)
    plt.close()
    plt.close(fig)