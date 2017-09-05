import pandas as pd
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import StratifiedKFold
import os
import subprocess, threading


def run_Rscipt_random_forest(tmplist,pathdict,set_,out_put_file_loc,logging):
    logging.info('begining to run random forest R code')
    Rscript_loc = set_["Rscript_dir"]
    Random_Forest_R_code_file=set_["Rcode"]
    train_data_file=os.path.join(set_["RF_loc"],"NoRedundPro/TRAINDATA68.csv")
    tmp_protein_name = set_["tm_protein_name"]
    tmp_protein_test_data = os.path.join(set_["RF_loc"], "TestData/%s/%s.mem.2gap.physipara.testdata.csv") % (set_["Datatype"],tmp_protein_name)
    #out_put_file_loc_handle=open(out_put_file_loc,"w")
    if os.path.isfile(tmp_protein_test_data):
        prediction_output_file = os.path.join(set_["RF_loc"],"%s.pred.out") % tmp_protein_name
        prediction_output_file = os.path.join("/home/students/zeng/workspace/test2/out", "%s.pred.out") % tmp_protein_name
        prediction_output_file_handle=open(prediction_output_file,"w")
        exect_str = "{Rscript} {Random_Forest_R_code} {train_data} {test_data} {output}".format(Rscript=Rscript_loc, Random_Forest_R_code=Random_Forest_R_code_file,train_data=train_data_file,test_data=tmp_protein_test_data,output=out_put_file_loc)
        print(exect_str)
        class Command(object):
            def __init__(self, cmd):
                self.cmd = cmd
                self.process = None

            def run(self, timeout):
                def target():
                    print('Thread started')
                    self.process = subprocess.Popen(self.cmd,shell=True,stdout=subprocess.PIPE)
                    #subprocess.call(self.cmd,shell=True,stdout=prediction_output_file_handle)
                    subprocess.call(self.cmd, shell=True)
                    #self.process.communicate()
                    # print(self.process.communicate())
                    print('Thread finished')

                thread = threading.Thread(target=target)
                thread.start()

                thread.join(timeout)
                if thread.is_alive():
                    print('Terminating process')
                    self.process.terminate()
                    thread.join()
                print(self.process.returncode)

        command = Command(exect_str)
        # command=Command(exect_str)
        command.run(timeout=2000)                       ###since hhblits requres more than 10 minutes to finish, maybe later we should consider using qsub to the server
        # command.run(timeout=1)
        # command=mtutils.Command(exect_str)
        # command.run(timeout=120)
        logging.info("Output file: %s\n" % out_put_file_loc)
        prediction_output_file_handle.close()







def RF_10flod_cross_validation(tmplist,pathdict,set_,logging):
    logging.info('10-fold cross validatation is running')
    data = pd.read_csv('/scratch2/zeng/homotypic_data/data/RandomForest/PsEnCo/TrainData2',delimiter="\s",engine='python')
    del data["Residue_id"]
    del data["Residue_name"]
    #print(data.as_matrix(data.columns))
    features=data.columns[0:28]
    X=data[features]
    y=data["Bind"]
    #n_samples, n_features = X.shape
    #random_state = np.random.RandomState(0)
    #X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]
    #print(X.iloc[[20,22]])
    cv = StratifiedKFold(y, n_folds=6)
    forest = RandomForestClassifier(n_estimators=100)
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []
    for i, (train, test) in enumerate(cv):
        probas_ = forest.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, 'k--',
             label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

def RF_variable_importance_calculate(tmplist,pathdic,set_,logging):
    data = pd.read_csv('/scratch2/zeng/homotypic_data/data/RandomForest/PsEnCo/TrainData2',delimiter="\s")
    del data["Residue_id"]
    del data["Residue_name"]
    #print(data.as_matrix(data.columns))
    features=data.columns[0:28]
    X=data[features]
    y=data["Bind"]
    forest = RandomForestClassifier(n_estimators=100)
    forest.fit(X, y)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    print("Feature ranking:")

    for f in range(X.shape[1]):
        print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
    plt.figure()
    plt.title("Feature importances")
    plt.bar(range(X.shape[1]), importances[indices],
           color="r", yerr=std[indices], align="center")
    #plt.xticks(range(x.shape[1]), indices)
    plt.xticks(range(X.shape[1]), indices)
    plt.xticks(range(X.shape[1]), indices)
    plt.xlim([-1, X.shape[1]])
    plt.show()