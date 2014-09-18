import numpy as np
import os.path
from sklearn import svm
import joblib
import pandas
#TODO replace with sklearns.externals

annotation_command = "python ~/code/cellulose_degraders/annotation/write_annot_array_script.py -d %(ORF_file)s -t %(annotation_out)s -s %(modes)s -f %(commands_out)s -h"

"""predict new samples"""
def normalize(test_data, norm_f):
    #scale row sums to 1 
    rs = np.sum(test_data, 1)
    test_data_rs = (test_data.T / rs).T
    #scale columns by normalization factors and return
    return test_data_rs/norm_f

def predict(pt, model_dir, test_data):
    """predict the class label for test_data based on the given model"""
    pass
    
def majority_predict(pt, model_dir, test_data):
    """predict the class label based on a committee vote of the models in models""" 
    #check if the model was based on binary input and load normalization factor or binarize elsewise
    nf = "%s%s_normf.dat"%(model_dir, pt)
    if os.path.exists(nf):
        norm_f = np.loadtxt(nf) 
        test_data_n = normalize(test_data, norm_f)
        #TODO test normalization
    else: 
        test_data_n = (test_data > 0).astype('int')
    predictors = joblib.load("%spickled/%s_predictors.pkl"%(model_dir, pt))
    #build prediction matrix
    preds = np.zeros((test_data.shape[0], len(predictors)))
    for i in range(len(predictors)):
        preds[:, i] = predictors[i].predict(test_data_n) 
        print predictors[i].decision_function(test_data_n)
    return preds

def annotate_and_predict(pt, pfam_f, model_dir, test_data_f, ORF_f, modes, out_dir):
    command = annotation_command % {"ORF_file":ORF_f, "annotation_out":"%s/annotation"%(out_dir), "commands_out":"%s/commands"%(out_dir), "modes":modes}
    os.system(command)
    #predict
    m = pandas.read_csv(test_data_f, sep="\t",)
    pfams = pandas.read_csv(pfam_f, header=None) 
    m_red = m.loc[:,pfams.iloc[:, 0]]
    return majority_predict(pt, model_dir, m_red)


if __name__ == "__main__":
    #test scenario with cow rumen bins
    test_data_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/annotation/pfam_27.0_hmmer_3.0_summary.txt"
    pfam_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfam_names.txt"
    #reduce to features in data set
    model_dir = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/ll_models/bin_phypat/"
    ORF_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/crb_ORFs.txt"
    out_dir = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/"
    modes = "2"
    print annotate_and_predict(8479, pfam_f, model_dir , test_data_f, ORF_f, modes, out_dir)
    






































