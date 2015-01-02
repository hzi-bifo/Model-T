import numpy as np
import os.path
from sklearn import svm
import joblib
#TODO replace with sklearns.externals
import pandas as ps
import sys

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

def filter_pred(scores):
    #print scores
    if (scores > 0).all():
        return scores.mean()
    #elif (scores < 0).all(): 
    #    return scores.mean()
    else:
        return ps.np.NaN

def aggregate(pred_df):
    """restrict to positive restrictions that are unique across all c params"""
    maj_pred_df = ps.DataFrame(np.zeros(shape = (pred_df.shape[0], pred_df.shape[1] / 5)))
    for i in range(0, pred_df.shape[1], 5):
        print i
        print "filtered preds", pred_df.iloc[:, i : i + 5].apply(filter_pred, axis = 1)
        maj_pred_df.iloc[:,i / 5] = pred_df.iloc[:, i: i + 5].apply(filter_pred, axis = 1)
        maj_pred_df.columns.values[i / 5] = pred_df.columns.values[i].split("_")[0] 
    return maj_pred_df

    
def majority_predict(pt, model_dir, test_data):
    """predict the class label based on a committee vote of the models in models""" 
    #check if the model was based on binary input and load normalization factor or binarize elsewise
    nf = "%s%s_normf.dat"%(model_dir, pt)
    if os.path.exists(nf):
        norm_f = np.loadtxt(nf) 
        test_data_n = normalize(test_data, norm_f)
        #TODO replace with the normalization used in the gideon2pts script 
    else: 
        #binarize
        test_data_n = (test_data > 0).astype('int')
    if os.path.exists("%spickled/%s_predictors.pkl"%(model_dir, pt)):
        predictors = joblib.load("%spickled/%s_predictors.pkl"%(model_dir, pt))
    else: 
        return ps.DataFrame()
    #build prediction matrix
    preds = np.zeros((test_data.shape[0], len(predictors)))
    for i in range(len(predictors)):
        preds[:, i] = predictors[i].decision_function(test_data_n) 
    pred_df = ps.DataFrame(preds, index = test_data.index)
    #print predictors[0].get_params()['C']
    pred_df.columns = ["%s_%s" %(pt, predictors[i].get_params()['C']) for i in range(len(predictors))]
    return pred_df

def annotate_and_predict(pt_range, pfam_f, model_dir, test_data_f, ORF_f, modes, out_dir):
    pred_df = ps.DataFrame()
    command = annotation_command % {"ORF_file":ORF_f, "annotation_out":"%s/annotation"%(out_dir), "commands_out":"%s/commands"%(out_dir), "modes":modes}
    print command
    #print os.system(command)
    for pt in range(pt_range[0], pt_range[1] + 1):
        #predict
        m = ps.read_csv(test_data_f, sep="\t",)
        pfams = ps.read_csv(pfam_f, header=None) 
        m_red = m.loc[:,pfams.iloc[:, 0]]
        preds = majority_predict(pt, model_dir, m_red)
        pred_df = ps.concat([pred_df, preds], axis = 1)
    pred_df.to_csv("%s/predictions.csv"%out_dir, sep = "\t", float_format='%.3f')
    #aggregate predictions
    aggregate(pred_df).to_csv("%s/predictions_aggr.csv"%out_dir, sep = "\t", float_format='%.3f')
    return pred_df


if __name__ == "__main__":
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-p <range of phenotypes> specify which phenotypes shall be predicted e.g. 8487-8500
-s <output file> for the annotation summary
-m <model dir> directory with the phenotype predictors
-t <test data input file>
-o <ORF file> with one fasta file per sample per line
-a <mode> for the annotation ususally 2 for pfam
-d <out dir> for the annotation
        """ % (sys.argv[0])
        sys.exit(2)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "p:s:m:t:o:a:d:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    pt1 = pt2 = sum_f = model_dir = test_data_f = ORF_f = modes = out_dir = None
    for o, a in optlist:
        if o == "-p":
            pt1, pt2 = [int(i) for i in a.split("-")]
        if o == "-s":
            sum_f = a
        if o == "-m":
            model_dir = a
        if o == "-t":
            test_data_f = a
        if o == "-o":
            ORF_f = a
        if o == "-a":
            modes = a
        if o == "-d":
            out_dir = a
    annotate_and_predict((pt1, pt2), sum_f, model_dir, test_data_f, ORF_f, modes, out_dir) 
    #test scenario with cow rumen bins
    #test_data_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/annotation/pfam_27.0_hmmer_3.0_summary.txt"
    #pfam_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfam_names.txt"
    #reduce to features in data set
    #model_dir = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/ll_models/bin_phypat/"
    #ORF_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/crb_ORFs.txt"
    #out_dir = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/cow_rumen_bins/"
    #modes = "2"
    #print annotate_and_predict(8479, pfam_f, model_dir , test_data_f, ORF_f, modes, out_dir)
    






































