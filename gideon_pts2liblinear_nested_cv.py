import nested_cv as ncv
import numpy as np
import os.path
import liblinearutil as lu
import sys

c_params = [0.0001, 0.001]

def cv_and_fs(data_f = None, bin=True,rec_based=False, model_out="stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g3_l1/", cv_outer=10, cv_inner=None):
    #read in data if not in anchestral character state prediction
    if data_f is not None:
        a = np.genfromtxt(data_f, missing_values=["?"], dtype=int)

    #f = open("%s/cv_acc.txt"%model_out, "w")
    #for pt in range(8476,8568):
    for pt in range(8486,8487):
        pt_out = pt
        print pt
        #feature matrix with missing data removed
        if rec_based:
            if not  os.path.exists("stol_2_NCBI20140115_candidatus_sample30/parsimony/input_g2_l1/pt%s.dat"%pt):
                print "skipping pt", pt, "no reconstruction matrix found, possibly due to no reconstruction events for this phenotype"
                continue
            a = np.genfromtxt("stol_2_NCBI20140115_candidatus_sample30/parsimony/input_g2_l1/pt%s.dat"%pt)
            #treat gains and losses the same
            a[a==-1]=1
            #discard the row names
            a = a[:,1:(a.shape[1])]
            #phenotype index
            pt = 8476
        x=a[np.logical_not(a[:,pt]==-1),0:8476]
        if bin:
            x = (x>0).astype('int')
        else:
            x = ncv.normalize(x.astype('double'))
            #print x
        #experiment: add inverse features
        #x_bin_inv = x_bin.copy()
        #plus =x_bin_inv==1
        #zero =x_bin_inv==0
        #x_bin_inv[plus] = 0
        #x_bin_inv[zero] = 1
        #print x_bin.shape, x_bin_inv.shape
        #x_bin = np.concatenate((x_bin,x_bin_inv),axis=1)
        #end experiment
        #target vector with missing data removed
        y=a[np.logical_not(a[:,pt]==-1),pt]
        y[y==0] = -1
        #experiment: reduce the number of negative examples to approximately the number of positive samples
        #neg_class = (y==-1)
        #pos_class = (y==1)
        #x_bin = np.concatenate((x_bin[pos_class,:] , x_bin[neg_class[0:30],:]), axis=0)
        #y = np.concatenate((y[pos_class], y[neg_class[0:30]]))
        #end experiment
        #skip phenotypes with too little samples
        if len(y) <= 30:
            print "skipping pt %s with only %s observations"%(pt_out,len(y))
            continue
        if sum(y==+1) < 10:
            print "skipping pt %s with only %s + observations"%(pt_out,sum(y==+1))
            continue
        if sum(y==-1) < 10:
            print "skipping pt %s with only %s - observations"%(pt_out,sum(y==-1))
            continue
        params_raw = '-s 5 -w+1 %s -w-1 %s'%(sum(y==-1),sum(y==+1))
        print "pt_out", pt, "with", sum(y==+1), "positive samples and", sum(y==-1), "negative samples"
        params = lu.parameter(params_raw)
        #get predictions for all values of the parameter C
        all_preds = ncv.outer_cv(x,y, params_raw, c_params, cv_outer, cv_inner, n_jobs=15)
        #compare predicted labels with true label and compute balanced accuracy (use np array for that)
        #print all_preds
        if not cv_inner is None:
            all_preds_a = np.array(all_preds)
            #accuracy +1 class
            pos_acc = ncv.recall_pos(y, all_preds_a)
            #accuracy -1 class
            neg_acc = ncv.recall_neg(y, all_preds_a)
            #balanced accuracy
            bacc = ncv.bacc(pos_acc, neg_acc)
            print "pos_acc", pos_acc
            print "neg_acc", neg_acc
            print "balanced acc", bacc
            #TODO write misclassified samples to disk
            f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc))
            f.flush()
        ncv.majority_feat_sel(x, y, all_preds, params_raw, c_params, 2, model_out, pt_out)
if __name__=="__main__":
    phy_p="stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
    cv_and_fs(data_f=phy_p, bin=True, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/bin_exp/", cv_outer=5, cv_inner=None)
    #cv_and_fs(data_f=phy_p, bin=False, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/norm/", cv_outer=10)
    #cv_and_fs( bin=True, rec_based=True, model_out = "stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g2_l1/", cv_inner = 5, cv_outer=5)
