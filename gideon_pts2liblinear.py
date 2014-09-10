"""convert input matrix of Pfams and phenotypes into liblinear input format"""



#params:
#indices for target variables (phenotypes)
#out
#import liblinear as ll
import liblinearutil as lu
#from liblinear import *
import numpy as np
import os.path


def normalize(array):
    """normalize a feature matrix
    1) by making the rows sum to 1
    2) by making the max value for each feature equal 1"""
    #by row
    rs = np.sum(array,1)
    #features with no annotation are discarded
    rs[rs==0] = 1000
    #print rs
    array_rs = (array.T/rs).T
    #print np.sum(array_rs,1)
    #by column
    cs = np.max(array_rs,0)
    cs[cs==0] = 1000
    array_rs_cs = array_rs/cs
    #print np.max(array_rs_cs,0)
    return array_rs_cs

def cv_and_fs(data_f = None, bin=True,rec_based=False, model_out="stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g3_l1/", cv_outer=10):
    #read in data if not in anchestral character state prediction
    if data_f is not None:
        a = np.genfromtxt(data_f, missing_values=["?"], dtype=int)

    f = open("%s/cv_acc.txt"%model_out, "w")
    for pt in range(8476,8568):
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
            x = normalize(x.astype('double'))
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
        params_raw = '-s 5 -c 0.001 -w+1 %s -w-1 %s'%(sum(y==-1),sum(y==+1))
        print "pt_out", pt, "with", sum(y==+1), "positive samples and", sum(y==-1), "negative samples"
        params = lu.parameter(params_raw)
        #do n fold cross validation
        #number of elements in each invidual fold
        print "y length", len(y)
        fold_l = len(y)/cv_outer
        print "fold_l", fold_l
        #number of folds with one more element
        fold_l_e = len(y)%cv_outer
        print "fold_l_e", fold_l_e
        #list of folds lists
        folds = []
        for i in range(fold_l_e):
            folds.append(range(i*(fold_l+1), i*(fold_l+1)+fold_l+1))
        n_e = fold_l_e * (fold_l+1)
        for i in range(cv_outer - fold_l_e):
            folds.append(range(n_e + fold_l*i, n_e + fold_l*i + fold_l))

        all_preds = []
        #TODO account for the case when folds exceeds number of samples
        #TODO think about adjusting the class weights in each run (the difference should only be marginal
        for i in range(len(folds)):
            #train on all folds execept ith fold
            #negative index array
            neg_fold = np.ones(len(y),dtype='bool')
            for k in folds[i]:
                neg_fold[k] = False
            x_bin_co = x[neg_fold,]
            y_bin_co = y[neg_fold]
            p = lu.problem(list(y_bin_co),[list(j) for j in x_bin_co])
            m = lu.train(p, params)
            #test on left out fold
            preds = lu.predict(list(y[folds[i]]), [list(j) for j in x[folds[i]]],m)
            #append to all predictions
            all_preds += preds[0]
        #compare predicted labels with true label and compute balanced accuracy (use np array for that)
        #print all_preds
        all_preds_a = np.array(all_preds)
        #accuracy +1 class
        pos_acc = sum(y[y == 1] == all_preds_a[y==1])/float(sum(y==+1))
        #accuracy -1 class
        neg_acc = sum(y[y == -1] == all_preds_a[y==-1])/float(sum(y==-1))
        #balanced accuracy
        b_acc = float(pos_acc + neg_acc)/2
        f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, b_acc))
        f.flush()
        print "pos_acc", pos_acc
        print "neg_acc", neg_acc
        print "balanced acc", b_acc
        #full model with all data points:
        #params_cv = lu.parameter('-s 5 -c 0.1 -v 10')
        x_bin_l = [list(i) for i in x]
        problem = lu.problem(list(y),x_bin_l)
        m = lu.train(problem,params)
        #pred = lu.predict(list(y),x_bin_l,m)
        #lu.save_model("stol_2_NCBI20140115_candidatus_sample30/ll_models/model_%s.txt"%pt,m)
        lu.save_model("%smodel_%s.txt"%(model_out,pt_out),m)
    f.close()
phy_p="stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
#cv_and_fs(data_f=phy_p, bin=True, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/bin/", cv_outer=10)
#cv_and_fs(data_f=phy_p, bin=False, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/norm/", cv_outer=10)
cv_and_fs( bin=True, rec_based=True, model_out = "stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g2_l1/", cv_outer=10)
