import nested_cv as ncv
import numpy as np
import os.path
import sys

#c_params = [0.1,1]
c_params = [0.03, 0.07, 0.1, 0.3, 0.7, 1]
#c_params = [1,5,10,50,100]

def write_miscl(miscl_plus,  model_out, pt_out):
    f = open("%s%s_miscl.txt"%(model_out,pt_out), 'w')
    gideon_f = open("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gideon_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_stol_20130904.sp_uq.taxid.txt", "r")
    id2sp = {}
    ls = gideon_f.readlines()
    for i in range(len(ls)):
        elems = ls[i].strip().split("\t")
        id2sp[i] = (elems[0],elems[1])
    for m in miscl_plus.T:
        print m
        f.write("%s\t%s\t%s\t%s\n"%(id2sp[m[0]][0],m[1], m[2], id2sp[m[0]][1]))
    f.close()



def cv_and_fs(data_f = None, bin=True, model_out="stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g3_l1/", cv_outer=10, cv_inner=None):
    #create output directory if not already existing, append input parameter specific suffices
    model_out = "%s_%s_%s"(model_out, "bin" if bin else "counts", "recbased" if os.path.isdir(data_f) else "phypat")
    if not os.path.exists(model_out):
        os.mkdir(model_out)
    else:
        sys.stderr("warning output directory %s already exists, exiting")
        sys.exit(1)
    #write values of the C parameter to disk
    with open("%s/c_params.txt" % model_out, 'w') as out_f:
        out_f.write("C_params\n")
        for c in c_params:
            out_f.write("%s\n"%c)

    #read in data if not in anchestral character state prediction
    if os.path.isfile(data_f):
        a = np.genfromtxt(data_f, missing_values=["?"], dtype=int)

    f = open("%s/cv_acc.txt"%model_out, "w")
    for pt in range(8476,8568):
    #for pt in range(8486,8487):
        pt_out = pt
        print pt
        #feature matrix with missing data removed
        if os.path.isdir(data_f):
            if not  os.path.exists("%s/pt%s.dat"%(data_f,pt)):
                print "skipping pt", pt, "no reconstruction matrix found, possibly due to no reconstruction events for this phenotype"
                continue
            a = np.genfromtxt("%s/pt%s.dat"%(data_f,pt))
            #treat gains and losses the same
            a[a==-1]=1
            #discard the row names
            #TODO change to panda dataframe to keep the row names
            a = a[:,1:(a.shape[1])]
            #phenotype index
            pt = 8476
        #discard samples which are missing for that phenotype and keep mapping of samples to orginal matrix a in x2a
        x2a =  np.ravel(np.logical_not(a[:,pt]==-1)).nonzero()
        x=a[np.logical_not(a[:,pt]==-1),0:8476]
        if bin:
            x = (x>0).astype('int')
        else:
            x = ncv.normalize(x.astype('double'))
            #print x
        #experiment: add inverse features
        #x_inv = x.copy()
        #plus =x_inv==1
        #zero =x_inv==0
        #x_inv[plus] = 0
        #x_inv[zero] = 1
        #print x.shape, x_inv.shape
        #x = np.concatenate((x,x_inv),axis=1)
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
        params = {'loss':'l2', 'tol':0.000001, 'penalty':'l1', 'dual':False, 'fit_intercept':True, 'intercept_scaling':1, 'class_weight':'auto', 'random_state':1}
        print "pt_out", pt, "with", sum(y==+1), "positive samples and", sum(y==-1), "negative samples"
        #check if we want to do nested cross validation accuracy estimation
        if not cv_inner is None:
            all_preds = np.array(ncv.outer_cv(x,y, params, c_params, cv_outer, cv_inner, n_jobs=15))
            #accuracy +1 class
            pos_acc = ncv.recall_pos(y, all_preds)
            #accuracy -1 class
            neg_acc = ncv.recall_neg(y, all_preds)
            #balanced accuracy
            bacc = ncv.bacc(pos_acc, neg_acc)
            print "pos_acc", pos_acc
            print "neg_acc", neg_acc
            print "balanced acc", bacc
            #TODO write misclassified samples to disk
            #get misclassified samples
            miscl = (all_preds != y).ravel().nonzero()[0]
            #bind actual labels and predictions
            miscl_plus = np.array([miscl, y[miscl], all_preds[miscl]])
            print miscl_plus
            write_miscl(miscl_plus, model_out, pt_out)
            f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc))
            f.flush()
        all_preds = ncv.outer_cv(x, y, params, c_params, cv_outer, cv_inner=None, n_jobs=15)
        ncv.majority_feat_sel(x, y, all_preds, params, c_params, 5, model_out, pt_out)
    f.close()
if __name__=="__main__":
    phy_p="stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
    pars_p="stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
    cv_and_fs(data_f = phy_p, bin=True,  model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/bin_exp_sklearn/", cv_outer=10, cv_inner=10)
    #cv_and_fs(data_f=phy_p, bin=False, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/norm/", cv_outer=10)
    #cv_and_fs( bin=True, rec_based=True, model_out = "stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g2_l1/", cv_inner = 5, cv_outer=5)
