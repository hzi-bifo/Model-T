import nested_cv as ncv
import numpy as np
import os.path
import sys
import shutil
import getopt

#c_params = [0.1,1]
c_params = [0.03, 0.07, 0.1, 0.3, 0.7, 1]
#c_params = [1,5,10,50,100]

def write_miscl(a, miscl_plus,  model_out, pt_out):
    gideon_f = open("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gideon_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_stol_20130904.sp_uq.taxid.txt", "r")
    ls = gideon_f.readlines()
    id2sp = {}
    f = open("%s/%s_miscl.txt"%(model_out,pt_out), 'w')
    for i in range(len(ls)):
        elems = ls[i].strip().split("\t")
        id2sp[i] = (elems[0],elems[1])
    for m in miscl_plus.T:
        print m
        try:
            f.write("%s\t%s\t%s\t%s\n"%(id2sp[m[0]][0],m[1], m[2], id2sp[m[0]][1]))
        except KeyError:
            break
    f.close()



def cv_and_fs(bin, model_out, gt_start, gt_end, pt_start, pt_end, data_f = None, cv_outer=10, cv_inner=None, n_jobs = 1):
    #create output directory if not already existing, append input parameter specific suffices
    #model_out = "%s/%s_%s"%(model_out, "bin" if bin else "counts", "recbased" if os.path.isdir(data_f) else "phypat")
    #if not os.path.exists(model_out):
    #    os.mkdir(model_out)
    #else:
    #    sys.stderr.write("warning output directory %s already exists\n"%model_out)
    #    sys.stderr.write("determine what you want to do\n")
    #    while True:
    #        c = raw_input("press 0 for abort and 1 for overwrite (delete existing directory)\n")
    #        if c == "0": 
    #            sys.exit(1)
    #        elif c == "1": 
    #            shutil.rmtree(model_out, ignore_errors=True)
    #            break 
    #write values of the C parameter to disk
    with open("%s/c_params.txt" % model_out, 'w') as out_f:
        out_f.write("C_params\n")
        for c in c_params:
            out_f.write("%s\n"%c)

    #read in data if not in anchestral character state prediction
    if os.path.isfile(data_f):
        a = np.genfromtxt(data_f, missing_values=["?"], dtype=int)

    f = open("%s/cv_acc.txt"%model_out, "w")
    #create a directory for storing the pickled models
    os.mkdir("%s/pickled"%model_out)
    for pt in range(pt_start, pt_end+1):
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
            pt = a.shape[1] - 1 
        #discard samples which are missing for that phenotype and keep mapping of samples to orginal matrix a in x2a
        x2a =  np.ravel(np.logical_not(a[:,pt]==-1)).nonzero()
        x=a[np.logical_not(a[:,pt]==-1),gt_start:gt_end+1]
        if bin:
            x = (x>0).astype('int')
        else:
            x, nf = ncv.normalize(x.astype('double'))
            #save the normalization factors to disk for later testing
            np.savetxt(fname="%s/%s_normf.dat"%(model_out, pt_out), X=nf)
             
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
            all_preds = np.array(ncv.outer_cv(x,y, params, c_params, cv_outer, cv_inner, n_jobs=n_jobs))
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
            write_miscl(a, miscl_plus, model_out, pt_out)
            f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc))
            f.flush()
        all_preds = ncv.outer_cv(x, y, params, c_params, cv_outer, cv_inner=None, n_jobs=n_jobs)
        ncv.majority_feat_sel(x, y, all_preds, params, c_params, 5, model_out, pt_out)
    f.close()
if __name__=="__main__":
    #only testing
    #phy_p="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
    #pars_p="/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
    #cv_and_fs( bin=True,  model_out = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/ll_models/", data_f = phy_p, cv_outer=10, cv_inner=None)
    #cv_and_fs(data_f=phy_p, bin=False, rec_based=False, model_out = "stol_2_NCBI20140115_candidatus_sample30/ll_models/norm/", cv_outer=10)
    #cv_and_fs( bin=True, rec_based=True, model_out = "stol_2_NCBI20140115_candidatus_sample30/parsimony/models_g2_l1/", cv_inner = 5, cv_outer=5)
    if len(sys.argv) == 1:
        print """USAGE: python %s
-d <input dir> with one matrix for each phenotype or in case of phyletic pattern classification one matrix with all the phenotypes 
-v <cross validation folds> the number of folds used for inner cross validation  
-i <inner cross validation folds> the number of folds used for inner cross validation if this option is given nested cross validation will be performed otherwise only feature selection routines will launched  
-o <out dir> for the models, selected features etc. 
-c <use counts> use raw data with normalization and don't binarize 
-j <number of jobs> that shall be used
-g <range of genotypes> e.g. 1-8400
-p <range of phenotypes> to consider e.g 8550-8560
        """ % (sys.argv[0])
        sys.exit(2)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "d:v:i:cj:g:p:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    in_data = None
    out = None
    cv_outer = None
    cv_inner = None
    binary = True 
    g1 = g2 = None
    pt1 = pt2 = None
    n_jobs = 1 
    for o, a in optlist:
        if o == "-d":
            in_data = a
        if o == "-v":
            cv_outer = int(a)
        if o == "-i":
            cv_inner = int(a)
        if o == "-c":
            binary = False
        if o == "-j":
            n_jobs = int(a)
        if o == "-g":
            g1, g2 = [int(i) for i in a.split("-")]
        if o == "-p":
            pt1, pt2 = [int(i) for i in a.split("-")]
        if o == "-o":
            out = a 
            if os.path.exists(out):
                sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
                sys.exit(1)
            else:
                os.mkdir(out)
    cv_and_fs(bin=binary, data_f = in_data, gt_start = g1, gt_end = g2, pt_start = pt1, pt_end = pt2, cv_inner = cv_inner, cv_outer = cv_outer, model_out = out, n_jobs = n_jobs)
