import nested_cv as ncv
import numpy as np
import os.path
import sys
import shutil
import getopt
import pandas as ps
import random
random.seed(0)

#c_params = [0.05, 0.05]
#c_params = [0.03, 0.07, 0.1, 0.3, 0.7, 1, 3, 7, 10]
#c_params = [0.02, 0.03, 0.05, 0.07, 0.1, 0.3, 0.5, 0.7, 1]
c_params = [0.001, 0.002, 0.005, 0.007,  0.01,  0.02,  0.05, 0.07, 0.1, 0.2, 0.5, 0.7, 1]
MIN_POS = 10
MIN_NEG = 10
MIN_SAMPLES = 20
SP2TXID = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gideon_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_stol_20130904.sp_uq.taxid.txt"
params = {'loss':'l2', 'tol':0.000001, 'penalty':'l1', 'dual':False, 'fit_intercept':True, 'intercept_scaling':1, 'class_weight':'auto', 'random_state':1}
#c_params = [1,5,10,50,100]

def write_miscl(miscl_plus,  model_out, pt_out):
    gideon_f = open(, "r")
    ls = gideon_f.readlines(SP2TXID)
    id2sp = {}
    f = open("%s/%s_miscl.txt"%(model_out,pt_out), 'w')
    for i in range(len(ls)):
        elems = ls[i].strip().split("\t")
        id2sp[int(elems[1])] = elems[0]
    for i in range(miscl_plus.shape[0]):
        f.write("%s\t%s\t%s\t%s\n"%(miscl_plus.index[i], miscl_plus.iloc[i,0], miscl_plus.iloc[i,1], id2sp[miscl_plus.index[i]]))
    f.close()



def cv_and_fs(model_out, gt_start, gt_end, pt_start, pt_end, phypat_f, rec_dir, likelihood_params, parsimony_params, is_phypat_and_rec,  cv_outer=10, cv_inner=None, n_jobs = 1, perc_samples = 1.0, perc_feats = 1.0, inverse_feats = False, do_normalization = False):
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
    with open("%s/configuration.txt" % model_out, 'w') as out_f:
        out_f.write("""MIN_POS\t%s\nMIN_NEG\t%s\nMIN_SAMPLES\t%s\ncv_inner\t%s\ncv_outer\t%s\ntol\t%s\nfit_intercept\t%s\nclass_weight\t%s\nrandom_state\t%s\n"""%(MIN_POS, MIN_NEG, MIN_SAMPLES, cv_outer, cv_inner, params['tol'], params['fit_intercept'], params['class_weight'], params['random_state']))
    #read in phyletic patterns
    p = ps.read_csv(phypat_f, sep = "\t", na_values = ["?"], index_col = 0, header = None)
    f = open("%s/cv_acc.txt"%model_out, "w")
    #create a directory for storing the pickled models
    os.mkdir("%s/pickled"%model_out)
    #iterate over all phenotypes
    for pt in range(pt_start, pt_end+1):
        pt_out = pt
        x_p = p[p.iloc[:,pt].notnull()].iloc[:, gt_start:gt_end+1]
        y_p = p[p.iloc[:,pt].notnull()].iloc[:, pt]
        #shuffle phyletic pattern samples to avoid biases in the cross validation
        rows = x_p.index.values.copy()
        random.shuffle(rows)
        x_p = x_p.loc[rows, ]
        y_p = y_p.loc[rows, ]
        if not rec_dir is None:
            if not os.path.exists("%s/pt%s.dat"%(rec_dir,pt)):
                print "skipping pt", pt, "no reconstruction matrix found, possibly due to no reconstruction events for this phenotype"
                continue
            a = ps.read_csv("%s/pt%s.dat"%(rec_dir,pt), index_col = 0, sep = "\t", header = None)
            #treat gains and losses the same
            #only relevant for the parsimony case
            if not parsimony_params is None :
                raise Exception("not yet implemented")
            #a[a==-1]=1
            #discard the row names
            #TODO change to panda dataframe to keep the row names
            a = a.iloc[:,0:(a.shape[1])]
            #phenotype index in the reconstruction matrix
            print pt
            pt = a.shape[1] - 1 
            x = a.iloc[:, gt_start:gt_end+1]
            y = a.iloc[:, pt]
        else:
            #we are in pure phyletic pattern mode
            x = x_p
            y = y_p
        if not do_normalization:
            x_p = (x_p > 0).astype('int')
            if likelihood_params is None:
                x = (x > 0).astype('int')
        #else:
            #x_p, nf = ncv.normalize(x_p.astype('double'))
            #np.savetxt(fname="%s/%s_normf_xp.dat"%(model_out, pt_out), X=nf)
            #x, nf = ncv.normalize(x.astype('double'))
            #np.savetxt(fname="%s/%s_normf_x.dat"%(model_out, pt_out), X=nf)
            #if is_phypat_and_rec:
            #    x_join, scaler = ncv.normalize(ps.concat([x.astype('double'), x_p.astype('double')], axis = 0))
            #    x = x_join.iloc[0: x.shape[0], ]
            #    x_p = x_join.iloc[x.shape[0]:x_join.shape[0], ]
            #    print x_p
            #    print x
            #    sys.exit(0)
            #else:
            #    x, scaler = ncv.normalize(x.astype('double'))
            #    x_p = ps.DataFrame(data = scaler.transform(x_p.astype('double')), index = x_p.index, columns = x_p.columns)
            #    print x_p
                
            #save the normalization factors to disk for later testing
            #TODO what about reconstruction case??
            #TODO if this goes alright
             
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
        y[(y==0)] = -1
        #if in likelihood reconstruction case set entries of the target vector with likely gain and loss to 1
        #TODO what about the other entries in x?
        if not rec_dir is None and not likelihood_params is None:
            y[(y==2)] = 1
        #experiment: reduce the number of negative examples to approximately the number of positive samples
        #neg_class = (y==-1)
        #pos_class = (y==1)
        #x_bin = np.concatenate((x_bin[pos_class,:] , x_bin[neg_class[0:30],:]), axis=0)
        #y = np.concatenate((y[pos_class], y[neg_class[0:30]]))
        #end experiment
        #check if we are in likelihood mode and y consists of target probabilities instead of discrete values
            #TODO sample restrictions don't hold anymore when doubling the data set
            #TODO all class 0 and class 1 samples are contingiuos in the input, shuffle or make 0,1,0,1
        #skip phenotypes with too little samples
        #in case of joint phyletic pattern and reconstruction-based classification, sum up observations
        print y
        y_join = y.copy()
        if is_phypat_and_rec:
            y_join =  y_join.append(y_p)
        if len(y_join) < MIN_SAMPLES:
            print "skipping pt %s with only %s observations"%(pt_out,len(y_join))
            continue
        if sum(y_join > 0) < MIN_POS:
            print "skipping pt %s with only %s + observations"%(pt_out,sum(y_join>0))
            continue
        if sum(y_join==-1) < MIN_NEG:
            print "skipping pt %s with only %s - observations"%(pt_out,sum(y_join==-1))
            continue
        print "pt_out", pt, "with", sum(y_join>0), "positive samples and", sum(y_join==-1), "negative samples"
        #check if we want to do nested cross validation accuracy estimation
        is_rec_based = False
        if not rec_dir is None:
            is_rec_based = True
        if not cv_inner is None:
            try:
                all_preds = ps.Series(np.array(ncv.outer_cv(x,y, params, c_params, cv_outer, cv_inner, n_jobs, is_rec_based, x_p, y_p, model_out, likelihood_params, parsimony_params, is_phypat_and_rec, perc_samples, perc_feats,  inverse_feats, do_normalization)))
            except ValueError as e:
                import traceback
                print traceback.print_exc(e)
                print e
                print "pt %s has too few samples in one of the classes possibly, due to a high threshold in the likelihood reconstruction; try to lower the number of cross validation fold and run again"
                continue
            #accuracy +1 class
            all_preds.index = y_p.index
            y_p_t = y_p.copy()
            y_p_t[y_p_t == 0] = -1
            #print all_preds, y_p_t, "y_p and all preds nested cv"
            pos_acc = ncv.recall_pos(y_p_t, all_preds)
            print "pos_acc", pos_acc
            #accuracy -1 class
            neg_acc = ncv.recall_neg(y_p_t, all_preds)
            print "neg_acc", neg_acc
            #balanced accuracy
            bacc = ncv.bacc(pos_acc, neg_acc)
            print "balanced acc", bacc
            #TODO get misclassified reconstructions samples
            miscl = y_p_t.index[(all_preds != y_p_t)]
            #bind actual labels and predictions
            #print miscl, y_p_t.loc[miscl], all_preds.loc[miscl]
            miscl_plus = ps.concat([y_p_t.loc[miscl], all_preds.loc[miscl]], axis = 1)
            #print miscl_plus
            write_miscl(miscl_plus, model_out, pt_out)
            f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc))
            f.flush()
        all_preds = ps.DataFrame(ncv.outer_cv(x,y, params, c_params, cv_outer, cv_inner = None, n_jobs = n_jobs, is_rec_based = is_rec_based, x_p = x_p, y_p = y_p, model_out = model_out, likelihood_params = likelihood_params, parsimony_params = parsimony_params, is_phypat_and_rec = is_phypat_and_rec, perc_samples = perc_samples, perc_feats = perc_feats, inverse_feats = inverse_feats, do_normalization = do_normalization))
        ncv.majority_feat_sel(x, y, x_p, y_p, all_preds, params, c_params, 5, model_out, pt_out, is_phypat_and_rec, perc_samples = perc_samples, perc_feats = perc_feats, likelihood_params = likelihood_params, inverse_feats = inverse_feats, do_normalization = do_normalization)
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
-y <data_f> phyletic patterns, i.e. one matrix with all the phenotypes 
-d <input dir> with one matrix for each phenotype or in case of phyletic pattern classification one matrix with all the phenotypes 
-l <threshold:x,mode:<gain, loss, gain_loss>> only use if in reconstruction classification i.e. option -d is set as well, in that case -l means that we have parsimony reconstruction, this is followed by the options for the likelihood-based reconstruction e.g. -l threshold:0.5,mode:gain 
-a <gain_costs:x,mode:<ACCTRAN, DELTRAN, RANDOM>> only use if in reconstruction classification i.e. option -d is set as well, in that case -a means that we have parsimony reconstruction, this is followed by the options for the parsimony-based reconstruction e.g. -l threshold:0.5,mode:gain 
-b only use if in reconstruction classification i.e. option -d is set as well, in that case -b means that training is done on both the phyletic patterns and the parsimony / likelihood reconstruction 
-v <inner cross validation folds> the number of folds used for inner cross validation  
-i <inner cross validation folds> the number of folds used for inner cross validation if this option is given nested cross validation will be performed otherwise only feature selection routines will launched  
-o <out dir> for the models, selected features etc. 
-c <use counts> use raw data with normalization and don't binarize 
-j <number of jobs> that shall be used
-g <range of genotypes> e.g. 1-8400
-p <range of phenotypes> to consider e.g 8550-8560
-r <percentage of features> to to use for the individual classifiers. 1.0 corresponds to the vanilla SVM with all features included
-c <inverse feats> if option set, extend the feature space in the classification by the inverse features (1-X) to get rid of noise features
        """ % (sys.argv[0])
        sys.exit(2)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "y:d:l:a:bv:i:cj:g:p:o:s:r:e")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    phypat_f = None
    rec_dir = None
    likelihood_params = None
    parsimony_params = None
    is_phypat_and_rec = False
    out = None
    cv_outer = None
    cv_inner = None
    do_normalization = False 
    g1 = g2 = None
    pt1 = pt2 = None
    n_jobs = 1 
    perc_samples = 1.0
    perc_feats = 1.0
    inverse_feats = False 

    for o, a in optlist:
        if o == "-d":
            rec_dir = a
        if o == "-y":
            phypat_f = a
        if o == "-l":
            #parse likelihood option
            likelihood_params =  dict(i.split(":") for i in a.strip().split(","))
        if o == "-a":
            #parse likelihood option
            parsimony_params =  dict(i.split(":") for i in a.strip().split(","))
        if o == "-b":
            is_phypat_and_rec = True
        if o == "-v":
            cv_outer = int(a)
        if o == "-i":
            cv_inner = int(a)
        if o == "-c":
            do_normalization = True 
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
        if o == "-s":
            perc_samples = float(a)
        if o == "-r":
            perc_feats = float(a)
        if o == "-e":
            inverse_feats = True
    cv_and_fs(phypat_f = phypat_f, gt_start = g1, gt_end = g2, pt_start = pt1, pt_end = pt2, rec_dir = rec_dir, likelihood_params = likelihood_params, parsimony_params = parsimony_params, is_phypat_and_rec = is_phypat_and_rec, cv_inner = cv_inner, cv_outer = cv_outer, model_out = out, n_jobs = n_jobs, perc_samples = perc_samples, perc_feats = perc_feats, inverse_feats = inverse_feats, do_normalization = do_normalization)
