import nested_cv as ncv
import numpy as np
import os.path
import sys
import shutil
import getopt
import pandas as ps
import random
random.seed(0)
import json

class pt_classification:

    def write_miscl(self, model_out, pt_out, miscl_plus):
        """map the misclassified sample ids to their scientific names and taxonomic ids"""
        gideon_f = open(self.config['sp2taxid'], "r")
        ls = gideon_f.readlines()
        id2sp = {}
        f = open("%s/%s_miscl.txt"%(model_out, pt_out), 'w')
        for i in range(len(ls)):
            elems = ls[i].strip().split("\t")
            id2sp[int(elems[1])] = elems[0]
        for i in range(miscl_plus.shape[0]):
            f.write("%s\t%s\t%s\t%s\n"%(miscl_plus.index[i], miscl_plus.iloc[i,0], miscl_plus.iloc[i,1], id2sp[miscl_plus.index[i]]))
        f.close()
    
    def setup_folds_synthetic(self, x, folds):
        #short folds length
        f = x / folds
        #number of extended folds
        f_e = x % folds
        return (folds - f_e, x - f, f_e, ((x - (f + 1))))
    

    
    def __init__(self, config_f, model_out, gt_start, gt_end, pt_start, pt_end, phypat_f, rec_dir, likelihood_params, parsimony_params, is_phypat_and_rec,  cv_outer=10, cv_inner=None, n_jobs = 1, perc_samples = 1.0, perc_feats = 1.0, inverse_feats = False, do_normalization = False, resume = False):
        """main routine to prepare data for classification and feature selection"""
        self.model_out = model_out
        self.resume = resume
        #read in configuration file
        c = [] 
        with open(config_f, 'r') as f:
            for l in f.readlines():
                #strip python style comments
                if not l.strip().startswith("#"):
                    c.append(l.strip())
        #parse pure json
        self.config = json.loads("\n".join(c))
        #random.seed(self.config["seed"])
        #check if we want to do nested cross validation accuracy estimation
        is_rec_based = False
        if not rec_dir is None:
            is_rec_based = True
        self.ncv = ncv.nested_cv(likelihood_params, parsimony_params, do_normalization, is_rec_based, is_phypat_and_rec, n_jobs, inverse_feats, self.config, perc_feats, perc_samples, model_out, cv_outer, resume)
        if not os.path.exists(model_out):
            os.mkdir(model_out)
        #write config to disk
        with open("%s/config.json" % self.model_out, 'w') as out_f:
            json.dump(self.config, out_f, indent=4, separators=(',', ': '))
        #read in phyletic patterns
        print phypat_f
        p = ps.read_csv(phypat_f, sep = "\t", na_values = ["?"], index_col = 0, header = None)
        #check if this is a new run or if the current run shall be recomputed
        if self.resume:
            f = open("%s/cv_acc.txt"%self.model_out, "a")
        else:
            f = open("%s/cv_acc.txt"%self.model_out, "w")
        #create a directory for storing the pickled models
        if not os.path.exists("%s/pickled"%self.model_out):
            os.mkdir("%s/pickled"%self.model_out)
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
                #binarize
                x_p = (x_p > 0).astype('int')
                if likelihood_params is None:
                    x = (x > 0).astype('int')
    
            #target vector with missing data removed
            y[(y==0)] = -1
            #if in likelihood reconstruction case set entries of the target vector with likely gain and loss to 1
            #TODO what about the other entries in x?
            if not rec_dir is None and not likelihood_params is None:
                y[(y==2)] = 1
            #check if we are in likelihood mode and y consists of target probabilities instead of discrete values
            #TODO sample restrictions don't hold anymore when doubling the data set
            #TODO all class 0 and class 1 samples are contingiuos in the input, shuffle or make 0,1,0,1
            #skip phenotypes with too little samples
            #in case of joint phyletic pattern and reconstruction-based classification, sum up observations
            y_join = y.copy()
            if is_phypat_and_rec:
                y_join =  y_join.append(y_p)
            if len(y_join) < self.config['min_samples']:
                print "skipping pt %s with only %s observations"%(pt_out,len(y_join))
                continue
            if sum(y_join > 0) < self.config['min_pos']:
                print "skipping pt %s with only %s + observations"%(pt_out,sum(y_join>0))
                continue
            if sum(y_join==-1) < self.config['min_neg']:
                print "skipping pt %s with only %s - observations"%(pt_out,sum(y_join==-1))
                continue
            print "pt_out", pt, "with", sum(y_join>0), "positive samples and", sum(y_join==-1), "negative samples"
            if not cv_inner is None:
                try:
                    all_preds = ps.Series(np.array(self.ncv.outer_cv(x,y, x_p, y_p, pt_out, cv_inner = cv_inner)))
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
                pos_acc = self.ncv.recall_pos(y_p_t, all_preds)
                print "pos_acc", pos_acc
                #accuracy -1 class
                neg_acc = self.ncv.recall_neg(y_p_t, all_preds)
                print "neg_acc", neg_acc
                #balanced accuracy
                bacc = self.ncv.bacc(pos_acc, neg_acc)
                print "balanced acc", bacc
                #TODO get misclassified reconstructions samples
                miscl = y_p_t.index[(all_preds != y_p_t)]
                #bind actual labels and predictions
                #print miscl, y_p_t.loc[miscl], all_preds.loc[miscl]
                miscl_plus = ps.concat([y_p_t.loc[miscl], all_preds.loc[miscl]], axis = 1)
                #print miscl_plus
                self.write_miscl(model_out, pt_out, miscl_plus)
                f.write('%s\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc))
                f.flush()
            all_preds = ps.DataFrame(self.ncv.outer_cv(x,y, x_p = x_p, y_p = y_p, pt_out = pt_out))
            #temporary hack to the get random generator into place
            folds = self.setup_folds_synthetic(len(x), cv_outer)
            #folds_inner = self.setup_folds_synthetic(folds[1], 10)
            #folds_inner_extd = self.setup_folds_synthetic(folds[3], 10)
            #outer_subsampling = folds[0] * folds[1] + folds[2] * folds[3]
            #inner_subsampling = folds[0] * (folds_inner[0] * folds_inner[1] + folds_inner[2] * folds_inner[3]) * (len(self.config['c_params'] )+ 1)
            #for i in range(2 * (outer_subsampling + inner_subsampling + inner_subsampling_extd + outer_subsampling)):
            #    random.randint(0,9)
            #inner_subsampling_extd = folds[2] * (folds_inner_extd[0] * folds_inner_extd[1] + folds_inner_extd[2] * folds_inner_extd[3]) * (len(self.config['c_params']) + 1)
            #additional c param
            for i in range(8476 * (10 + 10 * 10)):
                random.randint(0,9)
            #skip nested cv
            #for i in range(8476 * (10 + 10 + 10 * 10 * (len(self.config['c_params']) + 1)) ):
            #    random.randint(0,9)
            self.ncv.majority_feat_sel(x, y, x_p, y_p, all_preds, self.config['k'], pt_out)
        f.close()
if __name__=="__main__":
    #only testing
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
-h if set the program tries to resume the classification / feature selection i.e. the output directory already exists
        """ % (sys.argv[0])
        sys.exit(2)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "f:d:y:l:a:bv:i:cj:g:p:o:s:r:eh")
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
    resume = False
    config_f = "/net/metagenomics/projects/phenotypes_20130523/code/learning/config.json"

    for o, a in optlist:
        if o == "-f":
            config_f = a
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
        if o == "-s":
            perc_samples = float(a)
        if o == "-r":
            perc_feats = float(a)
        if o == "-e":
            inverse_feats = True
        if o == "-h":
            resume = True
    print config_f
    if os.path.exists(out) and not resume:
        sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
        sys.exit(1)
    elif not os.path.exists(out):
        os.mkdir(out)
    pt_cl = pt_classification(config_f = config_f, phypat_f = phypat_f, gt_start = g1, gt_end = g2, pt_start = pt1, pt_end = pt2, rec_dir = rec_dir, likelihood_params = likelihood_params, parsimony_params = parsimony_params, is_phypat_and_rec = is_phypat_and_rec, cv_inner = cv_inner, cv_outer = cv_outer, model_out = out, n_jobs = n_jobs, perc_samples = perc_samples, perc_feats = perc_feats, inverse_feats = inverse_feats, do_normalization = do_normalization, resume = resume) 
