import nested_cv as ncv
import os.path
import sys
import pandas as pd
import random
import json
import warnings
#avoid seeing the deprecation warnings for scikit learn
warnings.filterwarnings("ignore")

class pt_classification:

    def write_miscl(self, model_out, pt_out, miscl_plus):
        """map the misclassified sample ids to their scientific names and taxonomic ids"""
        miscl_m = pd.read_csv(self.config["phyn_f"], sep = "\t", index_col = 0)
        miscl_plus.index = miscl_plus.index.astype('string')
        miscl_plus.columns = ["ground_truth", "predicted"]
        miscl_plus = miscl_plus.astype('int')
        miscl_m.index = miscl_m.index.astype('string')
        miscl_m.columns = ["sample_names"]
        pd.concat([miscl_plus, miscl_m.loc[miscl_plus.index,]], axis = 1).to_csv("%s/%s_miscl.txt"%(model_out, pt_out), sep = "\t")

    
    def setup_folds_synthetic(self, x, folds):
        """helper method to manipulate the random generator"""
        #short folds length
        f = x / folds
        #number of extended folds
        f_e = x % folds
        return (folds - f_e, x - f, f_e, ((x - (f + 1))))
    

    
    def __init__(self, config_f, model_out, pf2acc_desc_f, pt2acc_f, phypat_f, ids2name, rec_dir, likelihood_params, is_phypat_and_rec, cv_outer, cv_inner, n_jobs, perc_samples, perc_feats, inverse_feats, do_normalization, resume, tree, tree_named, parsimony_params = None, consider_in_recon = None, with_seed = False, is_discrete_phenotype_with_continuous_features = False, block_cross_validation = None):
        """main routine to prepare data for classification and feature selection"""
        if with_seed:
            random.seed(0)
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
        #check if we want to do nested cross validation accuracy estimation
        is_rec_based = False
        if not rec_dir is None:
            is_rec_based = True
        self.ncv = ncv.nested_cv(likelihood_params, parsimony_params, do_normalization, is_rec_based, is_phypat_and_rec, n_jobs, inverse_feats, self.config, perc_feats, perc_samples, model_out, cv_outer, resume, pf2acc_desc_f, consider_in_recon, is_discrete_phenotype_with_continuous_features, block_cross_validation)
        #write config to disk
        #get version
        from traitarm._version import __version__
        self.config['version'] =  __version__
        with open("%s/config.json" % self.model_out, 'w') as out_f:
            json.dump(self.config, out_f, indent=4, separators=(',', ': '))
        self.config['tree_f'] = tree
        self.config['tree_l_f'] = tree_named
        self.config['phyn_f'] = ids2name 
        #read in phyletic patterns
        p = pd.read_csv(phypat_f, sep = "\t", na_values = ["?"], index_col = 0)
        p.index = p.index.astype('string')
        #iterate over all phenotypes
        pf_mapping = pd.read_csv(pf2acc_desc_f, index_col = 0, sep = "\t") 
        pt_mapping = pd.read_csv(pt2acc_f, index_col = 0, sep = "\t")
        pt_mapping.index = pt_mapping.index.values.astype('string')
        for pt in pt_mapping.index:
            pt_out = pt
            x_p = p.loc[p.loc[:,pt].notnull()].loc[:, pf_mapping.index]
            y_p = p.loc[p.loc[:,pt].notnull()].loc[:, pt]
            #shuffle phyletic pattern samples to avoid biases in the cross validation
            rows = x_p.index.values.copy()
            random.shuffle(rows)
            x_p = x_p.loc[rows, ]
            y_p = y_p.loc[rows, ]
            #determine if we should infer ancestral feature gains and losses
            if not rec_dir is None:
                if not os.path.exists("%s/pt%s.dat"%(rec_dir, pt)):
                    print "skipping pt", pt, "no reconstruction matrix found, possibly due to no reconstruction events for this phenotype"
                    continue
                a = pd.read_csv("%s/pt%s.dat"%(rec_dir, pt), index_col = 0, sep = "\t")
                a.index = a.index.astype('string')
                #only relevant for the parsimony case
                if not parsimony_params is None :
                    raise Exception("not yet implemented")
                a.columns = pf_mapping.index.tolist() + [pt]
                #x annotation, y phenotype table
                x = a.loc[:, pf_mapping.index]
                y = a.loc[:, pt]
            else:
                #we are in pure phyletic pattern mode
                x = x_p
                y = y_p
            if not do_normalization:
                #binarize
                x_p = (x_p > 0).astype('int')
                if likelihood_params is None:
                    x = (x > 0).astype('int')
            #set phenotype negative class to -1
            y[(y==0)] = -1
            #skip phenotypes with too little samples
            if len(y_p) < self.config['min_samples']:
                print "skipping pt %s with only %s observations"%(pt_out,len(y_p))
                continue
            if sum(y_p > 0) < self.config['min_pos']:
                print "skipping pt %s with only %s + observations"%(pt_out,sum(y_p>0))
                continue
            if sum(y_p <=0) < self.config['min_neg']:
                print "skipping pt %s with only %s - observations"%(pt_out,sum(y_p<=0))
                continue
            print "running phenotype", pt, "with", sum(y_p>0), "phenotype-positive samples and", sum(y_p<=0), "phenotype-negative samples"
            if not cv_inner is None:
                try:
                    outcome = self.ncv.outer_cv(x,y, x_p, y_p, pt_out, cv_inner = cv_inner)
                    all_preds, all_scores  = pd.Series(pd.np.array(outcome[0])), outcome[1]
                except ValueError as e:
                    import traceback
                    print traceback.print_exc(e)
                    print e
                    print "pt %s has too few samples in one of the classes possibly, due to a high threshold in the likelihood reconstruction; try to increase the number of cross validation fold and run again"
                    continue
                all_preds.index = y_p.index
                y_p_t = y_p.copy()
                y_p_t[y_p_t == 0] = -1
                pos_acc = self.ncv.recall_pos(y_p_t, all_preds)
                print "pos_acc", pos_acc
                #accuracy -1 class
                neg_acc = self.ncv.recall_neg(y_p_t, all_preds)
                print "neg_acc", neg_acc
                #balanced accuracy
                bacc = self.ncv.bacc(pos_acc, neg_acc)
                print "balanced acc", bacc
                precision = self.ncv.precision(y_p_t, all_preds)
                print "precision", precision 
                f1_score = self.ncv.f1_score(precision, pos_acc)
                print "f1_score", f1_score
                #TODO get misclassified reconstructions samples
                miscl = y_p_t.index[(all_preds != y_p_t)]
                #bind actual labels and predictions
                miscl_plus = pd.concat([y_p_t.loc[miscl], all_preds.loc[miscl]], axis = 1)
                self.write_miscl(model_out, pt_out, miscl_plus)
                #cv accuracy stats
                cv_out = "%s/cv_acc.txt"%self.model_out
                with open(cv_out, "a") as f:
                    if not os.path.exists(cv_out):
                        f.write('\tpos_recall\tneg_recall\tbalanced_accuracy\tprecision\tf1-score\n')
                with open(cv_out, "a") as f:
                    f.write('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (pt_out, pos_acc, neg_acc, bacc, precision, f1_score))
                    f.flush()
                #auc/roc 
                auc = self.ncv.roc_curve(y_p, all_scores, "%s/%s_roc_curve.png" % (self.model_out, pt_out), pos_acc, 1 - neg_acc)
                print "auc ", auc
            all_preds, all_scores  = self.ncv.outer_cv(x,y, x_p = x_p, y_p = y_p, pt_out = pt_out, do_calibration = True)
            all_preds = pd.DataFrame(all_preds)
            all_scores = pd.DataFrame(all_scores)
            folds = self.setup_folds_synthetic(len(x), cv_outer)
            self.ncv.majority_feat_sel(x, y, x_p, y_p, all_preds, all_scores, self.config['k'], pt_out)

