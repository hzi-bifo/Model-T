import numpy as np
import pandas as pd
#learning
import sklearn.svm as svm
from sklearn import tree 
#standardization
import sklearn.preprocessing as preprocessing
#ancestral reconstruction of phenotype
import cv_rec_helper as crh 
#cross validation
from sklearn.model_selection import KFold, GroupKFold
#feature importance
import scipy.stats
#auc/roc + probability calibration#auc/roc
import sklearn.calibration as clb
import matplotlib
#heatmap
#avoid using X display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
#make sure figures produced in this script are not cut
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})
from sklearn.metrics import roc_curve, auc

from operator import itemgetter
import math
import random
import sys
import itertools
import os


class nested_cv:

    def __init__(self, likelihood_params, parsimony_params, do_standardization, is_rec_based, is_phypat_and_rec, n_jobs, inverse_feats, config, perc_feats, perc_samples, model_out, cv_outer, resume, pf2desc_f, consider_in_recon, is_discrete_phenotype_with_continuous_features, block_cross_validation, opt_measure):
        self.config = config
        self.likelihood_params = likelihood_params
        self.parsimony_params = parsimony_params
        self.do_standardization = do_standardization
        self.is_rec_based = is_rec_based
        self.n_jobs = n_jobs
        self.inverse_feats = inverse_feats
        self.is_phypat_and_rec = is_phypat_and_rec
        self.perc_feats = perc_feats
        self.perc_samples = perc_samples
        self.model_out = model_out
        self.cv_outer = cv_outer
        self.resume = resume
        self.pf2desc_f = pf2desc_f
        self.consider_in_recon = consider_in_recon
        self.is_discrete_phenotype_with_continuous_features  = is_discrete_phenotype_with_continuous_features
        self.opt_measure = opt_measure
        if block_cross_validation is not None:
            self.block_cross_validation = pd.read_csv(block_cross_validation, index_col = 0, sep = "\t")
        else:
            self.block_cross_validation = None

    def transf_from_probs(self, x, y):
        """create a series of sample specific weights according to the probabilities"""
        y_extd = [] 
        x_index = []
        for yi in range(len(y)):
            if y[yi] == 0 or y[yi] == -1:
                y_extd.append(-1)
                x_index.append(yi)
            else:
                y_extd.append(1)
                x_index.append(yi)
                if "extend" in self.likelihood_params:
                    y_extd.append(-1)
                    x_index.append(yi)
        w = []
        for yi in range(len(y)):
            if y[yi] == 0 or y[yi] == -1:
                w.append(1)
            else:   
                w.append(y[yi])
                if "extend" in self.likelihood_params:
                    w.append(1 - y[yi])
        return x.iloc[x_index,], pd.Series(y_extd), pd.Series(w)
    
    def standardize(self, array):
        """standardize a feature matrix"""
        scaler = preprocessing.StandardScaler(with_mean = True, with_std = True).fit(array)
        transformed = pd.DataFrame(data = scaler.transform(array), index = array.index, columns = array.columns)
        return transformed, scaler
    
    @staticmethod 
    def confusion_m(y,y_pred):
        """get confusion matrix TN FP /n FN TP"""
        TP = (y[y == 1] == y_pred[y == 1]).sum()
        FN = (y[y == 1] != y_pred[y == 1]).sum()
        TN = (y[y == -1] == y_pred[y == -1]).sum()
        FP = (y[y == -1] != y_pred[y == -1]).sum()
        return pd.np.array([TN, FP, FN, TP])

    @staticmethod 
    def recall_pos(y,y_pred):
        """compute recall of the positive class"""
        return (y[y == 1] == y_pred[y==1]).sum()/float((y==+1).sum())
    
    @staticmethod 
    def recall_pos_conf(conf):
        """compute recall of the positive class"""
        TN, FP, FN, TP = conf
        if (TP + FN) == 0:
            float('nan') 
        return TP/float(TP+FN)
    
    @staticmethod 
    def recall_neg_conf(conf):
        """compute recall of the positive class"""
        TN, FP, FN, TP = conf
        if (TN + FP) == 0:
            return float('nan')
        return TN/float(TN+FP)

    @staticmethod 
    def recall_neg(y, y_pred):
        """compute recall of the negative class"""
        return (y[y == -1] == y_pred[y==-1]).sum()/float((y==-1).sum())
    
    @staticmethod 
    def precision(y, y_pred):
        """compute precision"""
        TP = (y[y == 1] == y_pred[y == 1]).sum()
        FP = (y[y == -1] != y_pred[y == -1]).sum()
        if (TP + FP) == 0:
            return 0
        return TP / float(TP + FP)   
    
    @staticmethod 
    def npv(y, y_pred):
        """compute precision"""
        FN = (y[y == 1] != y_pred[y == 1]).sum()
        TN = (y[y == -1] == y_pred[y == -1]).sum()
        if (TN + FN) == 0:
            return 0
        return TN / float(TN + FN)   

    @staticmethod 
    def f1_score(recall, precision):
        """compute f1-measure"""
        if (precision + recall) != 0:
            return 2 * (precision * recall) / (precision + recall)    
        return 0
    
    @staticmethod 
    def f1_score_neg(recall_neg, npv):
        """compute negative f1-measure"""
        if (npv + recall_neg) != 0:
            return 2 * (npv * recall_neg) / (npv + recall_neg)    
        return 0

    @staticmethod 
    def precision_conf(conf):
        """compute precision"""
        TN, FP, FN, TP = conf
        if (TP + FP) == 0:
            return float('nan')
        return TP / float(TP + FP)

    @staticmethod
    def bacc(pos_acc, neg_acc):
        """compute balanced accuracy"""
        return float(pos_acc + neg_acc)/2
    
    def roc_curve(self, y, all_scores, out, tpr_opt, fpr_opt):
        """plot roc curve and compute auc"""
        #auc / roc
        fpr, tpr, _ = roc_curve(y, all_scores)
        roc_auc = auc(fpr, tpr)
        plt.figure()
        lw = 2
        fig, ax = plt.subplots()
        plt.plot(fpr, tpr, color='darkorange',lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        ax.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        ax.scatter([fpr_opt],[tpr_opt])
        ax.annotate("bacc opt %.2f " % (((1 - fpr_opt) + tpr_opt)/2), (fpr_opt, tpr_opt)) 
        plt.title('Receiver operating characteristic')
        plt.legend(loc="lower right")
        plt.savefig(out)
        plt.close()
        return roc_auc
        
    def balance_weights(self, weights, y):
        """balance weights between pos/neg class in phyletic pattern and in the reconstruction based samples """
        weights.index = y.index
        weights[y == -1]  = weights[y == -1] / weights[y == -1].sum() 
        weights[y == 1]  = weights[y == 1] / weights[y == 1].sum()
        return weights
    
    
    def cv(self, x_train, y_train, xp_train, yp_train, x_test , C, no_classifier = 1, do_calibration = False):
        """train model on given training features and target and return the predicted labels for the left out samples"""
        if not self.likelihood_params is None and "continuous_target" in self.likelihood_params:
            #transform x and y
            x_train, y_train, w = self.transf_from_probs(x_train, y_train, self.likelihood_params)
        else:
            w = pd.Series(np.ones(shape = len(y_train)) )
        #set y negative labels to -1
                #EXPERIMENTAL BALANCE WEIGHTS
                #sample_weight = pd.concat([balance_weights(w, y_train_t), balance_weights(pd.Series(n, sep = "\t"p.ones(len(yp_train_t))), yp_train_t)])
                #print sample_weight.sum(), "sample weights total"
                #END EXPERIMENTAL BALANCE WEIGHTS
        y_train.index = x_train.index
        w.index = x_train.index
        y_train_t = y_train.copy() 
        y_train_t[y_train_t == 0] = -1
        yp_train_t = yp_train.copy() 
        yp_train_t[yp_train_t == 0] = -1
        predictor = svm.LinearSVC(C=C)
        predictor.set_params(**self.config['liblinear_params'])
        #check if we are in vanilla linear SVM and set no_classifiers to 1 if so or in subspace mode 
        if self.perc_feats == 1.0 and self.perc_samples == 1.0:
            no_classifier = 1
        all_preds = pd.DataFrame(np.zeros(shape = (x_test.shape[0], no_classifier)))
        all_scores = pd.DataFrame(np.zeros(shape = (x_test.shape[0], no_classifier)))
        for i in range(no_classifier):
            #select a subset of features for classification
            sample_feats = sorted(random.sample(x_train.columns, int(math.floor(x_train.shape[1] * self.perc_feats))))
            #select a subset of samples for classification
            sample_samples = sorted(random.sample(x_train.index, int(math.floor(x_train.shape[0] * self.perc_samples))))
            sample_samples_p = sorted(random.sample(xp_train.index, int(math.floor(xp_train.shape[0] * self.perc_samples))))
            #reduce feature space to the selected features
            x_train_sub = x_train.loc[sample_samples, sample_feats].copy()
            y_train_t_sub = y_train_t.loc[sample_samples]
            yp_train_t_sub = yp_train_t.loc[sample_samples_p]
            #print y_train_t_sub, yp_train_t_sub, w_sub
    
            if self.is_phypat_and_rec:
                #EXPERIMENTAL BALANCE WEIGHTS
                #sample_weight = pd.concat([balance_weights(w, y_train_t), balance_weights(pd.Series(n, sep = "\t"p.ones(len(yp_train_t))), yp_train_t)])
                #print sample_weight.sum(), "sample weights total"
                #END EXPERIMENTAL BALANCE WEIGHTS
                #reduce xp feature space to the selected features
                xp_train_sub = xp_train.loc[sample_samples_p, sample_feats]
                X = pd.concat([x_train_sub, xp_train_sub], axis = 0)
                if self.inverse_feats:
                    #add inverse features
                    X = pd.concat([X,1-X], axis = 1) 
                if self.do_standardization:
                    X, scaler = self.standardize(X)
                y = pd.concat([y_train_t_sub, yp_train_t_sub], axis = 0)
                predictor.fit(X = X, y = y)
                if do_calibration:
                    cclf = clb.CalibratedClassifierCV(predictor, method = 'sigmoid', cv = 'prefit')
                    cclf.fit(X, y)
            else: 
                if self.inverse_feats:
                    #add inverse features
                    x_train_sub = pd.concat([x_train_sub, 1 - x_train_sub], axis = 1) 
                if self.do_standardization:
                    #if reconstruction based use absolute change for prediction
                    if self.is_rec_based:
                        x_train_sub, scaler = self.standardize(x_train_sub.abs())
                    else:
                        x_train_sub, scaler = self.standardize(x_train_sub.abs())
                predictor.fit(x_train_sub, y_train_t_sub)
                #if learning from gene gains and losses, get selected feature and rebuild model based on phyletic patterns
                if self.is_rec_based:
                    models = pd.DataFrame(np.zeros(shape=(x_train.shape[1], 1)))
                    models.index = x_train_sub.columns
                    models.iloc[:, 0] = predictor.coef_[0]
                    xp_train_sub = xp_train.loc[sample_samples_p, sample_feats]
                    xp_train_sub_t = xp_train_sub.copy()
                    xp_train_sub_t.loc[:, ~models.apply(lambda x: (x > 0).sum() >= 1 or (x < 0).sum() >= 1, axis = 1) ] = 0
                    if self.do_standardization:
                        xp_train_sub_t, scaler = self.standardize(xp_train_sub_t)
                    predictor.fit(xp_train_sub_t, yp_train_t_sub)
                #probability calibration
                    if do_calibration:
                        cclf = clb.CalibratedClassifierCV(predictor, method = 'sigmoid', cv = 'prefit')
                        cclf.fit(xp_train_sub_t, yp_train_t_sub)
                else:
                    if do_calibration:
                        cclf = clb.CalibratedClassifierCV(predictor, method = 'sigmoid', cv = 'prefit')
                        cclf.fit(x_train_sub, y_train_t_sub)
            #copy test sample before modification
            x_test_sample = x_test.loc[:, sample_feats].copy()
            if self.inverse_feats:
                #add inverse features to test sample
                x_test_sample = pd.concat([x_test_sample, 1 - x_test_sample], axis = 1) 
            #standardize test sample
            if self.do_standardization:
                x_test_sample = pd.DataFrame(data = scaler.transform(x_test_sample), index = x_test_sample.index, columns = x_test_sample.columns)
            all_preds.iloc[:, i]  = predictor.predict(x_test_sample)
            #use calibration classifier
            if do_calibration:
                uncalibrated  = predictor.decision_function(x_test_sample)
                all_scores.iloc[:, i]  = pd.Series(cclf.predict_proba(x_test_sample)[:, 1])
        #do majority vote to aggregate predictions
        aggr_preds = all_preds.apply(lambda x: 1 if sum(x == 1) > sum(x == -1) else -1, axis = 1).astype('int')
        return aggr_preds, all_scores 
    
    def get_training_and_testing_mapping(self, train_pt_edges, all_pt_edges):
        """map phenotype edges that are not in the full reconstruction and thus must have been joined"""
        train2all = {}
        for e in train_pt_edges:
            nodes = e.split("_")
            for e2 in all_pt_edges:
                nodes2 = e2.split("_")
                maps = True
                for n2 in nodes2:
                    if n2 not in nodes:
                        maps = False
                if maps: 
                    if not e in train2all: 
                        train2all[e] = [e2]
                    else: 
                        train2all[e].append(e2)
        test_edges = all_pt_edges.intersection(set(itertools.chain.from_iterable(train2all.values())))
        return train2all, test_edges
    
    def join_edges(self, x_r, train2all, likelihood_params, parsimony_params):
        """aggregate edges in the full reconstruction matrix to fit the joint edges in train2all"""
        df = pd.DataFrame()
        for e in train2all:
            #aggregate rows by summing up
            if not self.likelihood_params is None:
                #if features have been discretized before just sum over the values of the branches
                if not 'gt_probs' in self.likelihood_params:
                    cs = x_r.loc[train2all[e], :].sum(axis=0)
                    df = pd.concat([df, cs], axis = 1)
                #in case of continuous features sum over the branches (the feature value of a branch represents the difference of the reconstructed value of the feature for the start and the end node of that branch)
                elif self.is_phenotype_and_continuous_features:
                    cs = np.zeros(shape = x_r.shape[1])
                    cs = cs + x_r.loc[train2all[e], :].iloc[i, :]
                    df = pd.concat([df, pd.Series(cs)], axis = 1)
                #sum up the probabilities
                else:
                    cs = np.zeros(shape = x_r.shape[1])
                    for i in range(x_r.loc[train2all[e], :].shape[0]):
                        #print "before aggregation", cs 
                        #print "being aggregated", x_r.loc[train2all[e], :].iloc[i, :]
                         cs = cs + (1 - cs) * x_r.loc[train2all[e], :].iloc[i, :]
                        #print "after aggregation", cs
                    df = pd.concat([df, pd.Series(cs)], axis = 1)
            else: 
                raise Exception("parsimony case not yet implemented")    
        #make sure there are no entries in the summed up likelihood matrix other than 0 and 1
        if not self.likelihood_params is None and 'gt_probs' in self.likelihood_params:
            df[x_r > 0] = 1
        #get df back into samples x features shape 
        df.columns = train2all.keys()
        return df.T
    
    
    def get_rec_samples(self, x_r, y_r, yp_train, pt_out, ofold, ifold):
        """compute the training samples by running parsimony/likelihood reconstruction for the pt for the training set without the held out samples"""
        #retrieve matrix like:
        #N1_N4 0 
        #N2_44 1
        if not self.parsimony_params is None:
            #This is just an implementation stump
            m = crh.reconstruct_pt_parsimony(yp_train, self.model_out, self.config, self.parsimony_params)
        else:
            m = crh.reconstruct_pt_likelihood(yp_train, self.model_out, self.config, self.likelihood_params, pt_out, ofold, ifold, self.consider_in_recon)
        all_pt_edges = set(x_r.index)
        train_pt_edges = set(m.index)
        #print len(all_pt_edges), "anzahl aller reconstruction samples"
        #print len(train_pt_edges), "anzahl aller reconstruction samples ohne test samples"
        #print yp_train.index.values, "training samples"
        #all edges only present in the full reconstruction
        df1 = all_pt_edges.difference(train_pt_edges)
        #print "held out edges including not yet joined edges", df1
        #all edges that are only in the training reconstruction and thus must have been joined
        df2 = train_pt_edges.difference(all_pt_edges)
        #print "edges that are only found in the training reconstruction", df2
        #print "training edges", train_pt_edges
        #map training phenotype edges which don't have a 1:1 equivalent in the all_pt_edges 
        train2all, test_edges = self.get_training_and_testing_mapping(df2, df1)
        #print "mapping from edges only in the training reconstruction to edges in the full reconstruction", train2all
        #print "held out edges from get training and testing mapping", test_edges
        #get the corresponding samples in the full reconstruction matrix and aggregate them
        joint_x_r = self.join_edges(x_r, train2all, self.likelihood_params, self.parsimony_params)
        #add the y labels to those rows
        joint_y_r = m.loc[train2all, :].iloc[:, 0] 
        #add all the other training samples
        training_edges = all_pt_edges.intersection(train_pt_edges)
        #print x_r.loc[training_edges,:].shape
        #print joint_x_r.shape, joint_x_r.index
        x_r_train = pd.concat([x_r.loc[training_edges,:], joint_x_r], axis = 0)
        #print "x_r_train shape", x_r_train.shape
        #print "x_r_train index", x_r_train.index
        y_r_train = pd.concat([m.loc[training_edges, :].iloc[:,0], joint_y_r], axis = 0)
        y_r_train[y_r_train == 2] = 1
        #print "y_r_train shape", y_r_train.shape
        #replace negative labels with -1
        x_r_test = x_r.loc[test_edges,]
        y_r_test = y_r.loc[test_edges,]
        return x_r_train, y_r_train, x_r_test, y_r_test
    
    
    def setup_folds(self, x, y, cv, x_p, y_p, pt_out, ofold = None, do_block_cross_validation = True):
        """prepare all combinations of training and test folds, either for the standard phyletic pattern or the likelihood case or for the combined case"""
        #divide the samples into cv folds and yield training and test fold
        folds = []
        #block cross validation
        if self.block_cross_validation is not None and do_block_cross_validation:
            kf = GroupKFold(n_splits = cv)
            kf_split = kf.split(x, groups = self.block_cross_validation.loc[x.index, "group_id"].tolist())
        #normal k fold cross-validation
        else:
            kf = KFold(n_splits = cv)
            kf_split = kf.split(x)
        #standard phyletic pattern cv
        for train, test in kf_split:
            x_train = x.iloc[train, :].copy()
            y_train = y.iloc[train].copy()
            x_test = x.iloc[test, :].copy()
            if not self.is_rec_based:
                yield ((x_train, y_train, x_test, None, x_train, y_train, x_test))
            #otherwise do max likelihood reconstruction
            else:
                xp_train, yp_train, xp_test = (x_train, y_train, x_test) 
                if ofold is None:
                    x_train, y_train, x_test, y_test = self.get_rec_samples(x, y, yp_train, pt_out, ofold = i, ifold = None)
                else:
                    x_train, y_train, x_test, y_test  = self.get_rec_samples(x, y, yp_train, pt_out, ofold = ofold, ifold = i)
                yield (x_train, y_train, x_test, y_test,  xp_train, yp_train, xp_test)
    
    
    def outer_cv(self, x, y, x_p, y_p, pt_out, cv_inner = None, do_calibration = True):
        """do a cross validation with cv_outer folds
        if only cv_outer is given do a cross validation for each value of the paramter C given
        if cv_inner is given, do a nested cross validation optimizing the paramter C in the inner loop"""
        #do n fold cross validation
        if not cv_inner is None:
            #DEBUG print "outer folds set up"
            #list of predictions
            ocv_preds = pd.Series(np.zeros(shape=[len(y_p)]))
            ocv_preds.index = x_p.index
            ocv_scores = ocv_preds.copy()
            i = 0
            for x_train, y_train, x_test, y_test,  xp_train, yp_train, xp_test in self.setup_folds(x, y, self.cv_outer, x_p, y_p, pt_out):
                #DEBUG print "inner folds set up"
                #set up prediction data frame
                all_preds = pd.DataFrame(np.zeros(shape=[len(yp_train), len(self.config['c_params'])]))
                all_preds.columns = self.config['c_params']
                all_preds.index = yp_train.index
                #do inner cross validation
                for x_train_train, y_train_train, x_train_test, y_train_test, xp_train_train, yp_train_train, xp_train_test in self.setup_folds(x_train, y_train, cv_inner, xp_train, yp_train, pt_out, i, do_block_cross_validation = False):
                    for c_param in self.config['c_params']:
                        preds, _ = self.cv(x_train_train, y_train_train, xp_train_train, yp_train_train, xp_train_test, c_param)
                        all_preds.loc[x_train_test.index, c_param] = list(preds)
                #determine C param with the best accuracy
                #copy yp train and change the negative labels from 0 to -1
                yp_t_t = yp_train.copy()
                yp_t_t[yp_t_t == 0] = -1
                perf_m = self.perf_evaluation(yp_t_t, all_preds)
                perf_m = perf_m.sort_values(by = self.opt_measure, axis = 1, ascending = False)
                c_opt = perf_m.columns[0]
                #DEBUG print c_opt, "optimal value of the C param is this fold"
                #use that C param to train a classifier to predict the left out samples in the outer cv
                p, score = self.cv(x_train, y_train, xp_train, yp_train, xp_test, c_opt, do_calibration = do_calibration)
                ocv_preds.loc[xp_test.index] = list(p)
                ocv_scores.loc[xp_test.index] = list(score.iloc[:, 0])
                i += 1
            return ocv_preds, ocv_scores
        #store predictions for each fold in all_preds
        all_preds = pd.DataFrame(np.zeros(shape=[len(y_p), len(self.config['c_params'])]))
        all_preds.index = x_p.index
        all_preds.columns = self.config['c_params']
        all_scores = all_preds.copy()
        #do a cross validation for each value of the C param
        for x_train, y_train, x_test, y_test, xp_train, yp_train, xp_test in self.setup_folds(x, y, self.cv_outer, x_p, y_p, pt_out):
            for c_param in self.config['c_params']:
                preds, scores = self.cv(x_train, y_train, xp_train, yp_train, xp_test, c_param, do_calibration = do_calibration)
                all_preds.loc[xp_test.index, c_param] = list(preds)
                all_scores.loc[xp_test.index, c_param] = list(scores.iloc[:, 0]) 
        return all_preds, all_scores
        
    def perf_evaluation(self, gold_standard, predictions):
        """evaluate different performance measures"""
        colnames = ['bacc', "pos-rec", "neg-rec", "precision", "F1-score", "neg-F1-score"]
        perf_m = pd.DataFrame(pd.np.zeros((len(colnames), len(self.config['c_params']))), index = colnames, columns = self.config['c_params'])
        for j in range(len(self.config['c_params'])):
            c_param = self.config['c_params'][j]
            #recall of pt positive class
            perf_m.loc['pos-rec', c_param] = self.recall_pos(gold_standard, predictions.iloc[:,j]) 
            #recall of pt negative class
            perf_m.loc['neg-rec', c_param] =  self.recall_neg(gold_standard, predictions.iloc[:,j]) 
            #balanced accuracy
            perf_m.loc['bacc', c_param] = self.bacc(perf_m.loc['pos-rec', c_param], perf_m.loc['neg-rec', c_param])
            #precision
            perf_m.loc['precision', c_param] = self.precision(gold_standard, predictions.iloc[:,j]) 
            #negative predictive value
            perf_m.loc['npv', c_param] = self.npv(gold_standard, predictions.iloc[:,j]) 
            #f1 score
            perf_m.loc['F1-score', c_param] = self.f1_score(perf_m.loc['pos-rec', c_param], perf_m.loc['precision', c_param]) 
            perf_m.loc['neg-F1-score', c_param] = self.f1_score(perf_m.loc['neg-rec', c_param], perf_m.loc['npv', c_param]) 
        return perf_m
        
    def majority_feat_sel(self, x, y, x_p, y_p, all_preds, all_scores, k, pt_out, no_classifier = 10):
        """determine the features occuring in the majority of the k best models"""
        #sanity check if number of classifiers selected for majority feature selection exceeds the number of c params
        if k > len(self.config['c_params']):
            sys.exit("number of selected classifiers for feature selection (%s) exceeds the number of c params (%s)" %(k, len(self.config['c_params'])))
        #retrieve weights from likelihood_params
        if not self.likelihood_params is None and "continuous_target" in self.likelihood_params:
            #transform x and y
            x, y, w = self.transf_from_probs(x, y)
        else:
            w = pd.Series(np.ones(shape = len(y)) )
        #discard non binary phenotype labels in the reconstruction matrix and target matrix
        if not self.likelihood_params is None:
            t = float(self.likelihood_params["threshold"])
            condition = ~((y != -1) & (y < t))
            y = y.loc[condition]
            x = x.loc[condition, :]
        #determine the k best classifiers
        #make sure the indices match
        all_preds.index = y_p.index
        x.columns = x_p.columns
        y.index = x.index
        #w.index = x.index
        y_p_t = y_p.copy()
        #make sure the negative class label -1 instead of 0
        y_p_t[y_p_t == 0] = -1
        y_t = y.copy()
        y_t[y_t == 0] = -1
        #prepare performance statistics
        perf_m = self.perf_evaluation(y_p_t, all_preds)
        #sort by optimal  measure
        perf_m_k = perf_m.sort_values(by = self.opt_measure, axis = 1, ascending = False).iloc[:, :k]
        #roc curves / auc
        all_scores.columns = self.config['c_params']
        auc_df = pd.DataFrame(pd.np.zeros(k), index = perf_m_k.columns, columns = ['auc'])
        for c_param, pos_rec, neg_rec in [(perf_m_k.columns[i], perf_m_k.loc['pos-rec'].iloc[i], perf_m_k.loc['neg-rec'].iloc[i])  for i in range(k)]:
            auc_df.loc[c_param] = self.roc_curve(y_p, all_scores.loc[:, c_param], "%s/%s_%s_roc_curve.png" %(self.model_out, pt_out, c_param), pos_rec, 1 - neg_rec)
        #add auc to performance summary
        perf_m_k = pd.concat([perf_m_k, auc_df.T], axis = 0)
        #write them to disk
        perf_m_k.to_csv("%s/%s_perf.txt"%(self.model_out,pt_out), sep="\t", float_format = '%.3f')
        predictors = []
        #check if we are in vanilla linear SVM and set no_classifiers to 1 if so or in subspace mode 
        if self.perc_feats == 1.0 and self.perc_samples == 1.0:
            no_classifier = 1
        models = pd.DataFrame(np.zeros(shape=(x.shape[1], k * no_classifier)))
        models.index = x.columns
        bias_cparams = []
        for i in range(k):
            predictor = svm.LinearSVC(C=perf_m_k.columns[i])
            predictor.set_params(**self.config["liblinear_params"])
            #subset the features
            sample_feats = sorted(random.sample(x.columns, int(math.floor(x.shape[1] * self.perc_feats))))
            #subset the samples
            sample_samples_p = sorted(random.sample(x_p.index, int(math.floor(x_p.shape[0] * self.perc_samples))))
            sample_samples = sorted(random.sample(x.index, int(math.floor(x.shape[0] * self.perc_samples))))
            #train the full model with the k best models
            for l in range(no_classifier):
                x_sub = x.loc[sample_samples, sample_feats]
                y_t_sub = y_t.loc[sample_samples]
                #check if using ancestral phenotype gains and losses and phyletic patterns combined
                if self.is_phypat_and_rec:
                    x_p_sub = x_p.loc[sample_samples_p, sample_feats]
                    y_p_t_sub = y_p_t.loc[sample_samples_p]
                    X = pd.concat([x_sub, x_p_sub], axis = 0)
                    #add inverse features if the corresponding option is set
                    if self.inverse_feats:
                        X = pd.concat([X, 1 - X], axis = 1)
                    #standardize if the corresponding option is set
                    if self.do_standardization:
                        X, _ = self.standardize(X)
                        pass
                    #Start EXPERIMENTAL this works only with the fork of scikit learn that supports sample weights
                    #predictor.fit(X, pd.concat([y_t, y_p_t], axis = 0), sample_weight = pd.concat([balance_weights(w, y_t), balance_weights(pd.Series(np.ones(shape = len(y_p))), y_p_t)]))
                    #END EXPERIMENTAL this works only with the fork of scikit learn that supports sample weights
                    predictor.fit(X, pd.concat([y_t_sub, y_p_t_sub], axis = 0))
                else:
                    if self.inverse_feats:
                        #add inverse features
                        x_sub =  pd.concat([x_sub, x_sub], axis = 1)
                    #standardize if , reduce = Truethe corresponding option is set
                    if self.do_standardization:
                        x_sub, _ = self.standardize(x_sub.abs())
                        #x_sub = x_sub.abs()
                        pass
                    predictor.fit(x_sub, y_t_sub)
                    if self.is_rec_based:
                        models_t = pd.DataFrame(np.zeros(shape=(x.shape[1], 1)))
                        models_t.index = x_sub.columns
                        models_t.iloc[:, 0] = predictor.coef_[0]
                        xp_sub = x_p.loc[sample_samples_p, sample_feats]
                        xp_sub_t = xp_sub.copy()
                        xp_sub_t.loc[:, ~models_t.apply(lambda x: (x > 0).sum() >= 1 or (x < 0).sum() >= 1, axis = 1) ] = 0
                        if self.do_standardization:
                            xp_sub_t, _ = self.standardize(xp_sub_t) 
                        y_p_t_sub = y_p_t.loc[sample_samples_p]
                        predictor.fit(xp_sub_t, y_p_t_sub)
                #add inverse features if the corresponding option is set
                if self.inverse_feats:
                    print "in inverse features mode"
                    #discard negative features in first half of the weight vector and negative features in the second half of the weight vector
                    pos_coef = predictor.coef_[0][0:(len(predictor.coef_[0]))/2]
                    pos_coef[pos_coef < 0] = 0
                    neg_coef = predictor.coef_[0][(len(predictor.coef_[0]))/2 : len(predictor.coef_[0])]
                    neg_coef[neg_coef > 0] = 0
                    rel_weights = pos_coef + neg_coef
                else:
                    rel_weights = predictor.coef_
                predictors.append(predictor)
                models.loc[sample_feats, i * no_classifier + l] = rel_weights 
            bias_cparams.append((predictor.C, predictor.intercept_[0]))
        pd.DataFrame(bias_cparams).to_csv("%s/%s_bias.txt"%(self.model_out,pt_out), sep = "\t", index = None, header = None)
        feats = []
        #determine features with non-zero weights
        for i in range(models.shape[0]):
            if sum(models.loc[models.index[i],:] != 0) > 1:
                feats.append(models.index[i])
        #write all the features with their weights to disk
        rownames_extd = [str(perf_m_k.columns[i]) + "_" + str(l) for i in range(k) for l in range(no_classifier)] 
        models_df = pd.DataFrame(models)
        models_df.columns = rownames_extd
        for i in range(k):
            models_df[models_df.columns[i]] = models_df[models_df.columns[i]].map(lambda x: '%.3f' % x)
        models_df.to_csv("%s/%s_feats.txt"%(self.model_out,pt_out), sep="\t")
        #standardize the features and write standardization parameters to disk
        if self.do_standardization:
            _, scaler = self.standardize(x_p)
            scale_df = pd.DataFrame(scaler.scale_, index = x_p.columns, columns = ["scale"]).to_csv("%s/%s_scale.txt" % (self.model_out, pt_out), sep = "\t")
            scale_df = pd.DataFrame(scaler.mean_, index = x_p.columns, columns = ["mean"]).to_csv("%s/%s_mean.txt" % (self.model_out, pt_out), sep = "\t")
        #get baseline classification models for each individual feature
        if not len(feats) == 0:
            #initiate, fit and predict with decision stump for each feature
            preds = [tree.DecisionTreeClassifier(max_depth = 1, class_weight = 'balanced').fit(pd.DataFrame(x_p.loc[:, i]), y_p_t).predict(pd.DataFrame(x_p.loc[:, i])) for i in feats]
            #get confusion matrix
            conf_per_feat = pd.DataFrame([self.confusion_m(y_p_t, pd.Series(p, index = y_p.index).T) for p in preds ]) 
            conf_per_feat.index = feats 
            conf_per_feat.columns = ["TN", "FP", "FN", "TP"]
            #get macro accuracy
            bacc = conf_per_feat.apply(lambda x: self.bacc(self.recall_pos_conf(x), self.recall_neg_conf(x)), axis = 1)
            perf_per_feat = pd.concat([conf_per_feat, bacc], 1)
            perf_per_feat.columns = ["TN", "FP", "FN", "TP"] + ["MACC"]
            #write majority features with their weights to disk
            pf2desc = pd.read_csv(self.pf2desc_f, sep = "\t", index_col = 0).iloc[:, 0]
            feat_df = pd.concat([pf2desc.loc[feats, ], models_df.loc[feats, ], perf_per_feat], axis = 1)
            feat_df.columns = ["description"] + rownames_extd + ["TN", "FP", "FN", "TP", "MACC"]
            columns_out = rownames_extd + ["description"] + ["TN", "FP", "FN", "TP", "MACC"]
            feat_df.sort_values(by = ["MACC"], axis = 0, ascending = False).to_csv("%s/%s_non-zero+weights.txt"%(self.model_out,pt_out), columns = columns_out, float_format='%.3f',  sep = "\t")
            #produce heatmaps
            self.feat_heatmap(models_df.loc[feats, ], perf_per_feat, x_p.loc[:, feats], y_p_t, pt_out)

    def feat_heatmap(self, models_df, perf_per_feat, feat_df, y, pt_out):
        """make a heatmap of feature occurence across samples and feature SVM weights"""
        #feature occurence heatmap
        feats_sorted = perf_per_feat.loc[:, "MACC"].sort_values()
        #sort by phenotype positive and negative samples
        y_sorted = y.sort_values()
        target_df = feat_df.loc[y_sorted.index, feats_sorted.index]
        target_df.index = target_df.index.astype('string')
        phn_f = pd.read_csv(self.config["phyn_f"], sep = "\t", index_col = 0)
        phn_f.index = phn_f.index.astype('string')
        phn_name2id = pd.read_csv(self.config["phyn_f"], sep = "\t", index_col = 1)
        phn_name2id.index = phn_name2id.index.astype('string')
        target_df.index = phn_f.loc[target_df.index, :].iloc[:, 0]
        target_df.index.name = ""
        seaborn.set()
        f, (ax1, ax2) = plt.subplots(2, sharex = True, gridspec_kw = {'height_ratios': [3, 1]})
        #restrict to top 20 features
        hm = seaborn.heatmap(target_df.iloc[:, -20:], ax = ax1, cbar_kws={'label': 'Feature occurence'})
        #feature weight heatmap
        models_df_sorted = models_df.T.loc[:, feats_sorted.index]
        hm2 = seaborn.heatmap(models_df_sorted.astype('float').iloc[:, -20:], ax = ax2, cbar_kws={'label': 'Feature SVM weight'})
        feats_sorted = feats_sorted.iloc[-20:]
        #perf per feat
        #seaborn.heatmap(pd.DataFrame(feats_sorted).T)
        plt.suptitle("Discriminatory features across samples and their contribution to the classifier")
        ax1.set_xlabel('Features')
        ax1.set_ylabel('Samples (colored by class)')
        ax2.set_ylabel('SVM C-parameter')
        #ax3_labels = [item.get_text() for item in ax3.get_xticklabels()]
        #for i in range(2,5):
        #    ax3_labels [i]= ""
        #ax3.set_yticklabels(ax3_labels)
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.tight_layout()

        #adjust label size to number of samples
        for label in (ax1.get_yticklabels()):
                if len(y) > 10:
                    label.set_fontsize(12)
                if len(y) > 30:
                    label.set_fontsize(8)
                if len(y) > 50:
                    label.set_fontsize(5)
                if len(y) > 70:
                    label.set_fontsize(3)
                if len(y) > 100:
                    label.set_fontsize(2)
                if len(y) > 200:
                    label.set_fontsize(1)
                if y.loc[str(phn_name2id.loc[label.get_text(), :].iloc[0])] == 1:
                    label.set_color('orange')
                else:
                    label.set_color('blue')
        plt.savefig("%s/%s_heatmap.png" %(self.model_out, pt_out), dpi = 300)
        plt.close()
