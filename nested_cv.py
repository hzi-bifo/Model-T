#params:
#indices for target variables (phenotypes)
#out
#import liblinear as ll
#from liblinear import *
import numpy as np
import pandas as ps
import sklearn.svm as svm
from joblib import Parallel, delayed, load, dump
#TODO use the sklearn version of joblib instead to remove the direct joblib dependency
from operator import itemgetter
import math
import random
import cv_rec_helper as crh 
import sys
import itertools
import sklearn.preprocessing as preprocessing

def normalize(array):
    #TODO check normalization procedure
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
    #scaler = preprocessing.StandardScaler(with_mean = False, with_std = True).fit(array)
    #transformed = ps.DataFrame(data = scaler.transform(array), index = array.index, columns = array.columns)
    #return transformed, scaler.std_
    return  array_rs_cs, cs


def get_pfam_names_and_descs(feats):
    pfam_f = open("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfam_pts_names_nl_desc.txt", "r")
    id2pf = {}
    ls = pfam_f.readlines()
    for i in range(len(ls)):
        elems = ls[i].strip().split("\t")
        #mapping counts from 1 
        id2pf[int(elems[0]) - 1] = (elems[1],elems[2])
    pfam_f.close()
    return id2pf

def recall_pos(y,y_pred):
    """compute recall of the positive class"""
    return (y[y == 1] == y_pred[y==1]).sum()/float((y==+1).sum())

def recall_neg(y, y_pred):
    """compute recall of the negative class"""
    return (y[y == -1] == y_pred[y==-1]).sum()/float((y==-1).sum())

def bacc(pos_acc, neg_acc):
    """compute balanced accuracy"""
    return float(pos_acc + neg_acc)/2

def cv(x_train, y_train, xp_train, yp_train, x_test, params , C, is_phypat_and_rec, perc_feats, no_classifier = 10):
    """train model on given training features and target and return the predicted labels for the left out samples"""
    #set y negative labels to -1
    y_train_t = y_train.copy() 
    y_train_t[y_train_t == 0] = -1
    yp_train_t = yp_train.copy() 
    yp_train_t[yp_train_t == 0] = -1
    predictor = svm.LinearSVC(C=C)
    predictor.set_params(**params)
    #check if we are in vanilla linear SVM and set no_classifiers to 1 if so or in subspace mode 
    if perc_feats == 1.0:
        no_classifier = 1
    all_preds = ps.DataFrame(np.zeros(shape = (x_test.shape[0], no_classifier)))
    for i in range(no_classifier):
        #select a subset of features for classification
        sample = sorted(random.sample(x_train.columns, int(math.floor(x_train.shape[1] * perc_feats))))
        #reduce feature space to the selected features
        x_train_sub = x_train.loc[:, sample]
        if is_phypat_and_rec:
            #reduce xp feature space to the selected features
            xp_train_sub = xp_train.loc[:, sample]
            predictor.fit(ps.concat([x_train_sub, xp_train_sub], axis = 0), ps.concat([y_train_t, yp_train_t], axis = 0))
        else: 
            predictor.fit(x_train_sub, y_train_t)
        
        all_preds.iloc[:, i]  = predictor.predict(x_test.loc[:, sample])
    #do majority vote to aggregate predictions
    aggr_preds = all_preds.apply(lambda x: 1 if sum(x == 1) > sum(x == -1) else -1, axis = 1).astype('int')
    return aggr_preds 

def get_training_and_testing_mapping(train_pt_edges, all_pt_edges):
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
            #all_pt_edges.remove(e2)
    test_edges = all_pt_edges.intersection(set(itertools.chain.from_iterable(train2all.values())))
    return train2all, test_edges

def join_edges(x_r, train2all, likelihood_params, parsimony_params):
    """aggregate edges in the full reconstruction matrix to fit the joint edges in train2all"""
    df = ps.DataFrame()
    for e in train2all:
        #aggregate rows by summing up
        if not likelihood_params is None:
            cs = x_r.loc[train2all[e], :].sum(axis=0)
            df = ps.concat([df, cs], axis = 1)
            #print "df shape", df.shape
        else: 
            raise Exception("parsimony case not yet implemented")    
    #make sure there are no entries in the summed up likelihood matrix other than 0 and 1
    if not likelihood_params is None:
        df[x_r==2] = 1
    #get df back into samples x features shape 
    df.columns = train2all.keys()
    return df.T


def get_rec_samples(x_r, y_r, yp_train, model_out, likelihood_params, parsimony_params):
    """compute the training samples by running parsimony/likelihood reconstruction for the pt for the training set without the held out samples"""
    #retrieve matrix like:
    #N1_N4 0 
    #N2_44 1
    if not parsimony_params is None:
        m = crh.reconstruct_pt_parsimony(yp_train, parsimony_params)
    else:
        m = crh.reconstruct_pt_likelihood(yp_train, model_out, likelihood_params)
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
    train2all, test_edges = get_training_and_testing_mapping(df2, df1)
    #print "mapping from edges only in the training reconstruction to edges in the full reconstruction", train2all
    #print "held out edges from get training and testing mapping", test_edges
    #get the corresponding samples in the full reconstruction matrix and aggregate them
    joint_x_r = join_edges(x_r, train2all, likelihood_params, parsimony_params)
    #add the y labels to those rows
    joint_y_r = m.loc[train2all, :].iloc[:, 0] 
    #add all the other training samples
    training_edges = all_pt_edges.intersection(train_pt_edges)
    #print x_r.loc[training_edges,:].shape
    #print joint_x_r.shape, joint_x_r.index
    x_r_train = ps.concat([x_r.loc[training_edges,:], joint_x_r], axis = 0)
    #print "x_r_train shape", x_r_train.shape
    #print "x_r_train index", x_r_train.index
    y_r_train = ps.concat([m.loc[training_edges, :].iloc[:,0], joint_y_r], axis = 0)
    y_r_train[y_r_train == 2] = 1
    #print "y_r_train shape", y_r_train.shape
    #replace negative labels with -1
    x_r_test = x_r.loc[test_edges,]
    y_r_test = y_r.loc[test_edges,]
    return x_r_train, y_r_train, x_r_test, y_r_test


def setup_folds(x, y, cv, is_rec_based, x_p, y_p, model_out, likelihood_params, parsimony_params, is_phypat_and_rec):
    """prepare all combinations of training and test folds, either for the standard phyletic pattern or the likelihood / parsimony reconstruction case or for the combined case"""
    pfolds = setup_folds_phypat(x_p, y_p, cv)
    if not is_rec_based:
        #standard phyletic pattern cv
        return pfolds
    else:
        folds = []
        for xp_train, yp_train, xp_test, _, _, _, _ in pfolds:
            x_train, y_train, x_test, y_test  = get_rec_samples(x, y, yp_train, model_out, likelihood_params, parsimony_params)
            folds.append((x_train, y_train, x_test, y_test,  xp_train, yp_train, xp_test))
    return folds

def setup_folds_phypat(x, y, cv):
    """divide the samples into cv folds and return a list of index lists"""
    print "y length", len(y)
    fold_l = len(y)/cv
    print "fold_l", fold_l
    #number of folds with one more element
    fold_l_e = len(y)%cv
    print "fold_l_e", fold_l_e
    #list of folds lists
    folds = []
    for i in range(fold_l_e):
        index = range(i*(fold_l+1), i*(fold_l+1)+fold_l+1)
        #print x.shape, "phyletic pattern shape before removing test samples"
        x_train = x.drop(x.index[index])
        #print x_train.shape, "phyletic pattern shape after removing test samples"
        x_test = x.iloc[index,:]
        y_train = y.drop(y.index[index])
        folds.append((x_train, y_train, x_test, None, x_train, y_train, x_test))
    n_e = fold_l_e * (fold_l+1)
    for i in range(cv - fold_l_e):
        index = range(n_e + fold_l*i, n_e + fold_l*i + fold_l)
        x_train = x.drop(x.index[index])
        x_test = x.iloc[index,:]
        y_train = y.drop(y.index[index])
        folds.append((x_train, y_train, x_test, None, x_train, y_train, x_test))
    return folds



def outer_cv(x, y, params, c_params, cv_outer, cv_inner, n_jobs, is_rec_based, x_p, y_p, model_out, likelihood_params, parsimony_params, is_phypat_and_rec, perc_feats):
    """do a cross validation with cv_outer folds
    if only cv_outer is given do a cross validation for each value of the paramter C given
    if cv_inner is given, do a nested cross validation optimizing the paramter C in the inner loop"""
    #do n fold cross validation
    #number of elements in each invidual fold
    if not cv_inner is None:
        ofolds = setup_folds(x, y, cv_outer, is_rec_based, x_p, y_p, model_out, likelihood_params, parsimony_params, is_phypat_and_rec)
        ocv_preds = []
        for x_train, y_train, x_test, y_test,  xp_train, yp_train, xp_test  in ofolds:
            ifolds = setup_folds(x_train, y_train, cv_inner,  is_rec_based, xp_train, yp_train, model_out, likelihood_params, parsimony_params, is_phypat_and_rec)
            all_preds = ps.DataFrame(np.zeros(shape=[len(yp_train), len(c_params)]))
            all_preds.columns = c_params
            all_preds.index = yp_train.index
            #do inner cross validation
            preds = []
            for pred in Parallel(n_jobs=n_jobs)(delayed(cv)(x_train_train, y_train_train, xp_train_train, yp_train_train, xp_train_test, params, c_params[j], is_phypat_and_rec, perc_feats)
                    for(x_train_train, y_train_train, x_train_test, y_train_test,  xp_train_train, yp_train_train, xp_train_test, j) 
                        in ((x_train_train, y_train_train, x_train_test, y_train_test,  xp_train_train, yp_train_train, xp_train_test, j) for x_train_train, y_train_train, x_train_test, y_train_test,  xp_train_train, yp_train_train, xp_train_test in ifolds for j in range(len(c_params)))):
                preds.append(list(pred))
            preds_ordered = [[] for i in range(len(c_params))]
            for i in range(len(preds)):
                preds_ordered[i % len(c_params)] += preds[i]
            all_preds = ps.DataFrame(np.array(preds_ordered), index = c_params, columns = yp_train.index).T
            #determine C param with the best accuracy
            #copy yp train and change the negative labels from 0 to -1
            yp_t_t = yp_train.copy()
            yp_t_t[yp_t_t == 0] = -1
            baccs = [bacc(recall_pos(yp_t_t, all_preds.iloc[:,j]), recall_neg(yp_t_t, all_preds.iloc[:,j])) for j in range(len(c_params))]
            print baccs, "baccs"
            c_opt = c_params[np.argmax(np.array(baccs))]
            print c_opt
            #use that C param to train a classifier to predict the left out samples in the outer cv
            p = cv(x_train, y_train, xp_train, yp_train, xp_test, params , c_opt, is_phypat_and_rec, perc_feats, no_classifier = 10)
            ocv_preds += list(p)
        return ocv_preds

    folds = setup_folds(x, y, cv_outer, is_rec_based, x_p, y_p, model_out, likelihood_params, parsimony_params,is_phypat_and_rec)
    #for i in folds:
    #    for j in i:
            #print j.shape
    all_preds = np.zeros(shape=[len(y_p), len(c_params)])
    #TODO account for the case when folds exceeds number of samples
    #TODO think about assigning the folds randomly to shuffle the input data
    for j in range(len(c_params)):
        preds = []
        for pred in Parallel(n_jobs=n_jobs)(delayed(cv)(x_train, y_train, xp_train, yp_train, xp_test, params, c_params[j], is_phypat_and_rec, perc_feats, no_classifier = 10)
                for  x_train, y_train, x_test, y_test,  xp_train, yp_train, xp_test in folds):
            preds += list(pred)
        all_preds[:,j] = preds
    return all_preds


def majority_feat_sel(x, y, x_p, y_p, all_preds, params, c_params, k, model_out, pt_out, is_phypat_and_rec, perc_feats, no_classifier = 10):
    """determine the features occuring in the majority of the k best models"""
    #determine the k best classifiers
    all_preds.index = y_p.index
    x.columns = x_p.columns
    y_p_t = y_p.copy()
    y_p_t[y_p_t == 0] = -1
    y_t = y.copy()
    y_t[y_t == 0] = -1
    baccs = [bacc(recall_pos(y_p_t, all_preds.iloc[:,j]), recall_neg(y_p_t, all_preds.iloc[:,j])) for j in range(len(c_params))]
    recps = [recall_pos(y_p_t, all_preds.iloc[:,j]) for j in range(len(c_params))]
    recns = [recall_neg(y_p_t, all_preds.iloc[:,j]) for j in range(len(c_params))]
    baccs_s = sorted(((baccs[i], recps[i], recns[i], c_params[i]) for i in range(len(c_params))), key=itemgetter(0), reverse=True )
    predictors = []
    #check if we are in vanilla linear SVM and set no_classifiers to 1 if so or in subspace mode 
    if perc_feats == 1.0:
        no_classifier = 1
    models = np.zeros(shape=(x.shape[1], k * no_classifier))
    for i in range(k):
        predictor = svm.LinearSVC(C=baccs_s[i][3])
        predictor.set_params(**params)
        sample = random.sample(x.columns, int(math.floor(x.shape[1] * perc_feats)))
        for l in range(no_classifier):
            x_sub = x.loc[:, sample]
            if is_phypat_and_rec:
                x_p_sub = x_p.loc[:, sample]
                predictor.fit(ps.concat([x_sub, x_p_sub], axis = 0), ps.concat([y_t, y_p_t], axis = 0))
            else:
                predictor.fit(x_sub, y_t)
                #save the model
            predictors.append(predictor)
            models[[np.array(sample) - 1], i * no_classifier + l] = predictor.coef_
    feats = []
    #determine the majority features
    for i in range(models.shape[0]):
        #print models.shape[1], 'number of classifiers'
        #print math.ceil(k/2.0), 'threshold'
        if sum(models[i,:] > 0) >= math.ceil(k/2.0):
            feats.append(i)
    rownames = [baccs_s[i][3] for i in range(k)] 
    colnames = ['bacc', "pos_rec", "neg_rec"]
    baccs_s_np = np.array(baccs_s)[0:k,0:3].T
    baccs_s_np_p = ps.DataFrame(baccs_s_np).rename(dict((i,colnames[i]) for i in range(3)))
    ps.DataFrame(baccs_s_np_p).to_csv("%s/%s_perf.txt"%(model_out,pt_out), sep="\t", float_format = '%.3f',  header=rownames)
    #write majority features with their weights to disk
    #print feats
    id2pf = ps.DataFrame(get_pfam_names_and_descs(feats))
    models_df = ps.DataFrame(models)
    for i in range(k):
        models_df[models_df.columns[i]] = models_df[models_df.columns[i]].map(lambda x: '%.3f' % x)
    feat_df = ps.concat([id2pf.loc[:, feats], models_df.T.loc[:, feats]], axis = 0).T
    rownames_extd = [str(baccs_s[i][3]) + "_" + str(l) for i in range(k) for l in range(no_classifier)] 
    feat_df.columns = ["Pfam_acc", "Pfam_desc"] + rownames_extd
    columns_out = ["Pfam_acc"] + rownames_extd + ["Pfam_desc"]
    feat_df.to_csv("%s/%s_majority_features+weights.txt"%(model_out,pt_out), columns = columns_out, float_format='%.3f',  sep = "\t")
    id2pf.index = ("Pfam_acc", "Pfam_desc")
    id2pf.loc[:, feats].T.to_csv("%s/%s_majority_features.txt"%(model_out,pt_out), float_format='%.3f',  sep = "\t")
    #write coefficient matrix to disk
    #put column names
    ps.DataFrame(models).to_csv("%s/%s_feats.txt"%(model_out,pt_out), sep="\t")
    #pickle the predictors
    dump(predictors, '%s/pickled/%s_predictors.pkl'%(model_out,pt_out))

