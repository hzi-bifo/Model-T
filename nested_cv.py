#params:
#indices for target variables (phenotypes)
#out
#import liblinear as ll
#from liblinear import *
import numpy as np
import pandas
import sklearn.svm as svm
from joblib import Parallel, delayed, load, dump
#TODO use the sklearn version of joblib instead to remove the direct joblib dependency
from operator import itemgetter
import math
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
    return  array_rs_cs, cs


def write_features( feats, model_out, pt_out,):
    f = open("%s/%s_majority_features.txt"%(model_out,pt_out), "w")
    pfam_f = open("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfam_pts_names_nl_desc.txt", "r")
    id2pf = {}
    ls = pfam_f.readlines()
    for i in range(len(ls)):
        elems = ls[i].strip().split("\t")
        id2pf[int(elems[0])-1] = (elems[1],elems[2])
    for feat in feats:
        f.write("%s\t%s\n"%(id2pf[int(feat)][0],id2pf[int(feat)][1]))
    f.close()

def recall_pos(y,y_pred):
    """compute recall of the positive class"""
    return sum(y[y == 1] == y_pred[y==1])/float(sum(y==+1))

def recall_neg(y, y_pred):
    """compute recall of the negative class"""
    return sum(y[y == -1] == y_pred[y==-1])/float(sum(y==-1))

def bacc(pos_acc, neg_acc):
    """compute balanced accuracy"""
    return float(pos_acc + neg_acc)/2

def cv(x,y,params , C, fold):
    """train model on given fold and return the predicted labels for the left out samples"""
    #print ">>>>>>>>>>>>>>>> fold", C
    neg_fold = np.ones(len(y),dtype='bool')
    for k in fold:
        neg_fold[k] = False
    x_bin_co = x[neg_fold,]
    y_bin_co = y[neg_fold]
    predictor = svm.LinearSVC(C=C)
    predictor.set_params(**params)
    predictor.fit(x_bin_co,y_bin_co)
    #print predictor.get_params()
    #print predictor.transform(x).shape
    #test on left out fold
    preds = predictor.predict(x[fold])
    return preds

def setup_folds(x, y, cv):
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
        folds.append(range(i*(fold_l+1), i*(fold_l+1)+fold_l+1))
    n_e = fold_l_e * (fold_l+1)
    for i in range(cv - fold_l_e):
        folds.append(range(n_e + fold_l*i, n_e + fold_l*i + fold_l))
    return folds

def outer_cv(x, y, params, c_params, cv_outer, cv_inner, n_jobs):
    """do a cross validation with cv_outer folds
    if only cv_outer is given do a cross validation for each value of the paramter C given
    if cv_inner is given, do a nested cross validation optimizing the paramter C in the inner loop"""
    #do n fold cross validation
    #number of elements in each invidual fold
    if not cv_inner is None:
        ofolds = setup_folds(x, y, cv_outer)
        ocv_preds = []
        for fold in ofolds:
            neg_fold = np.ones(len(y),dtype='bool')
            for k in fold:
                neg_fold[k] = False
            x_train =  x[neg_fold,]
            y_train =  y[neg_fold]
            #y_test =  y[fold,]
            x_test =  x[fold]
            ifolds = setup_folds(x_train, y_train, cv_inner)
            all_preds = np.zeros(shape=[len(y_train), len(c_params)])
            #do inner cross validation
            preds = []
            for pred in Parallel(n_jobs=n_jobs)(delayed(cv)(x_train, y_train, params, c_params[j], ifolds[i])
                    for  (i,j) in ((i,j) for i in range(len(ifolds)) for j in range(len(c_params)) )):
                preds += list(pred)
            all_preds =  np.array(preds)
            all_preds_r = np.resize(all_preds, (len(y_train), len(c_params)))
            #determine C param with the best accuracy
            baccs = [bacc(recall_pos(y_train, all_preds_r[:,j]), recall_neg(y_train, all_preds_r[:,j])) for j in range(len(c_params))]
            c_opt = c_params[np.argmax(np.array(baccs))]
            predictor = svm.LinearSVC(C=c_opt)
            predictor.set_params(**params)
            predictor.fit(x_train, y_train)
            p = list(predictor.predict(x_test))
            #print "current preds", p
            ocv_preds +=p
        return ocv_preds

    folds = setup_folds(x, y, cv_outer)
    all_preds = np.zeros(shape=[len(y), len(c_params)])
    #TODO account for the case when folds exceeds number of samples
    #TODO think about assigning the folds randomly to shuffle the input data
    for j in range(len(c_params)):
        preds = []
        for pred in Parallel(n_jobs=n_jobs)(delayed(cv)(x, y, params, c_params[j], folds[i])
                for  i in range(len(folds))):
            preds += list(pred)
        all_preds[:,j] = preds
    return all_preds


def majority_feat_sel(x, y, all_preds, params, c_params, k, model_out, pt_out):
    """determine the features occuring in the majority of the k best models"""
    #determine the k best classifiers
    baccs = [bacc(recall_pos(y, all_preds[:,j]), recall_neg(y, all_preds[:,j])) for j in range(len(c_params))]
    recps = [recall_pos(y, all_preds[:,j]) for j in range(len(c_params))]
    recns = [recall_neg(y, all_preds[:,j]) for j in range(len(c_params))]
    baccs_s = sorted(((baccs[i], recps[i], recns[i], c_params[i]) for i in range(len(c_params))), key=itemgetter(0), reverse=True )
    models = np.zeros(shape=(x.shape[1], k))
    predictors = []
    for i in range(k):
        predictor = svm.LinearSVC(C=baccs_s[i][3])
        predictor.set_params(**params)
        predictor.fit(x, y)
        predictors.append(predictor)
        #save the model
        models[:,i] = predictor.coef_
    feats = []
    #determine the majority features
    for i in range(models.shape[0]):
        if sum(models[i,:] > 0) >= math.ceil(k/2.0):
            feats.append(str(i))

    baccs_s_np = np.array(baccs_s)[0:k,0:3].T
    names = ['bacc', "pos_rec", "neg_rec"]
    baccs_s_np_p = pandas.DataFrame(baccs_s_np).rename(dict((i,names[i]) for i in range(3)))
    pandas.DataFrame(baccs_s_np_p).to_csv("%s/perf_%s.txt"%(model_out,pt_out), sep="\t", header=c_params[0:k], float_format='%.3f')
    #write majority features to disk
    write_features(feats, model_out,pt_out)
    #write coefficient matrix to disk
    pandas.DataFrame(models).to_csv("%s/feats_%s.txt"%(model_out,pt_out), sep="\t", header=c_params[0:k])
    #pickle the predictors
    dump(predictors, '%s/pickled/predictors_%s.pkl'%(model_out,pt_out))

