import pandas as pd
from sklearn import tree
import nested_cv

def base_classifier(traitar_model, phenotype_feature_table, features, phenotype, out, do_normalization): 
    """get base classifier for each feature"""
    model = pd.read_csv(traitar_model, sep = "\t", index_col = 0)
    sel_feats = model.index
    table = pd.read_csv(phenotype_feature_table, sep = "\t", index_col = 0)
    feats = pd.read_csv(features, sep = "\t", index_col = 0).index
    #target
    pt_notnull = pd.notnull(table.loc[:, phenotype])
    y_p = table.loc[:, phenotype].loc[pt_notnull,]
    y_p[y_p == 0] = -1
    #features
    x_p = table.loc[:, feats].loc[pt_notnull,]
    if do_normalization:
        scaler = preprocessing.StandardScaler(with_mean = True, with_std = True).fit(x_p)
        x_p = pd.DataFrame(data = scaler.transform(array), index = array.index, columns = array.columns)
    #train decision stump
    preds = [tree.DecisionTreeClassifier(max_depth = 1).fit(pd.DataFrame(x_p.loc[:, i]), y_p).predict(pd.DataFrame(x_p.loc[:, i])) for i in sel_feats] 
    conf_per_feat = pd.DataFrame([nested_cv.nested_cv.confusion_m(y_p, pd.Series(p, index = y_p.index).T) for p in preds ])
    conf_per_feat.index = sel_feats
    conf_per_feat.columns = ["TP", "FP", "FN", "TN"]
    #get macro accuracy
    bacc = conf_per_feat.apply(lambda x: nested_cv.nested_cv.bacc(nested_cv.nested_cv.recall_pos_conf(x), nested_cv.nested_cv.recall_neg_conf(x)), axis = 1)
    perf_per_feat = pd.concat([conf_per_feat, bacc], 1)
    perf_per_feat.columns = ["TP", "FP", "FN", "TN", "MACC"]
    feat_df = pd.concat([model.drop("cor", axis = 1, inplace = False), perf_per_feat], axis = 1)
    feat_df.sort(columns = ["MACC"], ascending = False).to_csv(out, float_format='%.3f',  sep = "\t")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("replace simple correlation with base classifier performanced per feature")
    parser.add_argument("traitar_model", help = "traitar model file")
    parser.add_argument("phenotype_feature_table", help = "feature table that was used to train the model")
    parser.add_argument("features", help = "features in the model")
    parser.add_argument("phenotype", help = "phenotype used")
    parser.add_argument("out", help = "out file")
    parser.add_argument("--do_normalization", action = 'store_true', help = "use in case of continuous inputs")
    args = parser.parse_args()
    base_classifier(**vars(args))
