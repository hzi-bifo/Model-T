#import NCBI_tree_of_life
import pandas as ps
import ete2
import traitar.evaluation as evaluate 
import sys

def macro_accuracy(numbers):
    POS, NEG, FN, FP = numbers
    if POS >= 5 and NEG >= 5:
        conf = [NEG - FP, FP, FN, POS - FN]
        pos_rec = evaluate.evaluate.recall_pos_conf(conf)
        neg_rec = evaluate.evaluate.recall_neg_conf(conf)
        bacc = evaluate.evaluate.bacc(pos_rec, neg_rec)
        return bacc
    else: 
        return float('NaN') 

def nanmean(x):
    return x[ps.notnull(x)].mean()

def map_miscl2taxonomy(miscl_m_dir, tree_f, out_dir):
    """read in table of misclassified species and propagate the misclassified statistics to the ancestral taxa"""
    #read in the matrix of misclassified species statistics
    miscl_m = ps.read_csv("%s/misclassified_overall.tsv"%miscl_m_dir, sep = "\t", index_col = 0) 
    miscl_m.index = [str(i) for i in miscl_m.index]
    keys = ["FN", "FP", "NEG", "POS",  "overall"]
    ncbi_tree = ete2.Tree(tree_f)
    per_pt_dict = dict((i, ps.read_csv("%s/misclassified_per-pt_%s.tsv"%(miscl_m_dir, i), sep = "\t", index_col = 0)) for i in ["FN", "FP", "NEG", "POS"])
    for i in per_pt_dict:
        per_pt_dict[i].index = [str(j) for j in per_pt_dict[i].index]
    per_pt_dict["overall"] = miscl_m
    gc = ps.Series(ps.np.ones(shape = (per_pt_dict["overall"].shape[0])), name = "gc")
    gc.index = per_pt_dict["overall"].index
    per_pt_dict["overall"]  = ps.concat([per_pt_dict["overall"], gc], axis = 1)
    per_pt_extd_dict = {}
    ncbi_ids = [j.name for j in ncbi_tree.traverse()]
    miscl_in_ncbi = [i in ncbi_ids for i in miscl_m.index]
    for i in per_pt_dict:
        per_pt_extd_dict[i] = ps.DataFrame(ps.np.zeros(shape = (len([j for j in ncbi_tree.traverse()]), len(per_pt_dict[i].columns))))
        per_pt_extd_dict[i].index = [j.name for j in ncbi_tree.traverse()]
        per_pt_extd_dict[i].columns = per_pt_dict[i].columns
        per_pt_extd_dict[i].loc[miscl_m.index[miscl_in_ncbi], :] = per_pt_dict[i]
    #propagate the number of samples / misclassified samples up to the root
    for n in ncbi_tree.iter_leaves(): 
        n_leave = n
        while n.name != "NoName":
            n = n.up
            for i in keys:
                per_pt_extd_dict[i].loc[n.name, ] += per_pt_extd_dict[i].loc[n_leave.name, ]

    #compute accuracy for single phenotypes 
    bacc_per_pt = ps.DataFrame(ps.np.zeros_like(per_pt_extd_dict["FN"]))
    bacc_per_pt.index = per_pt_extd_dict["overall"].index 
    bacc_per_pt.columns = per_pt_extd_dict["FN"].columns 
    for i in per_pt_extd_dict["overall"].index:
        for j in per_pt_extd_dict["FN"].columns:
            FP = per_pt_extd_dict["FP"].loc[i, j] 
            FN = per_pt_extd_dict["FN"].loc[i, j] 
            N = per_pt_extd_dict["NEG"].loc[i, j] 
            P = per_pt_extd_dict["POS"].loc[i, j] 
            TP = P - FN
            TN = N - FP
            if N >= 3 and P >= 3:
                macc = evaluate.evaluate.get_performance(ps.np.array([TN, FP, FN, TP])).loc["macro accuracy"]
            else:
                macc = float('nan')
            bacc_per_pt.loc[i, j] = macc
    print per_pt_extd_dict["overall"].columns
    bacc_per_pt.to_csv("%s/bacc_per_pt.txt" % out_dir, sep = "\t")
    perf = per_pt_extd_dict["overall"].apply(lambda x: evaluate.evaluate.get_performance(ps.np.array([x.loc["negative_samples"] - x.loc["false_positives"], x.loc["false_positives"], x.loc["false_negatives"], x.loc["positve_samples"] - x.loc["false_negatives"]])), axis = 1)
    per_pt_extd_dict["overall"] = ps.concat([per_pt_extd_dict["overall"], perf], axis = 1)
    #produce results files
    for i in keys:
        per_pt_extd_dict[i].to_csv("%s/%s_extd.tsv"%(out_dir, i), sep = "\t")
    #compute macro accuracy if there is at least 10 samples in positive and negative class for this taxon
    macro_df =  ps.DataFrame(ps.np.zeros(shape = (len([j for j in ncbi_tree.traverse()]), len(per_pt_dict["FN"].columns))))
    macro_df.index = per_pt_extd_dict["FN"].index
    macro_df.columns = per_pt_extd_dict["FN"].columns
    for i in range(len(per_pt_extd_dict["FN"].index)):
        for j in range(len(per_pt_extd_dict["FN"].columns)):
            #print evaluate.evaluate.get_performance( [per_pt_extd_dict[k].iloc[i, j] for k in ["POS", "NEG", "FN", "FP"]])
            TP = per_pt_extd_dict["POS"].iloc[i, j]  - per_pt_extd_dict["FN"].iloc[i, j]
            TN = per_pt_extd_dict["NEG"].iloc[i, j]  - per_pt_extd_dict["FP"].iloc[i, j]
            FP = per_pt_extd_dict["FP"].iloc[i, j]
            FN = per_pt_extd_dict["FN"].iloc[i, j]
            #print TP, TN, FP, FN
            #print evaluate.evaluate.get_performance([TN, FP, FN, TP]).loc["macro accuracy"]

            macro_df.iloc[i, j] = evaluate.evaluate.get_performance([TN, FP, FN, TP]).loc["macro accuracy"]

    macro_df.to_csv("%s/macro_acc_per_pt.tsv"%out_dir, sep = "\t")
    macro_df_mean = macro_df.apply(nanmean, axis = 1)
    macro_df_mean.index = macro_df.index
    macro_df_mean.columns = "bacc" 
    macro_df_mean.to_frame().to_csv("%s/macro_acc_per_sample.tsv"%out_dir, sep = "\t")



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='propagate misclassified samples to the ancestors')
    parser.add_argument("miscl_in_dir", help='the input directory with the extant misclassified statistics ')
    parser.add_argument("tree_f",help='tree with sample taxa')
    parser.add_argument("out_dir",help='the output directory')
    args = parser.parse_args()
    map_miscl2taxonomy(args.miscl_in_dir, args.tree_f, args.out_dir)

    #map_miscl2taxonomy("./", "./")
        

