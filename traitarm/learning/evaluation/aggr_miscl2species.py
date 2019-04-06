import os.path
import pandas as ps
def collect_miscl(miscl_dir, pf_pt_matrix_f, pt_id2acc, out_dir, sample_ids_f, allow_missing_samples):
    ncbi2miscl = {}
    #read in pt matrix
    pt_m = ps.read_csv(pf_pt_matrix_f, sep = '\t', index_col = 0)
    pt_id2acc = ps.read_csv(pt_id2acc, sep = "\t", index_col = 0)
    pt_id2acc.index = pt_id2acc.index.values.astype('string')
    #read in sample ids
    sample_ids = ps.read_csv(sample_ids_f, header = None, index_col = 0, sep = "\t")
    #reduce to samples that are actually being used
    pt_m = pt_m.loc[sample_ids.index,]
    #complete prediction matrix
    miscl_per_pt_FP, miscl_per_pt_FN, miscl_per_pt_POS, miscl_per_pt_NEG = [ps.DataFrame(ps.np.zeros( shape = (pt_m.shape[0], pt_id2acc.shape[0])), index = pt_m.index, columns = pt_id2acc.index) for i in range(4)]
    #matrix of #pos samples #neg -samples #FP #FN
    miscl_m = ps.DataFrame(ps.np.zeros(shape = (pt_m.shape[0], 4)))
    miscl_m.index = pt_m.index 
    miscl_m.columns = ["positve_samples", "negative_samples", "false_positives", "false_negatives"]
    pt_has_model = [] 
    for i in pt_id2acc.index:
        #read in misclassified samples
        #print "%s/%s_miscl.txt"%(miscl_dir, i)
        if os.path.exists("%s/%s_miscl.txt"%(miscl_dir, i)):
            #correct index when using GIDEON indexed phenotype matrix
            #pos = pt_m.index[pt_m.loc[:, str(int(i) + 1)].astype("string") == "1"]
            #neg = pt_m.index[pt_m.loc[:, str(int(i) + 1)].astype("string") == "0"]
            pos = pt_m.index[pt_m.loc[:, i].astype("string") == "1"]
            neg = pt_m.index[pt_m.loc[:, i].astype("string") == "0"]
            miscl_per_pt_POS.loc[pos, i] = 1
            miscl_per_pt_NEG.loc[neg, i] = 1
            try:
                miscl = ps.read_csv("%s/%s_miscl.txt"%(miscl_dir, i), sep = "\t", index_col = 0)
            except ValueError:
                print "phenotype", i, "has no misclassified samples"
                continue
        else: 
            continue


        for j in miscl.index:
            if j not in miscl_m.index:
                if allow_missing_samples:
                    pass
                else:
                    print j
                    print miscl.index
                    raise KeyError
            else:
                #mark as false positive if prediction is positive ("1") and therefore the ground truth is negative
                if miscl.loc[j, "1"] == 1.0:
                    miscl_per_pt_FP.loc[j, i] = 1
                    #miscl_per_pt_NEG.loc[j, i] = 1
                    miscl_m.loc[j, "false_positives"] += 1 
                else: 
                    miscl_per_pt_FN.loc[j, i] = 1
                    #miscl_per_pt_POS.loc[j, i] = 1
                    miscl_m.loc[j, "false_negatives"] += 1
    #correct index when using gideon indexed phenotype matrix
    #miscl_m.loc[:, ["positve_samples", "negative_samples"]] = pt_m.loc[miscl_m.index, (pt_id2acc.index.values.astype(int) + 1).astype('str')].apply(lambda x: ps.Series([(x.astype("string") == "1").sum(), (x.astype("string") == "0").sum()], index = ["positve_samples", "negative_samples"]), axis = 1) 
    miscl_m.loc[:, ["positve_samples", "negative_samples"]] = pt_m.loc[miscl_m.index, pt_id2acc.index].apply(lambda x: ps.Series([(x.astype("string") == "1").sum(), (x.astype("string") == "0").sum()], index = ["positve_samples", "negative_samples"]), axis = 1) 
        
    #compute bacc, precision 
    #bacc = miscl_m.apply(lambda x: (1 - (x.iloc[2]/float(x.iloc[0]) + x.iloc[3]/float(x.iloc[1]))/2), axis = 1)
    #bacc.name = "bacc"
    #joined = ps.concat([miscl_m, bacc], axis = 1)
    #joined.to_csv("%s/misclassified_overall.tsv"%out_dir, sep = "\t", index_label = False)
    #write single pt prediction results to disk
    miscl_m.to_csv("%s/misclassified_overall.tsv"%out_dir, sep = "\t",)
    miscl_per_pt_FP.to_csv("%s/misclassified_per-pt_FP.tsv"%out_dir, sep = "\t",)
    miscl_per_pt_FN.to_csv("%s/misclassified_per-pt_FN.tsv"%out_dir, sep = "\t",)
    miscl_per_pt_POS.to_csv("%s/misclassified_per-pt_POS.tsv"%out_dir, sep = "\t",)
    miscl_per_pt_NEG.to_csv("%s/misclassified_per-pt_NEG.tsv"%out_dir, sep = "\t",)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("combine the misclassified samples of different phenotypes into data matrices")
    parser.add_argument("miscl_in_dir", help='the input directory with the extant misclassified statistics ')
    parser.add_argument("out_dir",help='the output directory')
    parser.add_argument("pf_pt_matrix", help='combined matrix of pfams and phenotypes')
    parser.add_argument('pt2acc', help = "phenotype mapping file")
    parser.add_argument('sample_ids', help = "file with samples")
    parser.add_argument('--allow_missing_samples', action = "store_true", help = "set if samples ids only contains a subset of the overall samples")
    args = parser.parse_args()
    collect_miscl(args.miscl_in_dir, args.pf_pt_matrix, args.pt2acc, args.out_dir, args.sample_ids, args.allow_missing_samples)




