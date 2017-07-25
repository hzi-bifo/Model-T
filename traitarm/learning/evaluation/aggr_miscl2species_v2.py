import os.path
import pandas as pd
def collect_miscl(miscl_dir, sample_ids_f, pt_mapping_f, evaluation_f, out_dir, are_phenotype_ids = False, allow_missing_samples = False):
    ncbi2miscl = {}
    #read in pt mapping 
    pt_id2acc = pd.read_csv(pt_mapping_f, sep = "\t", index_col = 0)
    pt_id2acc.index = pt_id2acc.index.values.astype('string')
    #read in sample ids
    sample_ids = pd.read_csv(sample_ids_f, sep = "\t", index_col = 0, header = None)
    #read in gold standard evaluation
    pt_m = pd.read_csv(evaluation_f, sep = "\t", index_col = 0)
    if are_phenotype_ids:
        pt_m.columns = pt_id2acc.loc[pt_m.columns, "accession"]
        print pt_m.columns
    #restrict to samples in sample ids 
    pt_m = pt_m.loc[sample_ids.index, :]
    #complete prediction matrix
    miscl_per_pt_FP, miscl_per_pt_FN, miscl_per_pt_POS, miscl_per_pt_NEG = [pd.DataFrame(pd.np.zeros( shape = (sample_ids.shape[0], pt_id2acc.shape[0])), index = sample_ids.index, columns = pt_id2acc.loc[:, "accession"]) for i in range(4)]
    miscl_m = pd.DataFrame(pd.np.zeros(shape = (sample_ids.shape[0], 4)))
    miscl_m.index = sample_ids.index 
    miscl_m.columns = ["positve_samples", "negative_samples", "false_positives", "false_negatives"]
    for i in pt_id2acc.loc[:, "accession"]: 
        if os.path.exists("%s/%s.txt"%(miscl_dir, i)):
            pos = pt_m.index[pt_m.loc[:, i].astype("string") == "+"]
            neg = pt_m.index[pt_m.loc[:, i].astype("string") == "-"]
            miscl_per_pt_POS.loc[pos, i] = 1
            miscl_per_pt_NEG.loc[neg, i] = 1
            miscl = pd.read_csv("%s/%s.txt"%(miscl_dir, i), sep = "\t", index_col = 0)
            try:
                miscl = pd.read_csv("%s/%s.txt"%(miscl_dir, i), sep = "\t", index_col = 0)
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
                    raise KeyError
            else:
                if miscl.loc[j, "gold standard"] == 1.0:
                    miscl_per_pt_FN.loc[j, i] = 1
                    miscl_m.loc[j, "false_negatives"] += 1 
                else: 
                    #miscl_per_pt_FP.loc[j, "gold standard"] = 0.0
                    miscl_per_pt_FP.loc[j, i] = 1
                    miscl_m.loc[j, "false_positives"] += 1
    miscl_m.loc[:, ["positve_samples", "negative_samples"]] = pt_m.loc[miscl_m.index, ].apply(lambda x: pd.Series([(x.astype("string") == "1").sum(), (x.astype("string") == "0").sum()], index = ["positve_samples", "negative_samples"]), axis = 1) 
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
    parser.add_argument("sample_ids",help='sample ids')
    parser.add_argument("evaluation_f",help='evaluation_f')
    parser.add_argument("pt_mapping_f", help='the output directory')
    parser.add_argument("--are_pt_ids", action = "store_true", help = "set if gold standard is indexed via phenotype ids rather than accessions")
    parser.add_argument('--allow_missing_samples', action = "store_true", help = "set if samples ids only contains a subset of the overall samples")
    args = parser.parse_args()
    collect_miscl(args.miscl_in_dir, args.sample_ids, args.pt_mapping_f, args.evaluation_f, args.out_dir, args.are_pt_ids, args.allow_missing_samples)
