import os.path
import pandas as ps
def collect_miscl(miscl_dir, pt_matrix_f, pt1, pt2, out_dir, mask_pts):
    ncbi2miscl = {}
    #read in pt matrix
    pt_m = ps.read_csv(pt_matrix_f, sep = '\t', index_col = 0, header = None)
    pt_m.columns = range(0, pt_m.shape[1])
    #complete prediction matrix
    miscl_per_pt_FP, miscl_per_pt_FN, miscl_per_pt_POS, miscl_per_pt_NEG  = [ps.DataFrame(ps.np.zeros(shape = (234, pt2 + 1 - pt1)) , index = pt_m.index[0:234], columns = range(pt1, pt2 + 1)) for i in range(4)]
    #matrix of #pos samples #neg -samples #FP #FN
    miscl_m = ps.DataFrame(ps.np.zeros(shape = (234, 4)))
    miscl_m.index = pt_m.index[0:234]
    miscl_m.columns = ["positve_samples", "negative_samples", "false_positives", "false_negatives"]
    pt_has_model = [] 
    for i in range(pt1, pt2 + 1):
        #read in misclassified samples
        #print "%s/%s_miscl.txt"%(miscl_dir, i)
        if os.path.exists("%s/%s_miscl.txt"%(miscl_dir, i)) and i not in mask_pts:
            pt_has_model.append(i)
            pos = pt_m.index[pt_m.loc[:, i].astype("string") == "1"]
            neg = pt_m.index[pt_m.loc[:, i].astype("string") == "0"]
            miscl_per_pt_POS.loc[pos, i] = 1
            miscl_per_pt_NEG.loc[neg, i] = 1
            try:
                miscl = ps.read_csv("%s/%s_miscl.txt"%(miscl_dir, i), sep = "\t", header = None, index_col = 0)
            except ValueError:
                print "phenotype", i, "has no misclassified samples"
                continue
        else: 
            continue


        for j in miscl.index:
            if miscl.loc[j, 1] == 1.0:
                miscl_per_pt_FN.loc[j, i] = 1
                #miscl_per_pt_NEG.loc[j, i] = 1
                miscl_m.loc[j, "false_negatives"] += 1 
            else: 
                miscl_per_pt_FP.loc[j, i] = 1
                #miscl_per_pt_POS.loc[j, i] = 1
                miscl_m.loc[j, "false_positives"] += 1
    miscl_m.loc[:, ["positve_samples", "negative_samples"]] = pt_m.loc[miscl_m.index, pt_has_model].apply(lambda x: ps.Series([(x.astype("string") == "1").sum(), (x.astype("string") == "0").sum()], index = ["positve_samples", "negative_samples"]), axis = 1) 
        
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
    parser.add_argument('pt_start', help = "the index of the first pt column to consider", type = int)
    parser.add_argument('pt_stop',  help='the index of the last pt column to consider', type = int)
    parser.add_argument('-m', '--mask_pts', default = [],  help="a comma separated list of phenotype ids that should be discarded")
    args = parser.parse_args()
    if not args.mask_pts == []:
        args.mask_pts = [int(i) for i in args.mask_pts.split(",")]
    collect_miscl(args.miscl_in_dir, args.pf_pt_matrix, args.pt_start, args.pt_stop, args.out_dir, args.mask_pts)
    #collect_miscl(".", "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts_named.tsv", 8476, 8568, 8569, "./")
    #collect_miscl(".", "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_bin_named.tsv", 8476, 8568, 8569, "misclassified.tsv")






