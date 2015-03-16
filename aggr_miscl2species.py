import os.path
import pandas as ps
def collect_miscl(miscl_dir, pt_matrix_f, pt1, pt2, pt_start, out_f):
    ncbi2miscl = {}
    #read in pt matrix
    pt_m = ps.read_csv(pt_matrix_f, sep = '\t', index_col = 0, header = None)
    pt_m.columns = range(0, pt_start)
    print pt_m.columns
    miscl_m = ps.DataFrame(ps.np.zeros(shape = (234, 4)))
    miscl_m.index = pt_m.index[0:234]
    miscl_m.columns = ["positve_samples", "negative_samples", "false_positives", "false_negatives"]
    pt_has_model = [] 
    for i in range(pt1, pt2 + 1):
        #read in misclassified samples
        print "%s/%s_miscl.txt"%(miscl_dir, i)
        if os.path.exists("%s/%s_miscl.txt"%(miscl_dir, i)):
            pt_has_model.append(i)
            try:
                miscl = ps.read_csv("%s/%s_miscl.txt"%(miscl_dir, i), sep = "\t", header = None, index_col = 0)
            except ValueError:
                print "phenotype", i, "has no misclassified samples"
        else: 
            continue
        for j in miscl.index:
            if miscl.loc[j, 1] == 1.0:
                miscl_m.loc[j, "false_positives"] += 1 
            else: 
                miscl_m.loc[j, "false_negatives"] += 1
    miscl_m.loc[:, ["positve_samples", "negative_samples"]] = pt_m.loc[miscl_m.index, pt_has_model].apply(lambda x: ps.Series([(x.astype("string") == "1").sum(), (x.astype("string") == "0").sum()], index = ["positve_samples", "negative_samples"]), axis = 1) 
        
    #compute bacc, precision 
    bacc = miscl_m.apply(lambda x: (1 - (x.iloc[2]/float(x.iloc[0]) + x.iloc[3]/float(x.iloc[1]))/2), axis = 1)
    bacc.name = "bacc"
    print bacc
    joined = ps.concat([miscl_m, bacc], axis = 1)
    joined.to_csv(out_f, sep = "\t", index_label = False)


                

if __name__ == "__main__":
    collect_miscl(".", "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus/pfams_pts_counts_named.tsv", 8682, 8774, 8775, "misclassified.tsv")
    #collect_miscl(".", "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_bin_named.tsv", 8476, 8568, 8569, "misclassified.tsv")






