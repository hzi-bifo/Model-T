import pandas as ps
def process(phypat_dir, phypat_ggl_dir, out_dir, strategy):
    a = ["POS", "NEG", "FP", "FN"]
    p = [ps.read_csv("%s/misclassified_per-pt_%s.tsv" % (phypat_dir, i), sep = "\t", index_col = 0).astype('int').astype('bool') for i in a]
    p_ggl = [ps.read_csv("%s/misclassified_per-pt_%s.tsv" % (phypat_ggl_dir, i), sep = "\t", index_col = 0).astype('int').astype('bool') for i in a]
    #adjust index of phypat
    for i in p:
        i.columns = i.columns.values.astype('int') + 206
    for i in p_ggl:
        i.columns = i.columns.values.astype("int")
    #get rid of phypat+GGL phenotypes not in phypat
    pts = []
    for i in p_ggl[0].columns:
        if p[0].sum().loc[i] != 0:
           pts.append(i) 
    print pts
    #combine matrix
    p_comb = [p[0], p[1], p[2] | p_ggl[2].loc[:, pts]  if strategy == "majority" else p[2] & p_ggl[2].loc[:, pts], p[3] & p_ggl[3].loc[:, pts]  if strategy == "majority" else p[3] | p_ggl[3].loc[:, pts]  ]
    for i in range(len(p)):
        p_comb[i].astype('float').to_csv("%s/misclassified_per-pt_%s.tsv" %(out_dir, a[i]), sep = "\t")
    overall = ps.concat([i.sum(axis = 1) for i in p_comb], axis = 1)
    overall.columns = ["positive_samples", "negative_samples", "false_positives", "false_negatives"]
    overall.to_csv("%s/misclassified_overall.tsv" % out_dir, sep = "\t")
    #compute overall matrix
    

    
#read in phypat and phypat+GGL predictions
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("combine the misclassified samples of phypat and phypat+GGL classifier")
    parser.add_argument("phypat_in_dir", help='phypat directory')
    parser.add_argument("phypat_ggl_in_dir", help='phypat+GGL directory')
    parser.add_argument("out_dir", help='out_dir')
    parser.add_argument("strategy", choices = ['majority', 'conservative'], help='out_dir')
    args = parser.parse_args()
    process(args.phypat_in_dir, args.phypat_ggl_in_dir, args.out_dir, args.strategy)

