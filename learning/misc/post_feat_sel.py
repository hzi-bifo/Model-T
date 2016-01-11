import pandas as ps
import sys
import get_cor_feats

def feat_sel(weights_f, id2pf2desc_f, k = 5, strategy = "majority", include_negative_class = False, pt_correlation = False, phypat_f = None, target = None):
    w = ps.read_csv(weights_f, sep = "\t", index_col = 0)
    id2pf2desc = ps.read_csv(id2pf2desc_f, sep = "\t", index_col = 0, header = None)
    if strategy == "majority": 
        sel_feats_pos = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() > k/2, axis = 1) ]
        if include_negative_class:
            sel_feats = ps.concat([sel_feats_pos, w.loc[w.apply(lambda x: (x.iloc[0:k][x < 0] < 0).sum() > k/2, axis = 1) ]])
    elif strategy =="conservative":
        sel_feats_pos = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() == k, axis = 1), ]
        if include_negative_class:
            sel_feats = ps.concat([sel_feats_pos, w.loc[w.apply(lambda x: (x.iloc[0:k][x < 0] < 0).sum() == k, axis = 1), ]])
    elif strategy =="any":
        sel_feats_pos = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() >= 1, axis = 1), ]
        if include_negative_class:
            sel_feats = ps.concat([sel_feats_pos, w.loc[w.apply(lambda x: (x.iloc[0:k][x < 0] < 0).sum() >= 1, axis = 1), ]])
    else:
        raise IllegalArgumentException
    if not include_negative_class:
        sel_feats = sel_feats_pos
    id2pf2desc.index = [int(i) - 1 for i in id2pf2desc.index]
    #print id2pf2desc
    sel_feats[['Pfam_acc', 'Pfam_desc']] = id2pf2desc.loc[sel_feats.index,]
    if include_negative_class:
        class_col = ps.Series(ps.np.zeros(sel_feats.shape[0]), dtype = "string")
        class_col.iloc[:sel_feats_pos.shape[0]] = ps.np.repeat("+", sel_feats_pos.shape[0])
        class_col.iloc[sel_feats_pos.shape[0]:] = ps.np.repeat("-", sel_feats.shape[0] - sel_feats_pos.shape[0])
        class_col.index = sel_feats.index
        class_col.name = "class"
    else:
        class_col =  ps.Series(ps.np.repeat("+", sel_feats.shape[0]), name = "class", index = sel_feats.index)
    sel_feats = ps.concat([sel_feats, ps.DataFrame(class_col)], axis = 1)
    cols = sel_feats.columns.tolist()
    cols = cols[-3:-2] + cols[-1:] +  cols[:-3] + cols[-2:-1] 
    sel_feats = sel_feats[cols]
    if pt_correlation:
        phypat = ps.read_csv(phypat_f, index_col = 0, sep = "\t", header = None)
        pfams = ps.read_csv(id2pf2desc_f, sep = "\t", header = None)
        pfams.index = pfams.loc[:, 1]
        target_id = pfams.loc[target, 0] - 1 
        phypat_red = phypat.loc[phypat.iloc[:, target_id] != "?",]
        #phypat_red = phypat
        phypat_red_pt = phypat_red.iloc[:, id2pf2desc.iloc[range(w.shape[0]), :].index.tolist() + [target_id]].astype('int')
        phypat_red_pt.columns = pfams.iloc[range(w.shape[0]) + [target_id], 1]
        phypat_red_pt = (phypat_red_pt > 0).astype('int')
        cor1 = get_cor_feats.get_selected_cor(target, sel_feats.iloc[:sel_feats_pos.shape[0], 0], phypat_red_pt, pfams)
        cor2 = get_cor_feats.get_selected_cor(target, sel_feats.iloc[sel_feats_pos.shape[0]:, 0], phypat_red_pt, pfams)
        cor = ps.concat([cor1, cor2], axis = 0)
        cor.columns = ["cor", "_"]
        cor = cor.loc[sel_feats.loc[:, "Pfam_acc"], ].iloc[:, 0]
        cor.index = sel_feats.index
        sel_feats = ps.concat([sel_feats, cor] , axis = 1)
        sel_feats = sel_feats.reindex(sel_feats.loc[:, "cor"].abs().order(ascending = False).index)

    sel_feats.to_csv(sys.stdout, sep = "\t", float_format = '%.3f')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Extract important features from the SVM weight matrix and write to stdout")
    parser.add_argument("weight_f", help='SVM weight matrix, columns sorted by descreasing accuracy ')
    parser.add_argument("id2pf2desc_f",help='mapping from feature id to Pfam accession and description')
    parser.add_argument("-s", "--strategy", default = "majority", help='mapping from feature id to Pfam accession and description')
    parser.add_argument("-k", "--voters", default = 5, help='voters to be considered for the feature selection',type = int)
    parser.add_argument("-n", "--include_negative_class", action = 'store_true') 
    parser.add_argument("-c", "--pt_correlation", action = 'store_true', help = "compute correlation with phenotype") 
    parser.add_argument("-p", "--phypat_f",  help = "phyletic pattern matrix file") 
    parser.add_argument("-t", "--target",  help = "target phenotype") 
    args = parser.parse_args()
    feat_sel(args.weight_f,  args.id2pf2desc_f, k = args.voters ,strategy = args.strategy, include_negative_class = args.include_negative_class, pt_correlation = args.pt_correlation, phypat_f = args.phypat_f, target = args.target)
        
