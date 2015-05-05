import pandas as ps
import sys

def feat_sel(weights_f, id2pf2desc_f, k = 5, strategy = "majority"):
    w = ps.read_csv(weights_f, sep = "\t", index_col = 0)
    id2pf2desc = ps.read_csv(id2pf2desc_f, sep = "\t", index_col = 0)
    if strategy == "majority": 
        sel_feats = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() > k/2, axis = 1) ]
        #print sel_feats
    elif strategy =="conservative":
        sel_feats = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() == k, axis = 1), ]
    elif strategy =="any":
        sel_feats = w.loc[w.apply(lambda x: (x.iloc[0:k][x > 0] > 0).sum() >= 1, axis = 1), ]
    else:
        raise IllegalArgumentException
    id2pf2desc.index = [int(i) - 1 for i in id2pf2desc.index]
    sel_feats[['Pfam_acc', 'Pfam_desc']] = id2pf2desc.loc[sel_feats.index,]
    cols = sel_feats.columns.tolist()
    cols = cols[-2:-1] + cols[:-2] + cols[-1:]
    sel_feats = sel_feats[cols]
    sel_feats.to_csv(sys.stdout, sep = "\t", float_format = '%.3f')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Extract important features from the SVM weight matrix and write to stdout")
    parser.add_argument("weight_f", help='SVM weight matrix, columns sorted by descreasing accuracy ')
    parser.add_argument("id2pf2desc_f",help='mapping from feature id to Pfam accession and description')
    parser.add_argument("-s", "--strategy", default = "majority", help='mapping from feature id to Pfam accession and description')
    parser.add_argument("-k", "--voters", default = 5, help='voters to be considered for the feature selection',type = int)
    args = parser.parse_args()
    feat_sel(args.weight_f,  args.id2pf2desc_f, k = args.voters ,strategy = args.strategy)
        
