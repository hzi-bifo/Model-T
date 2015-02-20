"""compute information theoretical measures for the features selected from the svm"""
import mutual_info
import pandas as ps

#parameters: phenotype range, feature directory, pfam matrix, ncbi id to ancestry table,  output table


#output table: 1) majority features extended with extra columns for MI, CMI, CWMI
#              2) table with mean MI,..., for each pt
#              3) overall mean MI, CMI, CWMI
mi = mutual_info.MutualInfo()
def get_it_measures(feat_v, pt_v, tax_v):
    """compute the different mi measures"""
    MI = mi.MutualInfo.MI(feat_v, pt_v)
    CMI = mi.CMI(feat_v, pt_v, tax_v)
    CWMI = mi.CWMI(feat_v, pt_v, tax_v)
    return ps.np.array([MI, CMI, CWMI])

def mi(pt1, pt2, feat_dir, pfam_f, ncbi2anc_f, out_table, tax_level):
    pfam_m = ps.read_csv(pfam_f, sep = "\t", index = 0)
    ncbi2anc = ps.read_csv(ncbi2anc_f, sep = "\t", index = 0)
    #matrix for pt summary 
    pt_mi = ps.DataFrame(ps.np.zeros(shape = (pt2 - pt1, 3)))
    pt_mi.index = range(pt1, p2 + 1)
    for i in range(pt1, pt2 + 1):
        #read in features 
        feats = ps.read_csv("%s/%s_majority_features.txt", sep = "\t", index_col = 0)
        feat_mi = ps.DataFrame(ps.np.zeros(shape = (feats.shape[0], 3)))
        feat_mi.index = feats.index
        for feat in feats.index:
            #restrict the pf vector and the ncbi ancestry vector to non-missing pts
            pt_index = pfam_f.index[pfam_f.loc[:, pt]  == "?"]
            feat_mi.loc[feat, ] = get_it_measures(pfam_m.loc[pt_index, feat], pfam_m.loc[ncbi2anc_f.loc[:, tax_level]
        #save the feat mi matrix to disk
        ps.concat(feats, feat_mi, axis = 0).to_csv("%s/%s_majority_features.txt", sep = "\t")
        #summarize mi measures for that phenotype and fill respective matrix row
        pt_mi.loc[i, ] = feat_mi.apply(ps.np.mean, axis = 1)
    
             



