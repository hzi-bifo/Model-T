"""compute information theoretical measures for the features selected from the svm"""
import mutual_info
import pandas as ps
import sys
import getopt
import numpy as np
import os

#parameters: phenotype range, feature directory, pfam matrix, ncbi id to ancestry table,  output table


#output table: 1) majority features extended with extra columns for MI, CMI, CWMI
#              2) table with mean MI,..., for each pt
#              3) overall mean MI, CMI, CWMI
mutual_info = mutual_info.MutualInfo()

def get_factor(vector):
    tax2lev = dict([(np.unique(vector)[i], i) for i in range(len(np.unique(vector)))]) 
    factor = ps.Series(ps.np.zeros_like(vector))
    for i in range(vector.shape[0]):
        factor.iloc[i] = tax2lev[vector.iloc[i]]
    return factor

def get_cmi_measures(feat_v, pt_v, tax_v):
    """compute the different mi measures"""
    feat_f, pt_f, tax_f = [get_factor(i) for i in [feat_v, pt_v, tax_v]]
    CMI = mutual_info.CMI(feat_f, pt_f, tax_f)
    CWMI = mutual_info.CWMI(feat_f, pt_f, tax_f)
    return ps.np.array([CMI, CWMI])

def get_mi(feat_v, pt_v):
    """compute mi"""
    feat_f, pt_f = [get_factor(i) for i in [feat_v, pt_v]]
    MI = mutual_info.MI(feat_f, pt_f)
    return MI

def get_tax_factor(tax_m):
    """transform the taxonomy affiliation matrix into a factor per taxonomic level from 0, ...., n"""
    tax2lev = [dict([(np.unique(tax_m.iloc[:, j])[i], i) for i in range(len(np.unique(tax_m.iloc[:, j])))]) for j in range(tax_m.shape[1])]
    tax_factor = ps.DataFrame(ps.np.zeros_like(tax_m))
    tax_factor.index = tax_m.index
    tax_factor.columns = tax_m.columns
    for i in range(tax_m.shape[0]):
        for j in range(tax_m.shape[1]):
            tax_factor.iloc[i, j] = tax2lev[j][tax_m.iloc[i, j]]
    return tax_factor

def mi(pt1, pt2, feat_dir, pfam_f, ncbi2anc_f, out_table, tax_levels, last_feat_id):
    pfam_m = ps.read_csv(pfam_f, sep = "\t", index_col = 0, header = None)
    pfam_m.iloc[:, 0:last_feat_id] = (pfam_m.iloc[:, 0:last_feat_id] > 0).astype('int')
    ncbi2anc_text = ps.read_csv(ncbi2anc_f, sep = "\t", index_col = 0)
    ncbi2anc = get_tax_factor(ncbi2anc_text)
    #matrix for pt summary 
    pt_mi = ps.DataFrame(ps.np.zeros(shape = (pt2 - pt1 + 1, 2 * len(tax_levels) + 1)))
    pt_mi.index = range(pt1, pt2 + 1)
    pt_mi.columns = [ m + "_" + t for m in "CMI", "CMWI" for t in tax_levels] + ["MI"]
    pt_has_model = dict([(i, True) for i in range(pt1, pt2 + 1)])
    for i in range(pt1, pt2 + 1):
        #read in features 
        if not os.path.exists("%s/%s_majority_features.txt" % (feat_dir, i)):
            pt_has_model[i] = False 
            continue
        feats = ps.read_csv("%s/%s_majority_features.txt" % (feat_dir, i), sep = "\t", index_col = 0)
        feat_mi = ps.DataFrame(ps.np.zeros(shape = (feats.shape[0], 2 * len(tax_levels) + 1)))
        feat_mi.columns = [ m + "_" + t for m in "CMI", "CMWI" for t in tax_levels] + ["MI"]
        feat_mi.index = feats.index
        #restrict the pf vector and the ncbi ancestry vector to non-missing pts
        pt_index = pfam_m.index[~(pfam_m.iloc[:, i].astype("string")  == "?")]
        for feat in feats.index:
            feat_mi.loc[feat, "MI"] = get_mi(pfam_m.loc[pt_index, ].iloc[:, feat], pfam_m.loc[pt_index, ].iloc[:, i].astype('int'))
            for tax_level in tax_levels:
                tax_labels = [ m + "_" + tax_level for m in "CMI", "CMWI"]
                feat_mi.loc[feat, tax_labels] = get_cmi_measures(pfam_m.loc[pt_index, ].iloc[:, feat], pfam_m.loc[pt_index, ].iloc[:, i].astype('int'), ncbi2anc.loc[pt_index, tax_level])
        #save the feat mi matrix to disk
        ps.concat([feats, feat_mi], axis = 1).to_csv("%s/%s_majority_features_mi.txt"% (feat_dir, i), sep = "\t")
        #summarize mi measures for that phenotype and fill respective matrix row
        pt_mi.loc[i, ] = feat_mi.apply(ps.np.mean, axis = 0)
    #write summary table to disk
    pt_mi.to_csv(out_table, sep = "\t")
    #print summary to stdout
    print pt_mi.loc[[pt_has_model[i] for i in pt_mi.index],].apply(ps.np.mean)
    
             


if __name__=="__main__":
    #only testing
    if len(sys.argv) == 1:
        print """USAGE: python %s
-m <feat_dir> directory with the majority features 
-p <pt_range> range of the phentypes e.g. 8487-8490 
-f <pfam_f>  pfam / pt matrix 
-t <taxonomy_f> file with the taxnomoic affiliations of the samples in matrix form, one taxonomic level per column
-l <taxonomic level> one of genus, family, order, class, phylum
-i <id> of last feature e.g. 8475
-o <out_table> with mutual information summary over all phenotypes and all majority features
        """ % (sys.argv[0])
        sys.exit(2)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "m:p:f:t:l:i:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    pt1 = pt2 = feat_dir = pfam_f = ncbi2anc_f = out_table = tax_level = last_feat_id =  None
    for o, a in optlist:
        if o == "-m":
            feat_dir = a
        if o == "-p":
            pt1, pt2 = [int(i) for i in a.split("-")]
        if o == "-f":
            pfam_f = a
        if o == "-t":
           ncbi2anc_f = a
        if o == "-l":
            tax_levels = a.split(",")
        if o == "-i":
            last_feat_id = int(a)
        if o == "-o":
            out_table = a
    mi(pt1, pt2, feat_dir, pfam_f, ncbi2anc_f, out_table, tax_levels, int(last_feat_id))
        

