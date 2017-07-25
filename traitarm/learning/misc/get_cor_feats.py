import scipy.stats as stats
import pandas as ps

#phypat_f  = "../../pfams_pts_counts_named.tsv"
#ml_f = "../gain_loss_prob0.5/pt%s.dat"
#target = "Motile"
#feat_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus/likelihood/ll_gain_loss_threshold0.5_percfeat1.0_shuffled-b/%s_majority_features+weights.txt"
#pfam_mapping_f = "../../pfam_pts_names_nl_desc.txt"
def setup(phypat_f, ml_f, target, feat_f, pfam_mapping_f):
    phypat = ps.read_csv(phypat_f, index_col = 0, sep = "\t", header = None)
    pfams = ps.read_csv(pfam_mapping_f, sep = "\t", header = None)
    pfams.index = pfams.loc[:, 1]
    target_id = pfams.loc[target, 0] - 1
    ml = ps.read_csv(ml_f % target_id, sep = "\t", index_col = 0, header = None)
    phypat_red = phypat.loc[phypat.iloc[:, target_id] != "?",]
    phypat_red_pt = phypat_red.iloc[:, range(8682) + [target_id]]
    phypat_red_pt.columns = ml.columns
    phypat_red_pt.columns = pfams.iloc[range(8682) + [target_id], 1]
    ml.columns = phypat_red_pt.columns
    phypat_red_pt.loc[:, target]  = phypat_red_pt.loc[:, target].astype('int')
    phypat_red_pt = (phypat_red_pt > 0).astype("int")
    phypat_ml = ps.concat([phypat_red_pt, ml], axis = 0)
    #read in motility features 
    target_feats = ps.read_csv(feat_f%target_id, sep = '\t')
    return phypat_ml, phypat_red_pt, target_feats, pfams

def get_top_cor(pfam, df, pfam_mapping):
    a = df.apply(lambda x: stats.pearsonr(x, df.loc[:, pfam])[0], axis = 0)
    a.sort(ascending = False)
    return ps.concat([a, pfam_mapping.loc[a.index,2]], axis = 1)

def get_selected_cor(target, pfams, df, pfam_mapping):
    a = get_top_cor(target, df, pfam_mapping)
    return a.loc[pfams,].sort(columns = 0, ascending = False)

def get_cor_extd_feats(pfams, df, k, pfam_mapping):
    l = []
    for j in pfams:
        top = get_top_cor(j, df, pfam_mapping).iloc[:10,]
        l.append(top)
    #res = ps.DataFrame(l)
    #res.columns = pfams
    return l 

def get_top_cor_extd_feats(target, df, j, k, pfam_mapping):
    top = get_selected_cor(target, target_feats.loc[:, "Pfam_acc"], df, pfam_mapping)
    return get_cor_extd_feats(top.iloc[:j,].index, df, k, pfam_mapping)
 


#phypat_ml, phypat_red_pt, target_feats, pfam_mapping = setup(phypat_f, ml_f, target, feat_f, pfam_mapping_f) 
#cor = get_top_cor_extd_feats("Motile", phypat_red_pt,  10, 10, pfam_mapping)

