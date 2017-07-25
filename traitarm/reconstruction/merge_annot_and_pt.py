import pandas as pd
def merge_data(annotation, phenotypes, out_dir):
    #read in annotation from potentially different sources
    annots = [pd.read_csv(annot, index_col = 0, sep = "\t") for annot in annotation]
    #find samples common to all annotation tables
    common_ids = set(annots[0].index.tolist())
    for annot in annots:
        common_ids = common_ids.intersection(set(annot.index.tolist()))
    common_ids = list(common_ids)
    annot_m = pd.concat([annot.loc[common_ids, ] for annot in annots], axis = 1)
    pheno_m = pd.read_csv(phenotypes, index_col = 0, sep = "\t")
    is_in_annot_index = [i in pheno_m.index for i in annot_m.index]
    is_in_pheno_index = [i in annot_m.index for i in pheno_m.index]
    annot_m = annot_m.loc[is_in_annot_index,]
    pheno_m = pheno_m.loc[is_in_pheno_index,]
    #make sure pheno_m uses ? as missing characters
    pheno_m.where(pd.isnull(pheno_m), "?")
    #create dummy annotation to description mapping 
    anno_map = pd.DataFrame(annot_m.columns, index = annot_m.columns) 
    anno_map.columns = ["description"]
    anno_map.to_csv("{}/annot2desc.txt".format(out_dir), sep = '\t')
    #create dummy phenotype to description mapping 
    pheno_map = pd.DataFrame(pheno_m.columns, index = pheno_m.columns) 
    pheno_map.columns = ["accession"] 
    pheno_map.to_csv("{}/pt2desc.txt".format(out_dir), sep = '\t')
    #write phenotypes and feature list to disk
    annot_m.columns.name = "feats"
    annot_m.T.to_csv("{}/feats.txt".format(out_dir), columns = [], sep = '\t')
    pheno_m.columns.name = "phenotypes"
    pheno_m.T.to_csv("{}/phenotypes.txt".format(out_dir), columns = [], sep = '\t')
    
    #create a joint annotation phenotype table
    m = pd.concat([annot_m, pheno_m], axis = 1)
    m.to_csv("{}/annot_pheno.dat".format(out_dir), sep = '\t')
    m.index.name = "ids"
    m.to_csv("{}/ids.txt".format(out_dir), columns = [], sep = '\t')
    #create dummy ids mapping
    ids_map = pd.DataFrame(m.index, index = m.index) 
    ids_map.to_csv("{}/ids2name.txt".format(out_dir),  sep = '\t')
    return pheno_m

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("learn antibiotic resistance models")
    parser.add_argument("annotation", help= 'annotation table', nargs='+')
    parser.add_argument("phenotypes", help= 'phenotype table')
    parser.add_argument("out_dir", help= 'target')
    args = parser.parse_args()
    merge_data(**vars(args))
