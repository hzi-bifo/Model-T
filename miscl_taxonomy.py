import NCBI_tree_of_life
import pandas as ps
import ete2

def map_miscl2taxonomy(miscl_m_f, out_f):
    #ncbi_tree = NCBI_tree_of_life.ncbiTaxonomy(nodes_f = "/net/refdata/static/ncbi-taxonomy_20130503/nodes.dmp", names_f = "/net/refdata/static/ncbi-taxonomy_20130503/names.dmp")
    #ncbi_tree = NCBI_tree_of_life.ncbiTaxonomy(nodes_f = "nodes.dmp", names_f = "names.dmp")
    #read in the matrix of misclassified species statistics
    miscl_m = ps.read_csv(miscl_m_f, sep = "\t", index_col = 0) 
    miscl_m.index = [str(i) for i in miscl_m.index]
    print miscl_m
    #prune ncbi tree to the GIDEON species
    #ncbi_tree.tree.prune(miscl_m.index)
    #ncbi_tree.tree.write(features = ["name", "scientific_name", "rank"], outfile = "gideon_tree.nwck")
    ncbi_tree = ete2.Tree("/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gideon_tree.nwck")
    miscl_m_extd = ps.DataFrame(ps.np.zeros(shape = (len([i for i in ncbi_tree.traverse()]), 4)))
    miscl_m_extd.index = [i.name for i in ncbi_tree.traverse()]
    miscl_m_extd.columns = miscl_m.columns[0:4]
    miscl_m_extd.loc[miscl_m.index, miscl_m.columns[0:4]] = miscl_m.iloc[:, 0:4]
    #propagate the number of samples / misclassified samples up to the root
    for n in ncbi_tree.iter_leaves(): 
        n_leave = n
        while n.name != "NoName":
            print n.name
            n = n.up
            miscl_m_extd.loc[n.name, ] += miscl_m_extd.loc[n_leave.name, ]
    bacc = miscl_m_extd.apply(lambda x: (1 - (x.iloc[2]/float(x.iloc[0]) + x.iloc[3]/float(x.iloc[1]))/2), axis = 1)
    print ps.concat([miscl_m_extd, bacc], axis = 1).to_csv(out_f, sep = "\t")

if __name__ == "__main__":
    map_miscl2taxonomy("misclassified.tsv", "misclassified_extd.tsv")
    

