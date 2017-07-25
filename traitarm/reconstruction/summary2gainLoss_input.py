import pandas as pd
import sys
import argparse
import random
random.seed(0)

def split(m, out_f, parts, randomize = True, phenotypes = None, ids = None):
    """split matrix into <parts> parts"""
    #read in phenotypes
    #split data table into phenotypes and features
    if not phenotypes is None:
        pts = pd.read_csv(phenotypes, sep = "\t", index_col = 0)
        pts.index = pts.index.astype('string')
        idispheno = [i not in pts.index for i in m.columns]
        m_pheno = m.loc[:, pts.index] 
        m = m.loc[:, idispheno] 
    else: 
        m_pheno = None
    #randomize to avoid any biases
    if randomize:
        columns = m.columns.values.copy()
        random.shuffle(columns)
        m = m.loc[:, columns]
    #number of elements in each subset
    part_size = m.shape[1] / parts
    #number of subsets which have length part_size + 1 
    parts_extended = m.shape[1] % parts
    #process subsets with size part_size
    for i in range(parts - parts_extended):
        part = m.iloc[:,i * part_size : (i + 1) * part_size]
        #write features per subset
        part.T.to_csv("%s_%s_feats.txt" %(out_f, i), columns=[])
        if m_pheno is not None:
            #always add phenotypes
            part = pd.concat([part, m_pheno], axis = 1)
        to_fasta(part, out_f,  i)
    #process subsets with size part_size + 1
    for i in range(parts_extended):
        part = m.iloc[:, (parts - parts_extended) * part_size + i * (part_size + 1) : (parts - parts_extended) * part_size + (i + 1) * (part_size + 1)]
        #write features per subset
        part.T.to_csv("%s_%s_feats.txt" %(out_f, parts - parts_extended + i), columns=[])
        if m_pheno is not None:
            #always add phenotypes
            part = pd.concat([part, m_pheno], axis = 1)
        to_fasta(part, out_f, parts - parts_extended + i)

def to_fasta(m, out_n, index = None):
    """transform matrix into FASTA format"""
    #binarize and replace nas with ?
    m.fillna(0.5, inplace = True)
    m_isnan = (m == 0.5)
    m = (m > 0).astype('int')
    m[m_isnan] = "?"
    #write fast to disk
    with open(out_n if index is None else "%s_%s.fasta" %(out_n, index) , 'w') as out_f:
        for i in m.index:
            out_f.write(">%s\n" % i)
            out_f.write("%s\n" % "".join(m.loc[i,].astype('string').tolist()))
