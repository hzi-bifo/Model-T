import pandas as pd
import sys
import argparse



def split(m, out_f, parts, randomize = True):
    """split matrix into <parts> parts"""
    if randomize:
        m = m.sample(m.shape[0], index = 1)
    #number of elements in each subset
    part_size = m.shape[1] / parts
    print part_size
    #number of subsets which have length part_size + 1 
    parts_extended = m.shape[1] % parts
    #process subsets with size part_size
    for i in range(parts - parts_extended):
        part = m.iloc[:,i * part_size : (i + 1) * part_size]
        to_fasta(part, out_f,  i)
    #process subsets with size part_size + 1
    for i in range(parts_extended):
        part = m.iloc[:, (parts - parts_extended) * part_size + i * (part_size + 1) : (parts - parts_extended) * part_size + (i + 1) * (part_size + 1)]
        to_fasta(part, out_f, parts - parts_extended + i)

def to_fasta(m, out_n, index = None):
    """transform matrix into FASTA format"""
    #binarize
    m.fillna(0.5, inplace = True)
    m_isnan = (m == 0.5)
    m = (m > 0).astype('int')
    m[m_isnan] = "?"
    with open(out_n if index is None else "%s_%s.fasta" %(out_n, index) , 'w') as out_f:
        for i in m.index:
            out_f.write(">%s\n" % i)
            out_f.write("%s\n" % "".join(m.loc[i,].astype('string').tolist()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser("reconstruct likelihood matrix from gainLoss output")
    parser.add_argument("in_matrix", help='phylogenetic tree')
    parser.add_argument("out_f", help='phylogenetic tree')
    parser.add_argument("--parts", type = int, default = None, help='split matrix into <parts> subsets and output each sub matrix as fasta')
    parser.add_argument("--randomize", action = 'store_true', help='set if output should be randomized')
    a = parser.parse_args()
    m = pd.read_csv(a.in_matrix, sep = "\t", index_col = 0, na_values = "?")
    split(m, a.out_f, a.parts, a.randomize)


