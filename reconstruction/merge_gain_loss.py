import os
import pandas as pd
def reduce(in_dir, out_dir, parts, feats, phenotypes): 
    """concatenate matrices in in_dir"""
    pts = pd.read_csv(phenotypes, sep = "\t", index_col = 0).index.astype('string').tolist()
    feats = pd.read_csv(feats, sep = "\t", index_col = 0).index.tolist()
    for pt in pts:
        m_compl = None 
        for i in range(parts):
            fn = "%s_%s/pt%s.dat" %(in_dir, i, pt)       
            print fn
            if not os.path.exists(fn):
                break
            m = pd.read_csv(fn, sep = "\t", index_col = 0)
            if m_compl is None:
                m_compl = pd.DataFrame(pd.np.zeros((m.shape[0], len(feats) + 1)))
                m_compl.columns = feats + [pt]
                m_compl.index = m.index
            m_compl.loc[:, m.columns] = m
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if not m_compl is None:
            m_compl.to_csv("%s/pt%s.dat" %(out_dir, pt), sep = "\t")

if __name__ == '__main__':
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("out_dir", help='output directory')
    parser.add_argument("in_dir", help='input directory')
    parser.add_argument("feats", help='feature file')
    parser.add_argument("parts", type = int, help='number of splits')
    parser.add_argument("phenotypes", help='phenotype file')
    a = parser.parse_args()
    reduce(**vars(a))
