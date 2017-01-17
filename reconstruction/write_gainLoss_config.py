import argparse

def write_config(seq_f, tree_f, out_dir, result_f):
    with open(result_f, 'w') as f: 
        f.write("_seqFile\t%s\n" % seq_f)
        f.write("_treeFile\t%s\n" % tree_f)
        f.write("_calculateAncestralReconstruct\t1\n")
        f.write("_outDir\t%s\n" % out_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("write gainLoss param file")
    parser.add_argument("seq_f", help='binary character file in fasta format')
    parser.add_argument("tree_f", help='phylogenetic tree in newick format')
    parser.add_argument("out_dir", help='output directory for gainLoss run')
    parser.add_argument("result_f", help='gainLoss param file')
    a = parser.parse_args()
    write_config(a.seq_f, a.tree_f, a.out_dir, a.result_f)
