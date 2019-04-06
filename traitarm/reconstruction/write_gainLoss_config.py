def write_config(seq_f, tree_f, out_dir, result_f):
    with open(result_f, 'w') as f: 
        f.write("_seqFile\t%s\n" % seq_f)
        f.write("_treeFile\t%s\n" % tree_f)
        f.write("_calculateAncestralReconstruct\t1\n")
        f.write("_outDir\t%s\n" % out_dir)

