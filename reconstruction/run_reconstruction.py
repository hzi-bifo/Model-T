import sys
def reconstruction_cmds(out_dir, parts, tree, tree_named, phenotypes, features, annotation_table):
    split_str = "python ~/code/phenotypes/traitar-model/reconstruction/summary2gainLoss_input.py  %(annotation_table)s %(out_dir)s/in_ --parts %(parts)s"
    args_dict = {"tree_named" : tree_named, "tree" : tree, "out_dir" : out_dir, "features" : features, "phenotypes" : phenotypes, "annotation_table" : annotation_table, "parts": parts} 
    sys.stdout.write((split_str % args_dict) + "\n")
    for i in range(parts): 
        cmd_strings = ["mkdir %(out_dir)s/run_%(part)s",
                "python ~/code/phenotypes/traitar-model/reconstruction/write_gainLoss_config.py   %(out_dir)s/in_%(part)s.fasta %(tree)s %(out_dir)s/out_%(part)s %(out_dir)s/gainLoss_param%(part)s.txt",
                "~/code/phenotypes/traitar-model/gainLoss.VR01.266.dRep %(out_dir)s/gainLoss_param%(part)s.txt",
                "python ~/code/phenotypes/traitar-model/reconstruction/build_edge_matrix_likelihood.py %(tree_named)s %(annotation_table)s  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(features)s %(phenotypes)s  pgl_matrix_gain",
                "python ~/code/phenotypes/traitar-model/reconstruction/build_edge_matrix_likelihood.py %(tree_named)s %(annotation_table)s  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(features)s %(phenotypes)s pgl_matrix_gain"]
        args_dict["part"] = i
        cmd_strings_args = map(lambda x: x % args_dict, cmd_strings) 
        with open("%s/cmds%s.sh" % (out_dir, i), 'w') as f:
            f.write("\n".join(cmd_strings_args))
        sys.stdout.write("bash %s/cmds%s.sh\n" % (out_dir, i))
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("generate summary matrix from the filtered best hmmer annotation files")
    parser.add_argument("out_dir", help='summary matrix output file')
    parser.add_argument("parts", type = int, help='summary matrix output file')
    parser.add_argument("tree", help='summary matrix output file')
    parser.add_argument("tree_named", help='summary matrix output file')
    parser.add_argument("phenotypes", help='summary matrix output file')
    parser.add_argument("features", help='summary matrix output file')
    parser.add_argument("annotation_table", help='summary matrix output file')
    args = parser.parse_args()
    reconstruction_cmds(**vars(args))
