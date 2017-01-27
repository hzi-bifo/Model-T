import sys
def reconstruction_cmds(out_dir, parts, tree_unpruned, phenotypes, features, annotation_table, ids, code_dir):
    prune_str = "python %(code_dir)s/reconstruction/prune_ncbi_tree.py %(tree_unpruned)s --ncbi_ids %(ids)s %(out_dir)s/tree_named.tre --do_name_internal"
    prune_str2 = "python %(code_dir)s/reconstruction/prune_ncbi_tree.py %(tree_unpruned)s --ncbi_ids %(ids)s %(out_dir)s/tree.tre --format 5 --failed_to_map %(out_dir)s/failed_to_map_ids.txt" 
    split_str = "python %(code_dir)s/reconstruction/summary2gainLoss_input.py  %(annotation_table)s %(out_dir)s/in --parts %(parts)s --phenotypes %(phenotypes)s"
    args_dict = {"tree_unpruned" : tree_unpruned, "ids" : ids, "out_dir" : out_dir, "features" : features, "phenotypes" : phenotypes, "annotation_table" : annotation_table, "parts": parts, "code_dir" : code_dir} 
    for i in [prune_str, prune_str2, split_str]:
        sys.stdout.write((i % args_dict) + "\n")
    for i in range(parts): 
        cmd_strings = ["mkdir %(out_dir)s/run_%(part)s",
                "python %(code_dir)s/reconstruction/write_gainLoss_config.py   %(out_dir)s/in_%(part)s.fasta %(out_dir)s/tree.tre %(out_dir)s/out_%(part)s %(out_dir)s/gainLoss_param%(part)s.txt",
                "%(code_dir)s/gainLoss.VR01.266.dRep %(out_dir)s/gainLoss_param%(part)s.txt",
                "python %(code_dir)s/reconstruction/build_edge_matrix_likelihood.py %(out_dir)s/tree_named.tre %(annotation_table)s  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(phenotypes)s  %(out_dir)s/pgl_matrix_gain_%(part)s",
                "python %(code_dir)s//reconstruction/build_edge_matrix_likelihood.py %(out_dir)s/tree_named.tre %(annotation_table)s  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(phenotypes)s %(out_dir)s/pgl_matrix_loss_%(part)s --consider_loss_events"]
        args_dict["part"] = i
        cmd_strings_args = map(lambda x: x % args_dict, cmd_strings) 
        with open("%s/cmds%s.sh" % (out_dir, i), 'w') as f:
            f.write("\n".join(cmd_strings_args))
        sys.stdout.write("bash %s/cmds%s.sh\n" % (out_dir, i))
    merge_str1 = "python %(code_dir)s/reconstruction/merge_gain_loss.py %(out_dir)s/pgl_matrix_gain %(out_dir)s/pgl_matrix_gain %(features)s 2 %(phenotypes)s" 
    merge_str2 = "python %(code_dir)s/reconstruction/merge_gain_loss.py %(out_dir)s/pgl_matrix_loss %(out_dir)s/pgl_matrix_loss %(features)s 2 %(phenotypes)s" 
    discr_str = "python %(code_dir)s/learning/discretize_likelihood_recon.py -d %(out_dir)s/pgl_matrix_gain/ -t 0.5 -o %(out_dir)s/pgl_matrix_gain_loss -f %(out_dir)s/pgl_matrix_loss -g"
    for i in [merge_str1, merge_str2, discr_str]:
        sys.stdout.write((i % args_dict) + "\n")

        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("generate summary matrix from the filtered best hmmer annotation files")
    parser.add_argument("out_dir", help='summary matrix output file')
    parser.add_argument("parts", type = int, help='summary matrix output file')
    parser.add_argument("tree_unpruned", help='summary matrix output file')
    parser.add_argument("phenotypes", help='summary matrix output file')
    parser.add_argument("features", help='summary matrix output file')
    parser.add_argument("annotation_table", help='summary matrix output file')
    parser.add_argument("ids", help='ids to prune input tree')
    parser.add_argument("code_dir", help='traitar-model code directory')
    args = parser.parse_args()
    reconstruction_cmds(**vars(args))
