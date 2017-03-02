import sys
import pandas as pd
import os
def loop(parts, cmd_strings, args_dict, fns, out_dir, jobs):
    for fn in fns:
        with open("%s/%s" % (out_dir, fn), 'w') as f:
            pass 
    for i in parts: 
        args_dict["part"] = i
        cmd_strings_args = map(lambda x: x % args_dict, cmd_strings) 
        for fn, cmd in zip(fns, cmd_strings_args):
            with open("%s/%s" % (out_dir, fn), 'a') as f:
                f.write("%s\n" % cmd)
    for fn in fns:
        sys.stdout.write("cat %s/%s | parallel -j %s\n" % (out_dir, fn, len(jobs)))

def reconstruction_cmds(out_dir, parts, tree_unpruned, annotation_table, code_dir, phenotype_table, feature_mapping,  phenotype_mapping, sample_mapping):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    merge_str = "python /home/aweimann/code/phenotypes/traitar-model/reconstruction/merge_annot_and_pt.py %(annotation_table)s %(phenotype_table)s %(out_dir)s"
    prune_str = "python %(code_dir)s/reconstruction/prune_ncbi_tree.py %(tree_unpruned)s --ncbi_ids %(out_dir)s/ids.txt %(out_dir)s/tree_named.tre --do_name_internal --resolve_polytomy"
    prune_str2 = "python %(code_dir)s/reconstruction/prune_ncbi_tree.py %(tree_unpruned)s --ncbi_ids %(out_dir)s/ids.txt %(out_dir)s/tree.tre --format 5 --failed_to_map %(out_dir)s/mapped_ids.txt --resolve_polytomy" 
    split_str = "python %(code_dir)s/reconstruction/summary2gainLoss_input.py  %(out_dir)s/annot_pheno.dat %(out_dir)s/in --parts %(parts)s --phenotypes %(out_dir)s/%(phenotypes)s --ids %(out_dir)s/mapped_ids.txt --randomize"
    args_dict = {"tree_unpruned" : tree_unpruned, "out_dir" : out_dir, "features" : "feats.txt", "phenotypes" : "phenotypes.txt" , "annotation_table" : annotation_table, "parts": parts, "code_dir" : code_dir, "phenotype_table" : phenotype_table} 
    for i in [merge_str, prune_str, prune_str2, split_str]:
        sys.stdout.write((i % args_dict) + "\n")
    fns = ["mkdir.sh", "write_config.sh", "gainLoss.sh", "beml_gain.sh", "beml_loss.sh"]
    cmd_strings = ["mkdir %(out_dir)s/run_%(part)s",
                "python %(code_dir)s/reconstruction/write_gainLoss_config.py   %(out_dir)s/in_%(part)s.fasta %(out_dir)s/tree.tre %(out_dir)s/out_%(part)s %(out_dir)s/gainLoss_param%(part)s.txt",
                "%(code_dir)s/gainLoss.VR01.266.dRep %(out_dir)s/gainLoss_param%(part)s.txt",
                "python %(code_dir)s/reconstruction/build_edge_matrix_likelihood.py %(out_dir)s/tree_named.tre %(out_dir)s/annot_pheno.dat  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(out_dir)s/%(phenotypes)s  %(out_dir)s/pgl_matrix_gain_%(part)s",
                "python %(code_dir)s//reconstruction/build_edge_matrix_likelihood.py %(out_dir)s/tree_named.tre %(out_dir)s/annot_pheno.dat  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(out_dir)s/%(phenotypes)s %(out_dir)s/pgl_matrix_loss_%(part)s --consider_loss_events"]
    loop(range(parts), cmd_strings, args_dict, fns, out_dir, range(parts))
    pts = pd.read_csv("%s/phenotypes.txt" % out_dir, index_col = 0).index.tolist()
    merge_str1 = "python %(code_dir)s/reconstruction/merge_gain_loss.py %(out_dir)s/pgl_matrix_gain %(out_dir)s/pgl_matrix_gain %(out_dir)s/%(features)s %(parts)s <(echo $'\\taccession\\n%(part)s\\tbla')" 
    merge_str2 = "python %(code_dir)s/reconstruction/merge_gain_loss.py %(out_dir)s/pgl_matrix_loss %(out_dir)s/pgl_matrix_loss %(out_dir)s/%(features)s %(parts)s <(echo $'\\taccession\\n%(part)s\\tbla')" 
    discr_str = "python %(code_dir)s/learning/discretize_likelihood_recon.py  %(out_dir)s/pgl_matrix_gain/  %(out_dir)s/pgl_matrix_gain_loss <(echo $'\\taccession\\n%(part)s\\tbla') --in_dir_2 %(out_dir)s/pgl_matrix_loss --discretize_pt_only"
    fns = ["merge_gain.sh", "merge_loss.sh", "discretize.sh"]
    for fn, cmd in zip(fns, [merge_str1, merge_str2, discr_str]):
        loop(pts, [cmd], args_dict, [fn], out_dir, range(parts))    
    fns = ["learn.sh", "learn_phypat+pgl.sh", "learn_pgl.sh"]
    var_names= ["feature_mapping", "phenotype_mapping", "sample_mapping"]
    mappings = [feature_mapping, phenotype_mapping, sample_mapping]
    df_names = ["annot2desc.txt", "pt2desc.txt", "ids2name.txt"]
    for m, n, df in zip(mappings, var_names, df_names):
        if m is None:
            args_dict[n] = "%s/%s" % (out_dir, df)    
        else:
            args_dict[n] = m    
    cmd_strings = ["python ~/code/phenotypes/traitar-model/learning/learn.py %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_phypat_out/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(code_dir)s/learning/config.json %(sample_mapping)s --with_seed --resume ",
            "python ~/code/phenotypes/traitar-model/learning/learn.py %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_phypat_pgl_out/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(code_dir)s/learning/config.json %(sample_mapping)s --is_phypat_and_rec --rec_dir %(out_dir)s/pgl_matrix_gain_loss/ --likelihood_params threshold:0.5,mode:gain_loss --consider_in_recon  %(out_dir)s/mapped_ids.txt --tree_named %(out_dir)s/tree_named.tre --tree %(out_dir)s/tree.tre --with_seed --resume ", 
            "python ~/code/phenotypes/traitar-model/learning/learn.py %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_pgl_out/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(code_dir)s/learning/config.json %(sample_mapping)s  --rec_dir %(out_dir)s/pgl_matrix_gain_loss/ --likelihood_params threshold:0.5,mode:gain_loss --consider_in_recon  %(out_dir)s/mapped_ids.txt --tree_named %(out_dir)s/tree_named.tre --tree %(out_dir)s/tree.tre --with_seed --resume " ]
    loop(pts, cmd_strings, args_dict, fns, out_dir, range(parts))


        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("generate summary matrix from the filtered best hmmer annotation files")
    parser.add_argument("out_dir", help='summary matrix output file')
    parser.add_argument("parts", type = int, help='summary matrix output file')
    parser.add_argument("tree_unpruned", help='summary matrix output file')
    parser.add_argument("phenotype_table", help='summary matrix output file')
    parser.add_argument("--phenotype_mapping", help='summary matrix output file')
    parser.add_argument("--feature_mapping", help='summary matrix output file')
    parser.add_argument("annotation_table", help='summary matrix output file')
    parser.add_argument("--sample_mapping", help='ids to prune input tree')
    parser.add_argument("code_dir", help='traitar-model code directory')
    args = parser.parse_args()
    reconstruction_cmds(**vars(args))
