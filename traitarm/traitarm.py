import sys
import pandas as pd
import os
from .reconstruction import merge_annot_and_pt

def loop(parts, cmd_strings, args_dict, fns, out_dir, jobs, cpus):
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
        sys.stdout.write("cat %s/%s | parallel -j %s\n" % (out_dir, fn, cpus))

def reconstruction_cmds(out_dir, parts, tree_unpruned, annotation_tables, phenotype_table, feature_mapping, phenotype_mapping, sample_mapping, do_nested_cv, do_phypat_only, cpus, anno_source, do_normalization, block_cross_validation):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    #merge annotation/features and phenotype data
    pts = merge_annot_and_pt.merge_data(annotation_tables, phenotype_table, out_dir).index.tolist()
    args_dict = {"tree_unpruned" : tree_unpruned, "out_dir" : out_dir, "features" : "feats.txt", "phenotypes" : "phenotypes.txt" , "parts": parts, "phenotype_table" : phenotype_table, "anno_source": anno_source} 
    pts = pd.read_csv("%s/phenotypes.txt" % out_dir, index_col = 0).index.tolist()
    if not do_phypat_only:
        prune_str = "prune_ncbi_tree %(tree_unpruned)s --ncbi_ids %(out_dir)s/ids.txt %(out_dir)s/tree_named.tre --do_name_internal --resolve_polytomy"
        prune_str2 = "prune_ncbi_tree %(tree_unpruned)s --ncbi_ids %(out_dir)s/ids.txt %(out_dir)s/tree.tre --format 5 --failed_to_map %(out_dir)s/mapped_ids.txt --resolve_polytomy" 
        split_str = "summary2gainLoss_input  %(out_dir)s/annot_pheno.dat %(out_dir)s/in --parts %(parts)s --phenotypes %(out_dir)s/%(phenotypes)s --ids %(out_dir)s/mapped_ids.txt --randomize"
        args_dict["tree_unpruned"] = tree_unpruned 
        for i in [prune_str, prune_str2, split_str]:
            sys.stdout.write((i % args_dict) + "\n")
        fns = ["mkdir.sh", "write_config.sh", "gainLoss.sh", "beml_gain.sh", "beml_loss.sh"]
        cmd_strings = ["mkdir %(out_dir)s/run_%(part)s",
                    "write_gainLoss_config   %(out_dir)s/in_%(part)s.fasta %(out_dir)s/tree.tre %(out_dir)s/out_%(part)s %(out_dir)s/gainLoss_param%(part)s.txt",
                    "gainLoss.VR01.266.dRep %(out_dir)s/gainLoss_param%(part)s.txt",
                    "build_edge_matrix_likelihood %(out_dir)s/tree_named.tre %(out_dir)s/annot_pheno.dat  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(out_dir)s/%(phenotypes)s  %(out_dir)s/pgl_matrix_gain_%(part)s",
                    "build_edge_matrix_likelihood %(out_dir)s/tree_named.tre %(out_dir)s/annot_pheno.dat  %(out_dir)s/out_%(part)s/gainLossProbExpPerPosPerBranch.txt %(out_dir)s/in_%(part)s_feats.txt %(out_dir)s/%(phenotypes)s %(out_dir)s/pgl_matrix_loss_%(part)s --consider_loss_events"]
        loop(range(parts), cmd_strings, args_dict, fns, out_dir, range(parts), cpus)
        merge_str1 = "merge_gain_loss %(out_dir)s/pgl_matrix_gain %(out_dir)s/pgl_matrix_gain %(out_dir)s/%(features)s %(parts)s <(echo $'\\taccession\\n%(part)s\\tbla')" 
        merge_str2 = "merge_gain_loss %(out_dir)s/pgl_matrix_loss %(out_dir)s/pgl_matrix_loss %(out_dir)s/%(features)s %(parts)s <(echo $'\\taccession\\n%(part)s\\tbla')" 
        discr_str = "discretize_likelihood_recon %(out_dir)s/pgl_matrix_gain/  %(out_dir)s/pgl_matrix_gain_loss <(echo $'\\taccession\\n%(part)s\\tbla') --in_dir_2 %(out_dir)s/pgl_matrix_loss --discretize_pt_only"
        fns = ["merge_gain.sh", "merge_loss.sh", "discretize.sh"]
        for fn, cmd in zip(fns, [merge_str1, merge_str2, discr_str]):
            loop(pts, [cmd], args_dict, [fn], out_dir, range(parts), cpus)    
    fns = ["learn.sh"]
    args_dict["config"] = os.path.join(os.path.dirname(__file__), "config.json")
    var_names= ["feature_mapping", "phenotype_mapping", "sample_mapping"]
    mappings = [feature_mapping, phenotype_mapping, sample_mapping]
    df_names = ["annot2desc.txt", "pt2desc.txt", "ids2name.txt"]
    for m, n, df in zip(mappings, var_names, df_names):
        if m is None:
            args_dict[n] = "%s/%s" % (out_dir, df)    
        else:
            args_dict[n] = m    
    if do_normalization:
        args_dict['do_normalization'] = '--do_normalization'
    else:
        args_dict['do_normalization'] = ''
    if block_cross_validation:
        args_dict['block_cross_validation'] = '--block_cross_validation %s' % block_cross_validation
        args_dict['block_cross_validation_switch'] = "_block_cv" 
    else:
        args_dict['block_cross_validation'] = ''
        args_dict['block_cross_validation_switch'] = ''
    learn_phypat = "learn %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_phypat_%(block_cross_validation_switch)sout/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(config)s %(sample_mapping)s --with_seed --resume %(do_normalization)s %(block_cross_validation)s"
    cmd_strings = [learn_phypat] 
    #prediction mode
    modes = ["phypat"]
    if not do_phypat_only: 
        modes = modes + ["phypat_pgl", "pgl"]
        fns = fns + ["learn_phypat+pgl.sh", "learn_pgl.sh"]
        learn_phypat_pgl = "learn %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_phypat_pgl%(block_cross_validation_switch)s_out/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(config)s %(sample_mapping)s --is_phypat_and_rec --rec_dir %(out_dir)s/pgl_matrix_gain_loss/ --likelihood_params threshold:0.5,mode:gain_loss --consider_in_recon  %(out_dir)s/mapped_ids.txt --tree_named %(out_dir)s/tree_named.tre --tree %(out_dir)s/tree.tre --with_seed --resume %(block_cross_validation)s"
        learn_pgl =  "learn %(out_dir)s/annot_pheno.dat 10 %(out_dir)s/traitar-model_pgl_%(block_cross_validation_switch)sout/ %(feature_mapping)s <(echo $'\\taccession\\n%(part)s\\tbla') %(config)s %(sample_mapping)s  --rec_dir %(out_dir)s/pgl_matrix_gain_loss/ --likelihood_params threshold:0.5,mode:gain_loss --consider_in_recon  %(out_dir)s/mapped_ids.txt --tree_named %(out_dir)s/tree_named.tre --tree %(out_dir)s/tree.tre --with_seed --resume %(block_cross_validation)s" 
        cmd_strings.append(learn_phypat_pgl)
        cmd_strings.append(learn_pgl)
    if do_nested_cv:
        #add nested cv parameter
        cmd_strings = [i + " --cv_inner 10" for i in cmd_strings]
    loop(pts, cmd_strings, args_dict, fns, out_dir, range(parts), cpus)
    #TODO include traitar new
    model_names = ["%s" % i for i in modes] 
    traitar_new = "traitar new %(out_dir)s/traitar-model_%(mode)s%(block_cross_validation_switch)s_out %(feature_mapping)s %(phenotype_mapping)s %(anno_source)s %(archive_name)s"
    for mode, name in zip(modes, model_names):
        args_dict["mode"] = mode
        args_dict["archive_name"] = name
        args_dict["anno_source"] = anno_source 
        sys.stdout.write((traitar_new % args_dict) + "\n")
