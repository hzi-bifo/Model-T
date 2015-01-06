import os
import build_edge_matrix_likelihood as beml
import build_edge_matrix as bem
import join_likelihood_recon as jlr
import pandas as ps
import discretize_likelihood_recon as dlr
import shutil
import subprocess

def reconstruct_pt_likelihood(yp_train, model_out, config,   likelihood_params, pt_out, ofold, ifold = None):
    """run gainLoss reconstruction for yp train"""
    #output directory for the gainLoss output (tmp does not mean this directory is removed afterwards) 
    if 'gainLoss_ref' in config:
        model_out = config['gainLoss_ref']
    tmp_dir = os.path.join(model_out, "gainLoss_ofold%s" %ofold)
    if not ifold is None:
        tmp_dir = "%s_ifold%s"%(tmp_dir, ifold) 
    #only compute the ml reconstruction if the directory doesn't already exist or in the config a reference cv reconstruction is specified
    if not os.path.exists(tmp_dir):
        if 'gainLoss_ref' in config:
            sys.exit("fatal error; gainLoss reconstruction %s should already exist!" % tmp_dir)
        os.mkdir(tmp_dir)
        #write gainLoss param file to disk:
        f = open("%s/gainLoss_params.txt"%tmp_dir, 'w')
        f.write("_seqFile\t%s/pt_train.fasta\n"%tmp_dir)
        f.write("_treeFile\t%s\n"%config['tree_f'])
        f.write("_accountForMissingData\t0\n")
        f.write("_outDir\t%s/RESULTS\n"%tmp_dir)
        f.close()
        
        #write pt to disk
        f = open("%s/pt_train.fasta"%tmp_dir, 'w')
        g = open("%s/pt_train.tsv"%tmp_dir, 'w')
        h = open("%s/pt_train.names.tsv"%tmp_dir, 'w')
        for l in yp_train.index:
            f.write(">%s\n"%l)
            f.write("%i\n"%yp_train.loc[l])
            g.write("%i\n"%yp_train.loc[l])
            h.write("%s\n"%l)
        #gainLoss programm also needs the taxa with missing pt 
        all_taxa = ps.read_csv(config['phyn_f'], header = None) 
        for l in all_taxa.iloc[:, 0]:
            if not l in yp_train.index:
                f.write(">%s\n"%l)
                f.write("?\n")
                g.write("?\n")
                h.write("%s\n"%l)
        f.close()
        g.close()
        h.close()
        #execute gainLoss for that phenotype 
        command = ["%s %s/gainLoss_params.txt"%(config['gainLoss'], tmp_dir )]
        with open(os.devnull, "w") as fnull:
                result = subprocess.call(command,  stdout = fnull, stderr = fnull, shell = True)
    #map pt events back to the tree
    if likelihood_params["mode"] == "gain_loss" or likelihood_params["mode"] == "gain":
        b = beml.build_edge_matrix_likelihood(config['tree_l_f'], "newick", "%s/pt_train.names.tsv"%tmp_dir, "%s/pt_train.tsv"%tmp_dir, "%s/RESULTS/%s"%(tmp_dir,"gainLossProbExpPerPosPerBranch.txt"), use_likelihood = True, use_gain_events = True)
        #create gain output dir
        outdir_g = "%s/gain/"%tmp_dir
        if not os.path.exists(outdir_g):
            os.mkdir(outdir_g)
        if "continuous_target" in likelihood_params:
            m = b.get_all_edge_m(0,0,0,0, outdir_g, is_internal = True)
        else: 
            b.get_all_edge_m(0,0,0,0, outdir_g)
            outdir_dlr = "%s/discretized_gain"%tmp_dir
            if not os.path.exists(outdir_g):
                os.mkdir(outdir_dlr)
            #gain events only case
            if likelihood_params["mode"] == "gain":
                if not os.path.isfile(os.path.join(outdir_g, "pt0.dat")):
                    raise ValueError("Phenotype has no events associated with it")    
                m = dlr.threshold_matrix(outdir_g, float(likelihood_params["threshold"]), outdir_dlr, is_internal = True) 
    if likelihood_params["mode"] == "gain_loss" or likelihood_params["mode"] == "loss":
        b = beml.build_edge_matrix_likelihood(config['tree_l_f'], "newick", "%s/pt_train.names.tsv"%tmp_dir, "%s/pt_train.tsv"%tmp_dir, "%s/RESULTS/%s"%(tmp_dir,"gainLossProbExpPerPosPerBranch.txt"), use_likelihood = True, use_gain_events = False)
        #create gain output dir
        outdir_l = "%s/loss/"%tmp_dir
        if not os.path.exists(outdir_l):
            os.mkdir(outdir_l)
        if "continuous_target" in likelihood_params:
            m = b.get_all_edge_m(0,0,0,0, outdir_l, is_internal = True)
        else:
            b.get_all_edge_m(0,0,0,0, outdir_l)
            outdir_dlr = "%s/discretized_loss"%tmp_dir
            if not os.path.exists(outdir_dlr):
                os.mkdir(outdir_dlr)
            #loss events only case
            if likelihood_params["mode"] == "loss":
                if not os.path.isfile(os.path.join(outdir_l, "pt0.dat")):
                    raise ValueError("Phenotype has no events associated with it")    
                m = dlr.threshold_matrix(outdir_l, float(likelihood_params["threshold"]), outdir_dlr, is_internal = True) 
    if likelihood_params["mode"] == "gain_loss": 
        if not os.path.isfile(os.path.join(outdir_g, "pt0.dat")) or not os.path.isfile(os.path.join(outdir_l, "pt0.dat")):
            raise ValueError("Phenotype has no events associated with it")  
        #gain and loss events combined
        if "continuous_target" in likelihood_params:
            #join gain and loss event matrices
            m = jlr.threshold_matrix(outdir_g, "",  outdir_l, is_internal = True)
            
        else:
            if not os.path.isfile(os.path.join(outdir_g, "pt0.dat")) or not os.path.isfile(os.path.join(outdir_l, "pt0.dat")):
                raise ValueError("Phenotype has no events associated with it")  
            outdir_dlr = "%s/discretized_gain_loss"%tmp_dir
            if not os.path.exists(outdir_dlr):
                os.mkdir(outdir_dlr)
            m = dlr.threshold_matrix(outdir_g, float(likelihood_params["threshold"]), outdir_dlr, outdir_l, is_internal = True) 
    #drop first column that is due to genotype reconstruction and either empty or duplicate of phenotype column
    m = ps.DataFrame(m.iloc[:, 1] )
    return m



def reconstruct_pt_parsimony(yp_train, parsimony_params):
    raise Exception("parsimony not yet implemented")
