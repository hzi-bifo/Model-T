import os
import build_edge_matrix_likelihood as beml
import build_edge_matrix as bem
import pandas as ps
import discretize_likelihood_recon as dlr
import shutil
import subprocess

tree_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon.tre"
tree_l_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon.tre.internal_nodes_labeled.newick"
gainLoss = "~/software/GLOOME/gainLoss.VR01.266.dRep"
#TODO for parallelization purposes add a random tag or so to the directory
tmp_dir = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/proper_cv/output/tmp_rec"
phym_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_bin.tsv"
phyn_f = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon_sp+st.taxid_only.txt"
event_f = "gainLossProbExpPerPosPerBranch.txt"


def reconstruct_pt_likelihood(yp_train, model_out, likelihood_params):
    """run gainLoss reconstruction for yp train"""
    tmp_dir = os.path.join(model_out, "tmp_dir")
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)
    #write gainLoss param file to disk:
    f = open("%s/gainLoss_params.txt"%tmp_dir, 'w')
    f.write("_seqFile\t%s/pt_train.fasta\n"%tmp_dir)
    f.write("_treeFile\t%s\n"%tree_f)
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
    all_taxa = ps.read_csv(phyn_f, header = None) 
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
    #TODO consider setting gainLoss option for automatic tree prunning
    command = ["%s %s/gainLoss_params.txt"%(gainLoss, tmp_dir )]
    with open(os.devnull, "w") as fnull:
            result = subprocess.call(command,  stdout = fnull, stderr = fnull, shell = True)
    #map pt events back to the treee
    if likelihood_params["mode"] == "gain_loss" or likelihood_params["mode"] == "gain":
        b = beml.build_edge_matrix_likelihood(tree_l_f, "newick", "%s/pt_train.names.tsv"%tmp_dir, "%s/pt_train.tsv"%tmp_dir, "%s/RESULTS/%s"%(tmp_dir,event_f), use_likelihood = True, use_gain_events = True)
        #create gain output dir
        outdir_g = "%s/gain/"%tmp_dir
        os.mkdir(outdir_g)
        b.get_all_edge_m(0,0,0,0, outdir_g)
        outdir_dlr = "%s/discretized_gain"%tmp_dir
        os.mkdir(outdir_dlr)
        #gain events only case
        if likelihood_params["mode"] == "gain":
            if not os.path.isfile(os.path.join(outdir_g, "pt0.dat")):
                raise ValueError("Phenotype has no events associated with it")    
            m = dlr.threshold_matrix(outdir_g, float(likelihood_params["threshold"]), outdir_dlr, is_internal = True) 
    if likelihood_params["mode"] == "gain_loss" or likelihood_params["mode"] == "loss":
        b = beml.build_edge_matrix_likelihood(tree_l_f, "newick", "%s/pt_train.names.tsv"%tmp_dir, "%s/pt_train.tsv"%tmp_dir, "%s/RESULTS/%s"%(tmp_dir,event_f), use_likelihood = True, use_gain_events = False)
        #create gain output dir
        outdir_l = "%s/loss/"%tmp_dir
        os.mkdir(outdir_l)
        b.get_all_edge_m(0,0,0,0, outdir_l)
        outdir_dlr = "%s/discretized_loss"%tmp_dir
        os.mkdir(outdir_dlr)
        #loss events only case
        if likelihood_params["mode"] == "loss":
            if not os.path.isfile(os.path.join(outdir_l, "pt0.dat")):
                raise ValueError("Phenotype has no events associated with it")    
            m = dlr.threshold_matrix(outdir_l, float(likelihood_params["threshold"]), outdir_dlr, is_internal = True) 
    if likelihood_params["mode"] == "gain_loss": 
        #gain and loss events combined
        if not os.path.isfile(os.path.join(outdir_g, "pt0.dat")) or not os.path.isfile(os.path.join(outdir_l, "pt0.dat")):
            raise ValueError("Phenotype has no events associated with it")  
        outdir_dlr = "%s/discretized_gain_loss"%tmp_dir
        os.mkdir(outdir_dlr)
        m = dlr.threshold_matrix(outdir_g, float(likelihood_params["threshold"]), outdir_dlr, outdir_l, is_internal = True) 
    return m



def reconstruct_pt_parsimony(yp_train, parsimony_params):
    raise Exception("parsimony not yet implemented")
