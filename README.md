# traitar-model
Learning phenotype classification models from protein family phyletic patterns
## Example usage
Construct a plant biomass degradation phenotype model based on the example data

1. Convert table of dbCAN families (features) and phenotypes (targets) to FASTA input required for the gainLoss program
  
  ``python summary2gainLoss_input.py < dbcan_summary_stol_pt_named.txt > dbcan.fasta``
  
2. Prune sequenced tree of life to the plant biomass degrader and non-degraders; one version with named internal nodes and one unnamed

  `python prune_ncbi_tree.py  ncbi-bacteria.tre --ncbi_ids pbd_taxids_stol.txt pbd_named.tre  --do_name_internal`
  
  `python prune_ncbi_tree.py ncbi-bacteria.tre --ncbi_ids pbd_taxids_stol.txt pbd.tre --format 5  --failed_to_map     pbd_taxids_stol_mapped.txt`

3. Run reconstruction with gainLoss
  
  `gainLoss.VR01.266.dRep gainLoss_params.txt`
  
4. Reconstruct likelihood matrix
  
  Loss matrix
  
  `python build_edge_matrix_likelihood.py pbd_named.tre  dbcan_summary_stol_pt_named.txt   RESULTS/gainLossProbExpPerPosPerBranch.txt  dbcans.txt phenotypes.txt pgl_matrix_gain`
  
  Gain matrix
  
  `python build_edge_matrix_likelihood.py pbd_named.tre  dbcan_summary_stol_pt_named.txt   RESULTS/gainLossProbExpPerPosPerBranch.txt  dbcans.txt phenotypes.txt pgl_matrix_gain`
  
5. Discretize and join likelihood matrix
  
  `python discretize_likelihood_recon.py -d pgl_matrix_gain/ -t 0.5 -o pgl_matrix_gain_loss -f pgl_matrix_loss -g`
  
6. Build classficiation models
  
  Only with phyletic patterns
  
  `python learning.py  dbcan_summary_stol_pt_named.txt   10 ll_pbd dbcan2desc_name_full.txt pt_all_settings2acc.txt  config.json`
  
  With gain and loss events
  
 `python learning  dbcan_summary_stol_pt_named.txt   10 ll_pbd_pgl dbcan2desc_name_full.txt pt_all_settings2acc.txt  config.json   --is_phypat_and_rec --rec_dir pgl_matrix_gain_loss/  --likelihood_params threshold:0.5,mode:gain_loss`
