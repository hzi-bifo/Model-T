# TraitarM
Learning phenotype classification models from protein family phyletic patterns and other kinds of annotations 
# Basic usage
```
traitarm -h #show TraitarM help
traitarm <out_dir> <phenotype_table> <annotation_table>  <phenotype_name> --cpus <#CPUs>
```
Annotation table: a tab separated table of samples (one sample per row) vs. annotation (one feature per columns) as for example produced by ``traitar annotate`` (see example usage; https://github.com/aweimann/traitar-model/blob/master/example/dbcan_annot.txt).

Phenotype table: a tab separated table of samples vs. phenotypes (see example usage; https://github.com/aweimann/traitar-model/blob/master/example/pbd.txt). Use 0 to code the phenotype-negative class, 1 for the phenotype-positive class and ? for missing phenotype labels.

## Including a phylogenetic tree
TraitarM can also use the tree as an additional source of information by reconstructing the evovlutionary history of the genetic features and the phenotype. TraitarM will train two additional models in this case: The phypat+PGL model, which uses both the evolutionary history as well as the observable genotype and phenotype patterns and the PGL model, which is inferred only using the genotype and phenotype gains and losses.

```traitarm <out_dir> 10 <phenotype_table> <annotation_table>  <phenotype_name>   --cpus 10 --tree <tree>``` 

The tree should be provided in Newick format and should have been computed basaed for example on a set of aligned marker genes or core-genes for the input genomic bacterial samples.

## Re-using phenotype models
For each type of phenotype models trained (phypat, phypat+PGL and PGL) there will be a re-usable model available that can be used to predict new samples using Traitar (also see https://github.com/aweimann/traitar). Traitar by default only supports annotation using HMMSEARCH. 
```traitar phenotype <input_dir> <sample_table> <from_nucleotides/from_genes> -p <newly_computed_pt_model> --primary_hmm_db <hmm_db_file>```

If you're using some other kind of annotation for instance SNPs, gene expression data or Roary, you have to ensure that the annotation is generated the same way as the training data was generated e.g. using the same SNP calling pipeline. In this case you need to use Traitar with the from_annotation and the -a option to provide a pre-computed annotation table.
```traitar phenotype <input_dir> <sample_table> from_annotation -a <path_to_annotation> -p <newly_computed_pt_model> --primary_hmm_db <hmm_db_file>```

## Example usage
Construct a plant biomass degradation genotype-phenotype model based on the example data

1. Download amino acid FASTA files from PATRIC database for plant biomass degrading and non-degrading strains
  
  ```
  mkdir <example_dir>/faa
  cd <example_dir>
  cut -f1 sample_table.txt | while read i; do wget -P <example_dir>/faa ftp://ftp.patricbrc.org/patric2/current_release/faa/$i; done
  ```
  
2. Annotate downloaded CDS with CAZy families (using 10 cores)

   ```
   traitar annotate faa sample_table.txt from_genes traitar_out -p ../phypat.tar.gz --primary_hmm_db dbCAN-fam-HMMs.txt -c 10
   ```
3. Train phenotype model 
    ```
    traitarm pbd_out 10 example/pbd.txt example/traitar_out/annotation/dbcan/summary.dat  --feature_mapping  example/dbcan2desc_name_full.txt  10   --do_phypat_only --cpus 10
    ``` 
# Output
cv_acc.txt gives information about the overall performance of the models for each phenotype including the 
TPR = True positive rate, TNR = True negative rate, BACC = balanced accuracy, precision and F1-Score

heatmap.png shows in the upper panel the distribution of the most discriminatory feature found by TraitarM across the samples of the phenotype-positive class (black) and the phenotype-negative class (orange). In the bottom panel you can see the importance of each invidual feature for slightly different version of the classifier. Red indicates a feature relevant for the phenotype-positive class. 

non-zero+weights.txt lists all features, which were found to be relevant for the classification including some information on how good these would do discriminate the phenotype-positive  and negative-samples individually.

roc_curve.png shows different trade-offs of sensitivity vs. specificity that are achievable with the classifier.

miscl.txt provides all samples.txt that were assigned to either the phenotype-positive class although they belonng the phenotype-negative class or the vice versa.
