# Model-T
Learning phenotype classification models from protein family phyletic patterns and other kinds of annotations 
# Installation
Model-T is available via anaconda / bioconda: Install by 

```conda install model-t``` 

Be sure to add the bioconda channel as described on the bioconda website:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

# Basic usage
```
traitarm -h #show Model-T help
traitarm <out_dir> <phenotype_table> <annotation_table>  <phenotype_name> --cpus <#CPUs>
```
Annotation table: a tab separated table of samples (one sample per row) vs. annotation (one feature per columns) as for example produced by ``traitar annotate`` (see example usage; https://github.com/aweimann/traitar-model/blob/master/example/dbcan_annot.txt).

Phenotype table: a tab separated table of samples vs. phenotypes (see example usage; https://github.com/aweimann/traitar-model/blob/master/example/pbd.txt). Use 0 to code the phenotype-negative class, 1 for the phenotype-positive class and ? for missing phenotype labels.

## Including a phylogenetic tree
Model-T can also use the tree as an additional source of information by reconstructing the evovlutionary history of the genetic features and the phenotype. Model-T will train two additional models in this case: The evo+observed model, which uses both the evolutionary history as well as the observable genotype and phenotype patterns and the evo model, which is inferred only using the genotype and phenotype gains and losses.

```traitarm <out_dir> 10 <phenotype_table> <annotation_table>  <phenotype_name>   --cpus 10 --tree <tree>``` 

The tree should be provided in Newick format and should have been computed based for example on a set of aligned marker genes or core-genes for the input genomic bacterial samples.

## Re-using phenotype models
For each type of phenotype models trained (observed, evo, observed+evo) there will be a re-usable model available that can be used to predict new samples using Traitar (also see https://github.com/aweimann/traitar). Traitar by default only supports annotation using hmmsearch. 
```traitar phenotype <input_dir> <sample_table> <from_nucleotides/from_genes> -p <newly_computed_pt_model> --primary_hmm_db <hmm_db_file>```

If you're using some other kind of annotation for instance SNPs, gene expression data or gene family presence and absence information generated with Roary, you have to ensure that the annotation is generated the same way as the training data was generated e.g. using the same SNP calling pipeline or in case of Roary generating your own set of HMMs and annotating the new genomes. In this case you need to use Traitar with the from_annotation and the -a option to provide a pre-computed annotation table.
```traitar phenotype <input_dir> <sample_table> from_annotation -a <path_to_annotation> -p <newly_computed_pt_model> --primary_hmm_db <hmm_db_file>```

## Example usage
Construct a plant biomass degradation genotype-phenotype model based on the example data

1. Download amino acid FASTA files from PATRIC database for plant biomass degrading and non-degrading strains
  
  ```
  mkdir example/faa
  cd example
  cut -f1 sample_table.txt | while read i; do wget -P <example_dir>/faa ftp://ftp.patricbrc.org/patric2/current_release/faa/$i; done
  ```
  
2. Annotate downloaded CDS with Pfam families (using 10 cores)

   ```
   traitar annotate faa sample_table.txt from_genes traitar_out  -c 10
   ```
3. Train plant biomass degradation Pfam phenotype model 
    ```
    cd ..
    traitarm pbd_out example/pbd.txt example/traitar_out/annotation/dbcan/summary.dat  --feature_mapping  example/dbcan2desc_name_full.txt  10   --do_phypat_only --cpus 10
    ``` 
## Using an HMM database other than Pfam

1. Download HMM database

```
wget http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt.v4
```

2. Link dbCAN families with descriptions creating an empty phenotype model archive file using Traitar

```
traitar new  dbcan2desc_name_full.txt dbcan dbcan_v4
```

3. Annotate with Traitar using the new phenotype model archive using 10 CPUs

``` 
traitar annotate --primary_hmm_db dbCAN-fam-HMMs.txt.v4 -p dbcan_v4.tar.gz faa/ sample_table.txt from_genes traitar_dbcan_v4_out/ --cpus 10 
``` 

# Output
cv_acc.txt gives information about the overall performance of the models for each phenotype including the 
TPR = true positive rate, TNR = true negative rate, BACC = balanced accuracy, precision, F1-Score and AUC = Area under the curve

<pt_model>heatmap.png shows in the upper panel the distribution of the most discriminatory feature found by Model-T across the samples of the phenotype-positive class (black) and the phenotype-negative class (orange). In the bottom panel you can see the importance of each invidual feature for slightly different version of the classifier. Red indicates a feature relevant for the phenotype-positive class. 

<pt_model>_non-zero+weights.txt lists all features, which were found to be relevant for the classification including some information on how good these would do discriminate the phenotype-positive  and negative-samples individually.

<pt_model>_roc_curve.png shows different trade-offs of sensitivity vs. specificity that are achievable with the classifier.

<pt_model>_miscl.txt provides all samples.txt that were assigned to either the phenotype-positive class although they belonng the phenotype-negative class or vice versa.
