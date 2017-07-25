# TraitarM
Learning phenotype classification models from protein family phyletic patterns and other kinds of annotations 
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
    
