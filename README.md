DenoVar performs the genome level annotatation for the denovo sequenced variant peptides by searching the amino acid changes in the COSMIC, ClinVar and ICGC datasets. The program also fetchs the protein accession for the variant peptides when the protein accession is not available for the given peptides.    

### System requirments
  - Any operating system with conda enabled
  - Minimum of 16 GB of RAM
  - Minimum of 80GB storage space

## Dependencies
The package requires “conda” to be installed and add the following channels 
```bash
 conda config --add channels defaults
 conda config --add channels bioconda
 conda config --add channels conda-forge
```
ANNOVAR (https://annovar.openbioinformatics.org/en/latest/user-guide/download/) needs to be downloaded and placed in the tools directory (refer to directory_structure.txt). Other dependencies can be installed by running environment.yml
```bash
conda env create --file environment.yml
conda activate DenoVar
```

The backend database is created by merging COSMIC, ClinVar, and ICGC datasets. Due to the large size of the database and the license terms, the backend database is not included in the project. It can be manually created by following steps. 

### Refseq protein and geomic GTF file download
NCBI refseq database can be downloaded from https://www.ncbi.nlm.nih.gov/genome/?term=Humans. In this project latest release (109) has been used which can be downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/. The same link also provides GTF file GCF_000001405.39_GRCh38.p13_genomic.gtf.gz. Add the link and file names in config.yml file. These files can be downloaded using the wget command by the DenoVar tool.  

### ICGC database download
Download the README.txt file from https://dcc.icgc.org/api/v1/download?fn=/release_28/README.txt
It is recommended to use the latest release. Example the above link contains “release_28”. It can have updated version in coming future. 
Set the ICGC README.txt file path in config.yml

### ClinVar database download
ClinVar file name can be copied from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ copy the file name which has got large in file size. For Example, the current release is “clinvar_20220313.vcf.gz” updated on 03/14/2022. Add the file name to config.yml

### COSMIC database download
The COSMIC database requires a user to be registered. After the login, follow the below step
Navigate to menu bar -> Downloads and download the “CosmicGenomeScreensMutantExport.tsv.gz” from the COSMIC Mutation Data (Genome Screens) and set the file name in config.yml 
  
### Running the test dataset
Once the above configuration is set, the pipeline is ready to execute. The pipeline takes two types of input.

 With protein accession: (It is recommended to have "sequence" column name in the input file. In else case, feel free to edit the scripts/denovar.py at line XXX to the name of the peptide sequence column)
 Without protein accession: (Make sure "Sequece" and "Protein_accession" columns exist. Else, it can be edited at scripts/denovar at line XXX and XXX respectively) 
 
Set the protein accession availability status at config.yml. Any hassel in the input file format, please refere the test_datasets. The package is currently set to fetch the sequences / peptides based on the NCBI refseq accessions. The uniport accessions have not yet been implemented.

### Execution
```bash
snakemake --cores 10
```
