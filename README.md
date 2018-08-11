# COREugate - A snakemake pipeline for cgMLST
## From reads to cgMLST profile.

This is a simple pipeline that allows the user to input paired-end reads and will output a cgMLST profile for the isolates as well as statistics for allele calling.

1. Assemble - default using Shovill (implementing the latest version of Spades)
2. Call alleles using [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki).
3. Combine profiles and statisitics for the whole dataset.
### Get a cgMLST scheme
Download a cgMLST scheme from a public repository or use [Coreuscan](https://github.com/kristyhoran/coreuscan).

### Dependencies
```
Python <3.6
Biopython
Snakemake
```

### Biopython
Biopython is used here for quality control of assemblies, information about biopython can be found [here](https://biopython.org)
```
pip3 install biopython
```

### Snakemake
Ensure that you have Snakemake installed. Detailed instructions can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

```
pip3 install snakemake
```

### Run Coreugate

*1. Get Coreugate*
```
$ mkdir MyProjectDir
$ cd MyProjectDir
$ git clone https://github.com/kristyhoran/Coreugate.git
```
*2. Setup*
Input to Coreugate is a tab-delimited file (default name isolates.tab).
```
isolate_name	path/to/reads/R1.fq.gz	path/to/reads/R2.fq.gz
```
*3. Run Coreugate*
Basic
```
snakemake --config schemaPath=path/to/Schema --use-singularity
```
If you do not wish to use singualrity containers to run assembly and chewBBACA - proceed at own risk! and ensure that you have one of shovill, skesa or spades 3.7 and chewBBACA > 2.0.12 installed.

Customised - run Coreugate, assembling with skesa and only using samples with > 85 contigs for input into chewBBACA
```
snakemake --config schemaPath=path/to/Schema assembler=skesa min_contigs=85 --use-singularity
```

Other parameters that can be user-defined are
```
min_contig_size default=500
prepCPU default=36
chewieCPU default=36
isolate_file default=isolates.tab
```


### Things to note
* chewBBACA will update your schema with newly identified alleles for loci. Because of this chewBBACA currently can not be run in parallel as such COREugate limits the chewBBACA step to run consecutively
* If you interrupt the pipeline or want to go back and add isolates simply add new isolates to the end of the isolates.tab file
* Output of each step will be stored in a folder created by the pipeline called `coreugate`.

### Limitations of the pipeline
* At this stage, COREugate only supports reads as input. Future plans include being able to use a mix of reads and assemblies.
* Coreugate is only able to work with pre-exisiting schemas that have been prep as described above, to derive profiles for isolates.
* Possibly more, I just haven't found them yet!!

### Features to come
* Optionally use reads and/or assemblies
* Download schemas
* Call alleles without a pre-existing schema

* Will take requests!


