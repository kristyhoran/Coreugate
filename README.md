# COREugate - A Nextflow pipeline for cgMLST
## From reads to cgMLST profile.

This is a simple pipeline that allows the user to input paired-end reads and will output a cgMLST profile for the isolates as well as statistics for allele calling.

1. Assemble - default using Shovill (implementing the latest version of Spades)
2. Call alleles using [chewBBACA](https://github.com/B-UMMI/chewBBACA/wiki).
3. Combine profiles and statisitics for the whole dataset.
### Get a cgMLST scheme

Download a cgMLST scheme from a public repository or use [Coreuscan](https://github.com/kristyhoran/coreuscan) and use chewBBACA to prepare the scheme for allele calling (you may need to install chewBBACA)

```
chewBBACA.py PrepExternalSchema -i pathFolderSchemaFastaFiles/
```
where pathFolderSchemaFastaFiles/ is the path to downloaded scheme.

### Nextflow
Ensure that you have Nextflow installed. Detailed instructions can be found [here](https://www.nextflow.io/docs/latest/getstarted.html)

### Run Coreugate
Once you have a scheme and have installed Nextflow you can run Coreugate

```
REQUIRED

	--fastq	path to reads in the form of path/to/reads/*{1,2},
		where * is a sample name or identifier (note that this will be the identifier that is carried throughout the pipeline for that sample

	--schemaPath path/to/downloaded/scheme

OPTIONAL
	--assembler can choose: skesa, shovill (implementing latest Spades version), Spadesv0.3.7: default shovill 
		IMPORTANT: if you wish to use Spades v0.3.7, add `-profile spades` to the end of the command
	--contigs set the minimum contig size to discard assemblies with less than a set threshold: default no minumimun
	--contig_size the minimum contig size kept for analysis: default 200 bp

```

#### Example

```
nextflow run kristyhoran/Coreugate --fastq 'data/*_{1,2}*' --schemaPath 'schemaPath/' --contigs 85
```

### Things to note
* It is improtant to wrap the pipeline arguments in single quotes, unless they are supposed to be passed as a number as seen above for `--contigs 85`. Nextflow arguments can also be included, these do not need to be wrapped in quotes. An example of this can be seen above, when adding `-profile spades` to the command. `spades` does not need to be in quotes.

* If you interrupt the pipeline or want to go back and add isolates later on the flag `-resume` may be addedd to the end of the command, providing you run it from within the original folder.

* When setting up the file structure for a project, it is advisable to have reads named with the a meaningful name for later analysis, the sample name is used throughout to identify the samples at different steps.

* Output of each step will be stored in a folder created by the pipeline called `results`.

* The `work` folder will contain the singularity containers used in the pipeline and temporary files that may be used to troubleshoot if something fails, it may be removed at the completion of the pipeline if you desire.

### Limitations of the pipeline
* At this stage, Coreugate is only able to take reads as input, not assemblies.
* Coreugate is only able to work with pre-exisiting schemas that have been prep as described above, to derive profiles for isolates.
* Possibly more, I just haven't found them yet!!

### Features to come
* Optionally use reads and/or assemblies
* Download and prep schemas in pipeline
* Call alleles without a pre-existing schema

* Will take requests!


