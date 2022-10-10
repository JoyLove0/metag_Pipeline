# metag_Pipeline
This pipeline evaluates the quality of genomic data before and after trimming and aligns the quality data to a target genome. The reads that don't align to the target undergo metagenomic analysis that idetify the unmapped communities. Within the pipeline, several programs are used and therefore must be preloaded for the code to work. The programs used include: FastQC, MultiQC, Trimmomatic, Burrows-Wheeler Aligner (BWA), samtools, bedtools, seqtk, and Kraken2.

See the figure below for a summary:

![Metagenomics Pipeline](https://user-images.githubusercontent.com/108104001/194805896-1e1e923c-bb15-4c56-aac0-5d702cf97b00.jpg)

# Usage

`python metag_Pipeline.py`

# Included Tools and Packages
This tool comes in four parts: 
	(1)metag_Pipeline.py 
	(2)help_me.py
	(3)parameter_file
	(4)sample.parameter.file
 
# Dependecies
This pipeline requires Python, FASTQC, Trimmomatic, bwa, samtools, bedtools, seqtk, and kraken2 to be installed in order to run. The environment.yml file contains the conda environment recipe that includes these programs. The help_me.py file contain the option and/or of these programs.

### Genome Indexing
Important: the user must index the reference genome before running metag_Pipeline.py. Use the command below. 

`bwa index -a bwtsw ref_genome.fa`

Here is a link to the manuel for more information: https://bio-bwa.sourceforge.net/bwa.shtml

### Parameter File
In order for this code to work this file MUST be called "parameter_file". The only things that should be modified are the words AFTER the equals sign. Make sure there is a space before and after the equals sign once you input your data.
