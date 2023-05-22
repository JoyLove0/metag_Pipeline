# metag_Pipeline
This pipeline evaluates the quality of genomic data before and after trimming and aligns the quality data to a target genome. The reads that don't align to the target undergo metagenomic analysis that identify the unmapped communities. Within the pipeline, several programs are used and therefore must be preloaded for the code to work. The programs used include: FastQC, MultiQC, Trimmomatic, Burrows-Wheeler Aligner (BWA), samtools, bedtools, seqtk, and Kraken2.

See the figure below for a summary:

![Metagenomics Pipeline (3)](https://user-images.githubusercontent.com/108104001/235510791-75a25ad2-e2a7-4805-b2fe-9acfe8137e54.png)


# Usage

`python metag_Pipeline.py`

# Included Tools and Packages
This tool comes in several parts: 
1. metag_Pipeline.py 
2. help_me.py
3. parameter_file
4. sample.parameter.file
5. environment.yml
6. R scripts for further data processing
 
# Dependecies
This pipeline requires Python, FASTQC, Trimmomatic, bwa, samtools, bedtools, seqtk, and kraken2 to be installed in order to run. The environment.yml file contains the conda environment recipe that includes these programs. The help_me.py file contain the option and/or of these programs.

Note that this pipeline was written in the Python version 3.6.10. For for more information about the versions of all of these dependecies check the environment.yml.

For more information on each program, see the below:

| Program       | Manuel and/or Useful Link   |
| ------------- | ------------- |
| FastQC        | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ or https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf |
| MultiQC       | https://multiqc.info/docs/                                                                                       |
| Trimmomatic   | http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf                       |
| BWA           | https://bio-bwa.sourceforge.net/bwa.shtml                                                                        |
| samtools      | https://www.htslib.org/doc/samtools-sort.html                                                                    |
| Kraken2       | https://github.com/DerrickWood/kraken2/wiki/Manual                                                               |
| R             | https://cran.r-project.org/manuals.html.                                                                         |

### Genome Indexing
Important: the user must index the reference genome before running metag_Pipeline.py. This will take several hour. Use the command below. 

`bwa index -a bwtsw ref_genome.fa`

Here is a link to the manuel for more information: https://bio-bwa.sourceforge.net/bwa.shtml

### Parameter File
In order for this code to work this file MUST be called "parameter_file". The only things that should be modified are the words AFTER the equals sign. Make sure there is a space before and after the equals sign once you input your data. The following table gives a brief explanation of each line in the parameter file.

| Option        | Details    |
| ------------- | ------------- |
| threads          | The number the user chooses will be use for all steps. 24 is recommended. |
| steps            | This option includes each "step" or option used for FASTQC. ILLUMINACLIP:adapters:2:30:10 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:100 is recommended. Note the "adapters" in this examples. This file can be named anything as longs as it specified in the parameter file. This adapter file must be in the same working directory and has to be in fasta format.|
| stem             | The is the working directory. Note that the "stem" in the parameter file must be the directory that has your raw data.   | 
| filenames        | This is the name of the text file including a list of all the data files to be evaluated. Note, only the foward read files should be in the list. Use \*R1\*.fq.gz" as a name template. Input files can be compressed. Also, note that this file must be in the stem directory.|   
| samtools_quality | This gives you the option to delete mapped and unmapped reads that fall below a set quaility score. It is recommended that you use 20. |
| genome_path      | This is the path to the directory where the reference genome is. This is needed for alignment. An example of this is: /home/dbs/mouse-genome/GCF_000001635.27_GRCm39_genomic.fna |
| database_path    | This is the path to the directory where the database is. This needed for idenifying communities with KRAKEN2. |  
| output_directory_name | This is the name of the directory where all the final ouput will be in. This directory should contain the directories for Kraken, Multiqc, and Alignments. |

See the sample parameter file for more examples. 

### Addition Data Processing using R script Provided
Note that the R script that are provided were created using R version 4.1.2. 

While the Kraken report are very useful, there are several more ways to make the information more accessible. The R scripts provided will make a compiled chart of the Kraken reports of all of your samples, a depth bar graph, and Rarifaction curves. The latter two will be generated as pdf documnets. The complied kraken report will be generated as a csv document that can be opened in excel. 

In order for the R scripts to run properly, the file data_processing.R and the function.R must both be dowloaded in the Kraken directory. 
