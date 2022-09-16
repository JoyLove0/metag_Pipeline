# metag_Pipeline
This pipeline evaluates the quality of genomic data before and after trimming and aligns the quality data to a target genome. This allows for the visualization communities within a sample, which is important to metagenomic analysis. Within the pipeline, several programs are used and therefore must be preloaded for the code to work. The programs used include: FastQC, MultiQC, Trimmomatic, Burrows-Wheeler Aligner (BWA), samtools, and Kraken2. 

# Included Tools and Packages
This tool comes in four parts: (1)metag_Pipeline.py, (2)help_me.ppy, (3)config_file, and (4)sample_config_file
  
# Usage

`python metag_Pipeline.py`
 
# Dependecies
This pipeline requires Python, FASTQC, Trimmomatic, bwa, samtools, and  to be installed in order to run. 

The following is the conda environment recipe that includes these programs and what version this pipeline was developed in. 

| Command | Description |
| --- | --- |
| git status | List all new or modified files |
| git diff | Show file differences that haven't been staged |

Note: channel ordering might be required to download these programs to a new or existing conda environment. 

Along with a configuration file, a help henu program should be included with this Pipeline. This help menu, called Help_Me.py, contains all help menus of all the programs used, recommended parameters, how to overcome common program problems, and additional sources of information. While the information has less detail than this document, the Help_Me file is an incredibly useful tool.

This information is taken directly from the help menu or the manual of each program. 

### Fastqc

FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

	fastqc seqfile1 seqfile2 .. seqfileN
  
  fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
         [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of 
    which will help to identify a different potential type of problem in your
    data.
    
    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.
    
    The options for the program as as follows:

| Options | Description |
| --- | --- |
| `-h --help` | Print this help file and exit |
| `-v --version`| Print the version of the program and exit |
| `-o --outdir` | Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it.  If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.|
| ` --casava` | Files come from raw casava output. Files in the same sample |
|               |    |
|               |     |
|               |        |
                                                          
 BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
                   
### Multiqc
Usage: multiqc [OPTIONS] <analysis directory>

  MultiQC aggregates results from bioinformatics analyses across many samples
  into a single report.

  It searches a given directory for analysis logs and compiles a HTML report.
  It's a general use tool, perfect for summarising the output from numerous
  bioinformatics tools.

  To run, supply with one or more directory to scan for analysis results. To
  run here, use 'multiqc .'

  See http://multiqc.info for more details.

  Author: Phil Ewels (http://phil.ewels.co.uk)

| Options | Description |
| --- | --- |
| `-h --help` | Print this help file and exit |
  
# Trimmomatic
Usage:
  
  PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
  
  SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
 
Parameter Options:
| Options | Description |
| --- | --- |
| `ILLUMINACLIP` | Cut adapter and other illumina-specific sequences from the read. |
| `SLIDINGWINDOW` | Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below threshold.|
| `MAXINFO` | An adaptive quality trimmer which balances read length and error rate to maximise the value of each read. |
| `LEADING` | Cut bases off the start of a read, if below a threshold qualityTRAILING: Cut bases off the end of a read, if below. |
| `CROP` | Cut the read to a specified length by removing bases from the end. |
| `HEADCROP` | Cut the specified number of bases from the start of the read. |
| `MINLEN` | Drop the read if it is below a specified length. |  
| `AVGQUAL` | Drop the read if the average quality is below the specified level. |
| `TOPHRED33` | Convert quality scores to Phred-33. |
| `TOPHRED64` | Convert quality scores to Phred-64. |

# Burrows-Wheeler Alignment

# Samtools
  
# Kraken2
  

# Configuration File
In order for this code to work this file MUST be called "congfig_file". The only things that should be modified are the words AFTER the equals sign. Make sure there is a space before and after the equals sign once you input your data.

 
  
  
  
  
  
