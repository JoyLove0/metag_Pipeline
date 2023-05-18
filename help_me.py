#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 19:26:10 2022

@author: joy

cd Desktop/EatsWorms/Spyder/Raw_Data_copy
"""

def help_me():
    print("[1] FASTQC")
    print("[2] MULTIQC")
    print("[3] TRIMMOMATIC")
    print("[4] BWA INDEX")
    print("[5] BWA MEM")
    print("[6] SAMTOOLS SORT")
    print("[7] BEDTOOLS")
    print("[8] SEQTK")
    print("[9] KRAKEN2")
    print("[0] EXIT HELP MENU.")
    
###########################################################################
    
while True:
    help_me()
    program = int(input("Enter the program you are having trouble with >> "))
    if program == 0:
        break
    if program == 1:
        print("-------------------")
        print("Link to manuel: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf\n")
        print("COMMAND: fastqc R1_file.fastq R2_file.fastq\n")
        print("OPTIONS:\n",
              "\t -h --help:         Print help file and exit.\n",
              "\t -v --version:      Print the version of the program and exit.\n",
              "\t -o --outdir:       Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it. If this option\n",
              "\t                    is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.\n",
              "\t --casava:          Files come from raw casava output. Files in the same sample group (differing only by the group number) will be analysed as a set rather than\n",
              "\t                    individually. Sequences with  the  filter flag set in the header will be excluded from the analysis. Files must have the same names given to them by\n", 
              "\t                    casava (including being gzipped and ending with .gz) otherwise they won't be grouped together correctly. \n",
              "\t --nano:            Files come from nanopore sequences and are in fast5 format. In this mode you can pass in directories to process and the program will take in all fast5\n",
              "\t                    files within those directories and produce a single output file from the sequences found in all files.\n", 
              "\t --nofilter:        If running with --casava then don't remove read flagged by casava as poor quality when performing the QC analysis.\n",
              "\t --extract:         If set then the zipped output file will be uncompressed in the same directory after it has been created.  By default this option will be set if fastqc\n",
              "\t                    is run in non-interactive mode.\n",
              "\t -j --java:         Provides the full path to the java binary you want to use to launch fastqc. If not supplied then java is assumed to be inyour path.\n",
              "\t --noextract:       Do not uncompress the output file after creating it. You should set this option if you do not wish to uncompress the output when running in \n",
              "\t                    non-interactive mode.\n",
              "\t --nogroup:         Disable grouping of bases for reads >50bp. All reports will show data for every base in the read. WARNING: Using this option will cause fastqc to crash \n",
              "\t                    and burn if you use it on really long reads, and your plots may end up a ridiculous size. You have been warned!\n",
              "\t --min_length:      Sets an artificial lower limit on the length of the sequence to be shown in the report. As long as you set this to a value greater or equal to your\n",
              "\t                    longest read length then this will be the sequence length used to create your read groups. This can be useful for making directly comaparable statistics\n",
              "\t                    from datasets with somewhat variable read lengths.\n",
              "\t -f --format:       Bypasses the normal sequence file format detection and forces the program to use the specified format. Valid formats are bam,sam,bam_mapped,sam_mapped\n,"
              "\t                    and fastq.\n",
              "\t -t --threads:      Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads\n",
              "\t                    than your available memory will cope with, and not more than 6 threads on a 32 bit machine.\n",
              "\t -c --contaminants: Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against. The file must contain sets of named\n",
              "\t                    contaminants in the form name[tab]sequence. Lines prefixed with a hash will be ignored.\n",
              "\t -a --adapters:     Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library. The file must contain\n",
              "\t                    sets of named adapters in the form name[tab]sequence. Lines prefixed with a hash will be ignored.\n",
              "\t -l --limits:       Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules. This file\n",
              "\t                    can also be used to selectively remove some modules from the output all together. The format needs to mirror the default limits.txt file found in the\n",
              "\t                    Configuration folder.\n",
              "\t -k --kmers:        Specifies the length of Kmer when looking in the Kmer content module. Specified Kmer length must be between 2 and 10. Default length is 5.\n"
              "\t -q --quiet:        Supress all progress messages on stdout and only report errors.\n",
              "\t -d --dir:          Selects a directory to be used for temporary files written when generating report images. Defaults to system temp directory if not specified.\n")
        print("-------------------")
    elif program == 2:
        print("-------------------")
        print("Link to manuel: https://multiqc.info/docs/ \n")
        print("COMMAND: multiqc -o output_dir --title file_title multiqc_path \n")
        print("OPTIONS:\n",
              "\t -f, --force:                   Overwrite any existing reports.\n",
              "\t -d, --dirs:                    Prepend directory to sample names.\n",
              "\t -dd, --dirs-depth INTEGER:     Prepend [INT] directories to sample names.\n",
              "\t -s, --fullnames:               Do not clean the sample names (leave as full file name).\n",
              "\t -i, --title TEXT:              Report title. Printed as page header, used for filename if not otherwise specified.\n",
              "\t -b, --comment TEXT:            Custom comment, will be printed at the top of the report.\n",
              "\t -n, --filename TEXT:           Report filename. Use 'stdout' to print to standard out.\n",
              "\t -o, --outdir TEXT:             Create report in the specified output directory.\n",
              "\t -t, --template:                [default|default_dev|gathered|geo|sections|simple] Report template to use.\n",
              "\t --tag TEXT:                    Use only modules which tagged with this keyword, eg. RNA.\n",
              "\t --view-tags, --view_tags:      View the available tags and which modules they load.\n",
              "\t -x, --ignore TEXT:             Ignore analysis files (glob expression).\n",
              "\t --ignore-samples TEXT:         Ignore sample names (glob expression).\n",
              "\t --ignore-symlinks:             Ignore symlinked directories and files.\n",
              "\t --replace-names PATH:          Path to TSV file to rename sample names during report generation.\n",
              "\t --sample-names PATH:           Path to TSV file containing alternative sample names for renaming buttons in the report.\n",
              "\t --sample-filters PATH:         Path to TSV file containing show/hide patterns for the report.\n",
              "\t -l, --file-list:               Supply a file containing a list of file paths to be searched, one per row.\n",
              "\t -e, --exclude [module name]:   Do not use this module. Can specify multiple times.\n",
              "\t -m, --module [module name]:    Use only this module. Can specify multiple times.\n",
              "\t --data-dir:                    Force the parsed data directory to be created.\n",
              "\t --no-data-dir:                 Prevent the parsed data directory from being created.\n",
              "\t -k, --data-format:             [tsv|json|yaml] Output parsed data in a different format. Default: tsv\n",
              "\t -z, --zip-data-dir:            Compress the data directory.\n",
              "\t --no-report:                   Do not generate a report, only export data and plots.\n",
              "\t -p, --export:                  Export plots as static images in addition to the report.\n",
              "\t -fp, --flat:                   Use only flat plots (static images).\n",
              "\t -ip, --interactive:            Use only interactive plots (HighCharts Javascript).\n",
              "\t --lint:                        Use strict linting (validation) to help code development.\n",
              "\t --pdf:                         Creates PDF report with 'simple' template. Requires Pandoc to be installed.\n",
              "\t --no-megaqc-upload:            Don't upload generated report to MegaQC, even if MegaQC options are found.\n",
              "\t -c, --config PATH:             Specific config file to load, after those in MultiQC dir / home dir / working dir.\n",
              "\t --cl-config, --cl_config TEXT: Specify MultiQC config YAML on the commandline.\n",
              "\t -v, --verbose:                 Increase output verbosity.\n",
              "\t -q, --quiet:                   Only show log warnings. \n",
              "\t --profile-runtime:             Add analysis of how long MultiQC takes to run to the report.\n",
              "\t --no-ansi:                     Disable coloured log output.\n",
              "\t --custom-css-file PATH:        Custom CSS file to add to the final report.\n",
              "\t --version:                     Show the version and exit.\n",
              "\t -h, --help:                    Show help message and exit.\n",)
        print("-------------------")
    elif program == 3:
        print("-------------------")
        print("Link to manuel: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf\n")
        print("USAGE:\n",
              "PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> \n <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...\n \t                  \n or \n\n -version")
        print("\nPARAMETERS OPTIONS:\n", 
              "\t ILLUMINACLIP:  Cut adapter and other illumina-specific sequences from the read.\n",
              "\t SLIDINGWINDOW: Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a\n",
              "\t                threshold.\n",
              "\t MAXINFO:       An adaptive quality trimmer which balances read length and error rate to maximise the value of each read.\n",
              "\t LEADING:       Cut bases off the start of a read, if below a threshold quality""TRAILING: Cut bases off the end of a read, if below\n",
              "\t                a threshold quality.\n",
              "\t CROP:          Cut the read to a specified length by removing bases from the end HEADCROP: Cut the specified number of bases from the start of the read MINLEN: Drop the\n",
              "\t                read if it is below a specified length.\n",
              "\t AVGQUAL:       Drop the read if the average quality is below the specified level TOPHRED33: Convert quality scores to Phred-33.\n",
              "\t TOPHRED64:     Convert quality scores to Phred-64.\n")
        print("RECOMMENDED PARAMETERS: ILLUMINACLIP:adapters:2:30:10 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:75")
        print("-------------------")
    elif program == 4:
        print("-------------------")
        print("Link to manuel: http://bio-bwa.sourceforge.net/bwa.shtml\n")
        print("COMMAND: bwa index [-p prefix bwtsw|is] [-a algoType] <in.db.fasta>\n")
        print("OPTIONS:\n",
             "-p STR: Prefix of the output database [same as db filename]\n",
             "-a STR: Algorithm for constructing BWT index and can have two values:\n",
                 "\t is:    Does not work for long genomes. IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is\n",
                 "\t        moderately fast, but does not work with database larger than 2GB.\n",
                 "\t bwtsw: Does not work for short genomes. Algorithm implemented in BWT-SW. This method works with the whole human genome.\n")
        print("RECOMMENDED PARAMETERS: -a bwtsw")
        print("-------------------")
    elif program == 5:
        print("-------------------")
        print("Link to manuel: http://bio-bwa.sourceforge.net/bwa.shtml\n")
        print("COMMAND: bwa mem ref.fa read1.fq read2.fq > aln-pe.sam \n",
              "\nCOMMAND (with options):\n",
              "bwa mem [-aCHMpP] [-t nThreads] [-k minSeedLen] [-w bandWidth] [-d zDropoff] [-r seedSplitRatio] [-c maxOcc] [-A matchScore] [-B mmPenalty] [-O gapOpenPen]\n",
              "[-E gapExtPen] [-L clipPen] [-U unpairPen] [-R RGline] [-v verboseLevel] db.prefix reads.fq [mates.fq]\n")
        print("OPTIONS:\n",
              "\t -t INT:   Number of threads. \n",
              "\t -k INT:   Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. \n",
              "\t -w INT:   Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not\n",
              "\t           solely determined by this option. \n" ,
              "\t -d INT:   Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are\n",
              "\t           the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t\n",
              "\t           penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good\n",
              "\t           alignment. \n",
              "\t -r FLOAT: Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which\n",
              "\t           leads to faster alignment speed but lower accuracy. \n",
              "\t -c INT:   Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. \n",
              "\t -P:       In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair. \n",
              "\t -A INT:   Matching score. \n",
              "\t -B INT:   Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. \n",
              "\t -O INT:   Gap open penalty. \n",
              "\t -E INT:   Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). \n",
              "\t -L INT:   Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score\n",
              "\t           minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. \n",
              "\t -U INT:   Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty.\n",
              "\t           It compares these two scores to determine whether we should force pairing. \n",
              "\t -p:       Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details. \n",
              "\t -R STR:   Complete read group header line. Backslash t can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read\n",
              "\t           in the output. An example is ’@RG\tID:foo\tSM:bar’. \n",
              "\t -T INT:   Don’t output alignment with score lower than INT. This option only affects output. \n",
              "\t -a:       Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments. \n",
              "\t -C:       Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the\n",
              "\t           FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output. \n",
              "\t -H:       Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences. \n",
              "\t -M:       Mark shorter split hits as secondary (for Picard compatibility). \n",
              "\t -v INT:   Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value 0 for disabling all the output to stderr; 1\n",
              "\t           for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM\n")
        print("PARAMETERS USED: -M -t 32 ")
        print("-------------------")
    elif program == 6:
        print("-------------------")
        print("Link to manuel: https://www.htslib.org/doc/samtools-sort.html\n")
        print("COMMAND: samtools sort [-l level] [-u] [-m maxMem] [-o out.bam] [-O format] [-M] [-K kmerLen] [-n] [-t tag] [-T tmpprefix] [-@ threads] [in.sam|in.bam|in.cram]\n")
        print("NOTE: If you are having trouble with dowloading samtools to a conda environment, explicit channel ordering may be needed. Here is a link for more info/help:\n",
              "https://github.com/bioconda/bioconda-recipes/issues/22351. Here is an example from that website: conda install -n samtools_foo samtools=1.10 -c bioconda -c conda-forge -c defaults\n")
        print("OPTIONS:\n",
              "\t -K INT:        Sets the kmer size to be used in the -M option. \n",
              "\t -l INT:        Set the desired compression level for the final output file, ranging from 0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best\n",
              "\t                compression but slowest to write), similarly to gzip(1)'s compression level setting. \n",
              "\t -u:            Set the compression level to 0, for uncompressed output. This is a synonym for -l 0. \n",
              "\t -m INT:        Approximately the maximum required memory per thread, specified either in bytes or with a K, M, or G suffix. To prevent sort from creating a huge number\n",
              "\t                of temporary files, it enforces a minimum value of 1M for this setting. \n",
              "\t -M:            Sort unmapped reads (those in chromosome '*') by their sequence minimiser (Schleimer et al., 2003; Roberts et al., 2004), also reverse complementing\n",
              "\t                as appropriate. This has the effect of collating some similar data together, improving the compressibility of the unmapped sequence. The minimiser kmer\n",
              "\t                size is adjusted using the -K option. Note data compressed in this manner may need to be name collated prior to conversion back to fastq. Mapped sequences\n",
              "\t                are sorted by chromosome and position.\n",
              "\t -n:            Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates. \n",
              "\t -t TAG:        Sort first by the value in the alignment tag TAG, then by position or name (if also using -n). \n",
              "\t -o FILE:       Write the final sorted output to FILE, rather than to standard output. \n",
              "\t -O FORMAT:     Write the final output as sam, bam, or cram. By default, samtools tries to select a format based on the -o filename extension; if output is to standard\n",
              "\t                output or no format can be deduced, bam is selected.\n",
              "\t -T PREFIX:     Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique\n",
              "\t                to this invocation of the sort command. By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or if output is to\n",
              "\t                standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam. \n",
              "\t -@ INT:        Set number of sorting and compression threads. By default, operation is single-threaded. \n",
              "\t --no-PG:       Do not add a @PG line to the header of the output file.\n",
              "\nOrdering Rules: If option -t is in use, records are first sorted by the value of the given alignment tag, and then by position or name (if using -n). For example, “-t RG” will make read group the primary sort key. The rules for ordering by tag are:\n",
              "\t -Records that do not have the tag are sorted before ones that do. \n",
              "\t -If the types of the tags are different, they will be sorted so that single character tags (type A) come before array tags (type B), then\n",
              "\t string tags (types H and Z), then numeric tags (types f and i).\n",
              "\t -Numeric tags (types f and i) are compared by value. Note that comparisons of floating-point values are subject to issues of rounding and precision.\n",
              "\t -String tags (types H and Z) are compared based on the binary contents of the tag using the C strcmp(3) function.\n",
              "\t -Character tags (type A) are compared by binary character value.\n",
              "\t -No attempt is made to compare tags of other types — notably type B array values will not be compared.\n")       
        print("PARAMETERS USED: samtools sort -@64 -o sorted_alignment.bam")
        print("-------------------")
    elif program == 7:
        print("-------------------")
        print("Link to Manuel: https://bedtools.readthedocs.io/en/latest/")
        print("PARAMETERS USED:\n",
              "bedtools bamtofastq -i unmapped.bam -fq %sunmapped.fq")
        print("-------------------")
    elif program == 8:
        print("-------------------")
        print("Link to Manuel: https://github.com/lh3/seqtk \n")
        print("PARAMETERS USED:\n",
              "seqtk seq -1 unmapped_reads.fq > unmapped_reads_R1.fq \n",
              "seqtk seq -2 unmapped_reads.fq > unmapped_reads_R2.fq")
        print("-------------------")
    elif program == 9:
        print("-------------------")
        print("Link to manuel: https://github.com/DerrickWood/kraken2/wiki/Manual \n")
        print("If you are having trouble with installing kraken2 to your conda environment, explicit channel ordering may be needed. Try:\n",
              "conda install kraken2 --channel conda-forge --channel bioconda -c defaults\n")
        print("OPTIONS:\n",
            "\t --db NAME:                   Name for Kraken 2 DB (default: none).\n",
            "\t --threads NUM:               Number of threads (default: 1).\n",
            "\t --quick:                     Quick operation (use first hit or hits).\n",
            "\t --unclassified-out FILENAME: Print unclassified sequences to filename.\n",
            "\t --classified-out FILENAME:   Print classified sequences to filename.\n",
            "\t --output FILENAME:           Print output to filename (default: stdout); '-' will suppress normal output\n",
            "\t --confidence FLOAT:          Confidence score threshold (default: 0.0); must be in [0, 1].\n",
            "\t --minimum-base-quality NUM:  Minimum base quality used in classification (def: 0, only effective with FASTQ input).\n",
            "\t --report FILENAME:           Print a report with aggregrate counts/clade to file.\n",
            "\t --use-mpa-style:             With --report, format report output like Kraken 1's kraken-mpa-report.\n",
            "\t --report-zero-counts:        With --report, report counts for ALL taxa, even if counts are zero.\n",
            "\t --report-minimizer-data:     With --report , report minimizer and distinct minimizer count information in addition to normal Kraken report.\n",
            "\t --memory-mapping:            Avoids loading database into RAM.\n",
            "\t --paired:                    The filenames provided have paired-end reads.\n",
            "\t --use-names:                 Print scientific names instead of just taxids\n",
            "\t --gzip-compressed:           Input files are compressed with gzip.\n",
            "\t --bzip2-compressed:          Input files are compressed with bzip2\n",
            "\t --minimum-hit-groups NUM:    Minimum number of hit groups (overlapping k-mers sharing the same minimizer) needed to make a call (default: 2).\n",
            "\t --help:                      Print this message\n",
            "\t --version:                   Print version information\n",
            "If none of the *-compressed flags are specified, and the filename provided is a regular file, automatic format detection is attempted.\n")
        print("PARAMETERS USED:\n",
              "kraken2 --db DATABASE.fna --paired R1_file R2_files --threads 10 --output /KrakenDirectory/Kraken_Outputfile --report /KrakenDirectory/Kraken_Report --use-name")
        print("-------------------")
    else:
        print("Invalid option. Try again.") 
        
        
print('"Help Me" program terminated')
