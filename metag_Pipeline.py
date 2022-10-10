#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Joy
tile : Metagenomic Pipeline
"""
#################################### SET-UP ###################################

### Built-ins
import subprocess
import os

### Program Functions
## FASTQC function
#This code creates a function that needs the path and output directory to do fastqc
def do_FAStQC(fastq_path, output_dir, options=['-o',]):
    command = "fastqc {} {} {}".format(' '.join(options), output_dir, fastq_path)
    subprocess.check_call(command, shell=True)
## Multiqc function
#This code creates a function that needs the path, output directory, and the title 
#of the created file to do multiqc
def do_Multiqc(output_dir, file_title, multiqc_path):
    command = "multiqc -o " + output_dir + " --title " + file_title + " " + multiqc_path
    subprocess.check_call(command, shell=True)
## Pause function
#This code creat a function that will pause the pipeline so the user can verify data
def pause():
    input("Press the <ENTER> key to continue...")

### Paramter file 
#This code convert the the parameter file into a dictionary
param = {}
file1 = open("parameter_file", "r+")
for line in file1:
    print(line)
    (key, value) = line.strip().split('=')
    param[key.strip()] = value.strip()
file1.close()
#Double check dictionary
print("Parameters Dictionary:", param)

### Variables
## Paramter Variblaes
#This code uses the param dictionary to create variables
QC_threads = param["QC_threads"]
steps = param["steps"]
stem = param["stem"] #has to be Raw_data directory
filenames = param["filenames"]
trimmed_filenames = param["trimmed_filenames"]
ali_threads = param["ali_threads"]
genome_path = param["genome_path"]
samtools_threads = param["samtools_threads"]
database_path = param["database_path"]
kraken_threads = param["kraken_threads"]
#Double check variables
print("-------------------------------------------CHECK YOUR VARIABLES--------------------------------------------")
print("Thread for Fastqc and Multiqc:", QC_threads,
      "Steps", steps,
      "The stem:", stem,
      "Filenames:", filenames,
      "Trimmed Filenames:", trimmed_filenames,
      "Threads for bwa:", ali_threads,
      "Pathway to Reference Genome:", genome_path,
      "Threads for samtools:", samtools_threads,
      "Pathway to Database:", database_path)
## Directory Variables
#This code creates the nessecary directories if they don't already exist
checkpoint = (os.path.exists(stem + "/PE"),
              os.path.exists(stem + "/U"), 
              os.path.exists(stem + "/Pre_Trim"), 
              os.path.exists(stem + "/PE/Post_Trim"), 
              os.path.exists(stem + "/Multiqc_dir"), 
              os.path.exists(stem + "/Kraken"), 
              os.path.exists(stem + "/Alignment"))
if checkpoint[0] == True:
    print("PE directory already exsists.")
elif checkpoint[0] == False:
    os.mkdir(stem + "/PE")
if checkpoint[1] == True:
    print("U directory already exsists.")
elif checkpoint[1] == False:
    os.mkdir(stem + "/U")
if checkpoint[2] == True:
    print("Pre_Trim directory already exsists.")
elif checkpoint[2] == False:
    os.mkdir(stem + "/Pre_Trim")
if checkpoint[3] == True:
    print("Post_Trim directory already exsists.")
elif checkpoint[3] == False:
    os.mkdir(stem + "/PE/Post_Trim")
if checkpoint[4] == True:
    print("Multiqc directory already exsists.")
elif checkpoint[4] == False:
    os.mkdir(stem + "/Multiqc_dir")
if checkpoint[5] == True:
    print("Kraken directory already exsists.")
elif checkpoint[5] == False:
    os.mkdir(stem + "/Kraken")
if checkpoint[6] == True:
        print("Alignment directory already exsists.")
elif checkpoint[6] == False:
        os.mkdir(stem + "/Alignment")
#This code creates variables for directories
PE = stem + "/PE"
U = stem + "/U"
pre_t = stem + "/Pre_Trim"
post_t = PE + "/Post_Trim"
Multiqc_dir = stem + "/Multiqc_dir"
Kraken_dir = stem + "/Kraken"
Alignment_dir = stem + "/Alignment"
print("Paired-end directory:", PE, 
      "Unpaired directory:", U, 
      "Pre-Trimming directory:", pre_t, 
      "Post-Trimming:", post_t,
      "Multiqc directory", Multiqc_dir, 
      "Kraken Directory", Kraken_dir,
      "Alignment Directory", Alignment_dir)
#This code double checks working directories
wd = os.getcwd() #Get the current working directory
print("Current working directory:", wd) #Print the current working directory
print("-----------------------------------------------------------------------------------------------------------")

################################### PIPELINE ##################################

################################ QUALITY CONTROL ##############################

### FASTQC before trimming 
#This code does fastqc on a list of file names
file = open(filenames)
for line in file:
    x = line.replace("R1", "R2")
    do_FAStQC(line, pre_t)
    do_FAStQC(x, pre_t)
print("Done with first round of FASTQC")

### Multiqc before Trimming
#This code compiles FASTQC data using Multiqc 
do_Multiqc(Multiqc_dir, "Pre_Trim", stem)

### Trimming with Trimmomatic 
#This code split names around R1 and runs Trimmomatic
file = open(filenames)
for line in file:
    a = line.split("R1")[0].strip()
    b = line.split("R1")[-1].strip()
    trim_cmd = "trimmomatic PE -threads %s -phred33 %sR1%s %sR2%s %s/%strimmed_R1%s %s/%strimmed_U1%s %s/%strimmed_R2%s %s/%strimmed_U2%s %s" % (QC_threads, a, b, a, b, PE, a, b, U, a, b, PE, a, b, U, a, b, steps)
    print(a, b, trim_cmd)
    subprocess.call([trim_cmd], shell=True)
print("Done with trimming.")

### Manipulating the location of trimmed files 
#This code makes a list of the new trimmed files for final FASTQC
os.chdir(PE) #Change the current working directory
pwd = os.getcwd() #Get the current working directory
print("Current working directory:", pwd) #Print the current working directory
make_file = "ls *R1* > trimmed.file.names" 
subprocess.check_output(make_file, shell=True) #Making file from trimmed names

### FASTQC after Trimming
#This code does fastqc on a list of file names
file = open("trimmed.file.names")
for line in file:
    do_FAStQC(line, post_t)
    do_FAStQC(x, post_t)
print("Done with final round of FASTQC")

### Multiqc after Trimming
#This code compiles FASTQC data using Multiqc 
do_Multiqc(Multiqc_dir, "Post_Trim", PE)

### PAUSE pipeline to verify trimming was completed correctly
print("!!! The program has been paused. Please check the quality of your data by opening a new terminal tab and typing: open post_trim_Multiqc_multiqc_report.html. Come back to this tab to continue running this pipeline OR edit the 'steps' in the congfig_file.")
pause()

########################## ALIGNMENT AND EXTRACTION ###########################

### Aligning with Burrows-Wheeler Alignment Tools
#This code puts user in the Paired end directory
os.chdir(PE) #Change the current working directory
pwd = os.getcwd() #Get the current working directory
print("Current working directory:", pwd) #Print the current working directory
#This code aligns trimmed reads to a provided genome using bwa mem
file = open(trimmed_filenames)
for line in file:
    a = line.split("R1")[0].strip()
    b = line.split("R1")[-1].strip()
#This code extracts the unmapped and mapped alignment files using samtools
    unmapped_samtools_sort = "bwa mem -M -t %s %s %sR1%s %sR2%s | samtools view -@%s -b -f 4 -o %sunmapped%s.bam" % (QC_threads, genome_path, a, b, a, b, samtools_threads, a, b)
    subprocess.call([unmapped_samtools_sort], shell=True)
    mapped_samtools_sort = "bwa mem -M -t %s %s %sR1%s %sR2%s | samtools view -@%s -b -F 4 -o %smapped%s.bam" % (QC_threads, genome_path, a, b, a, b, samtools_threads, a, b)
#This code converts the unmapped alignment files to fastqc files
    bam_to_fastqc = "bedtools bamtofastq -i %sunmapped%s.bam -fq %sunmapped_R1%s.fq" % (a, b, a, b)
    subprocess.call([bam_to_fastqc], shell=True)
    foward_unmapped = "seqtk seq -1 %s_unmapped%s.fq > %sunmapped_R1%s.fq" % (a, b, a, b)
    subprocess.call([foward_unmapped], shell=True)
    reverse_unmapped = "seqtk seq -2 %s_unmapped%s.fq > %sunmapped_R2%s.fq" % (a, b, a, b)
    subprocess.call([reverse_unmapped], shell=True)
    
### Manipulating the location of unmapped files 
#This code makes a list of unmapped files in preparation for Kraken
os.chdir(PE) #Change the current working directory
pwd = os.getcwd() #Get the current working directory
make_file2 = "ls *R1* > unmapped.filenames" 
subprocess.check_output(make_file2, shell=True) #Making file from unmapped names
#This code manipulates the location of the alignment files to Alignment directory
copy_file1 = "mv *bam %s" % (Alignment_dir)
subprocess.check_output(copy_file1, shell=True) #Moving bam files to Alignment directory
copy_file2 = "mv *unmapped* %s" % (Alignment_dir)
subprocess.check_output(copy_file2, shell=True) #Moving unmapped files to Alignment directory

############################ SPECIES IDENTIFICATION ###########################

### Identifying unmapped reads with Kraken 
#This code puts user in the Paired end directory
os.chdir(Alignment_dir) #Change the current working directory
pwd = os.getcwd() #Get the current working directory
print("Current working directory:", pwd) #Print the current working directory
#This code identifies unmapped reads with Kraken
file = open(unmapped.filenames)
for line in file:
    a = line.split("R1")[0].strip()
    b = line.split("R1")[-1].strip()
    kraken_cmd = "kraken2 --db %s --paired %sR1%s %sR2%s --threads %s --output %s/Kraken_Outputfile --report %s/Kraken_Report --use-name" % (database_path, a, b, a, b, kraken_threads, Kraken_dir, Kraken_dir)
    subprocess.check_output([kraken_cmd], shell=True)

