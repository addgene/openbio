# Research
Scripts and stuff that's useful and public from Addgene Research

## Deconvolution
Detect contaminating sequences in FASTQ files

This tool finds unique sequences within a dictionary of AAV plasmid samples. These unique sequences, or chunks, are then compared to the reads in specified FASTQ files. The program counts the number of times that a unique chunk matches a read within a FASTQ file and calculates the percent of matches found for each AAV sample. This output allows the user to determine whether their sample is contaminated with other samples and how much contaminating material is in the sample.

This script was originally written by Benjie Chen and has been edited to assist in Addgene's Viral Genome Sequencing pipeline.


## Deconvolution v2

This version of the deconvolution script uses a csv file as input rather than hardcoded sequences within the script. The first column of the csv file should contain sequence ids, and the second column should contain the associated sequences. The command to run this script should look like:
./deconvolution_v2.py /path/to/FASTQ /path/to/csv
