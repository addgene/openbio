# The Addgene Open Tools Repository
[Intro on what this repository is]

## Prerequisites
A Python environment

## Installation
Clone this repository and run the following command:
`pip install -r requirements.txt`

Invoke the Addgene Toolkit by issuing the following command in a terminal:
`python atk.py -- help`

## Serotypes
Command to detect specific sequences in VGS data. It reads all FASTQ files from an input folder, extracts the reads and counts the occurrences of each sequence. The Research Team is using this software to detect specific sequences in the capsid genes, thereby differentiating between various serotypes.
The script’s parameters can be modified by editing the file parameters.py. To change the parameter values, locate the block named Serotype_Report_Parameters and follow the examples in the file, paying special attention to the use of double quotes for all text. The parameters are:
input_folder - folder where the FASTQ files are. Full path or relative to the tools folder.
output_folder - folder for the output CSV files. Full path or relative to the tools folder.
signatures - name and sequence of the signatures to look for. Add as many as you want, following the example syntax and separated by commas.
