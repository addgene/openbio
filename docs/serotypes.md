[addgene/openbio/docs](https://addgene.github.io/openbio)
# The Serotypes Command
The __serotypes__ command helps Addgene’s Research team determine the serotype of a viral vector prep using Next Generation Sequencing (NGS) of the prep. Addgene has identified unique sequences for commonly used capsids. The command reads all FASTQ files from an input folder, extracts the reads and counts the occurrences of each sequence. A small amount of RepCap plasmid is packaged within the AAV, so the capsid used for production will be represented in the FASTQ data. For example, an AAV2 vector should return counts for the AAV2 sequence, but not for the AAV5 or other sequences.

### Notes: 
* Occasionally, you will have spurious matches for other serotypes, but one sequence should be the clear majority. 
* Only a small amount of RepCap plasmid is packaged, so if you don’t have any matches it may be that you need a higher number of NGS reads. 
* If you are using a different capsid sequence (for instance, you may have the same amino acid sequence, but have a different DNA sequence than the one Addgene uses in our RepCap plasmids), you will need to adjust the capsid sequences that the program is searching for such that your capsid sequence matches the one in the program. See the Configuration section below for how to add your own signatures.
* If you include the name of the serotype you expect in the FASTQ file name, the command will report if the top match corresponds to this expectation.

## Configuration
The command’s parameters can be modified by editing the file `parameters.py` (using your favorite text editor) To change the parameter values, locate the block named `Serotypes_Parameters` and follow the examples in the file, paying special attention to the use of double quotes for all text. The parameters are:
* __input_folder__: the folder where the FASTQ files are. Enter the full path or a path relative to the toolkit folder.
* __output_folder__: the folder for the output CSV files. Enter the full path or a path relative to the toolkit folder.
* __signatures__: name and sequences of the signatures that the command will look for. This parameter is pre-populated with the signatures that Addgene has identified for commonly used capsids.  Note that it is possible to specify more than one sequence to match for a given capsid. You may add other sequences if you need to, just follow the example syntax. 

## Procedure
1. Make sure you have [downloaded](https://github.com/addgene/openbio/archive/master.zip) and expanded the latest code into your Home folder
1. Adjust the parameters for the script by editing the file `parameters.py` as described above.
1. Open a Terminal window and activate your Python environment:
    ```
    workon openbio
    ```
1. Navigate to the toolkit folder:
    ```
    cd openbio-master/toolkit
    ```
1. Issue the following command:
    ```
    python atk.py serotypes
    ```
1. Once the command finishes, you will find the output CSV files in the folder you selected. Two files will be generated: a full report and a summary. The file names will contain the date when the report was generated.
