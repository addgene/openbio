# The Recombination Command
The __recombination__ command helps Addgene's Research team detect recombination in VGS data. The command reads all FASTQ files from an input folder, extracts triples of the form (sequence before, seed sequence, sequence after) and counts the occurrences of each triple. In addition, it flags the likely correct sequences and potential recombinations and computes a recombination percentage based on these assumptions. The output is written to CSV files, one per FASTQ file and per seed sequence.

The command’s parameters can be modified by editing the file `parameters.py` in your favorite text editor.

To change the parameter values, locate the block named `Recombination_Parameters` and follow the examples in the file, paying special attention to the use of double quotes for all text. The parameters are:
* __input_folder__: the folder where the FASTQ files are. Enter the full path or a path relative to the toolkit folder.
* __output_folder__: the folder for the output CSV files. Enter the full path or a path relative to the toolkit folder.
* __HEAD__: the number of bases before the seed sequence
* __TAIL__: the number of bases after the seed sequence
* __seed_sequences__: the name and sequence of the seeds to look for. Add as many as you want, following the example syntax and separated by commas.

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
    python atk.py recombination
    ```
1. When the command finishes, you will find the output CSV files in the folder you selected. The file names reflect the parameters used. For example, the output file:

    `A11984_sW0158_H09_R_001-60-lox2272-60.csv`

    Is the output corresponding to:

    * FASTQ file: `A11984_sW0158_H09_R_001.fastq`
    * HEAD = 60, Seed sequence = lox2272, TAIL = 60
