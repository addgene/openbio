[addgene/openbio/docs](https://addgene.github.io/openbio)
# The Serotypes Command
The __serotypes__ command helps Addgene's Research team detect specific sequences in Viral Next Generation Sequencing (VGS) data. The command reads all FASTQ files from an input folder, extracts the reads and counts the occurrences of each sequence. The Research team uses this command to detect specific sequences in the capsid genes, thereby differentiating between various serotypes.

The command’s parameters can be modified by editing the file `parameters.py` using your favorite text editor) To change the parameter values, locate the block named `Serotypes_Parameters` and follow the examples in the file, paying special attention to the use of double quotes for all text. The parameters are:
* __input_folder__: the folder where the FASTQ files are. Enter the full path or a path relative to the toolkit folder.
* __output_folder__: the folder for the output CSV files. Enter the full path or a path relative to the toolkit folder.
* __signatures__: name and sequence of the signatures that the command will look for. The Research team has identified a number of signatures that correlate to specific serotypes, and this parameter is pre-populated with them.  Add as many as you want, following the example syntax, and separated them by commas. 

## Procedure
1. Make sure you have downloaded and expanded the latest code into your Home folder
1. Adjust the parameters for the command as described above, by editing the file `parameters.py`.
1. Open a Terminal window and activate your Python environment:
    ```
    workon openbio
    ```
1. Navigate to the toolkit folder:
    ```
    cd openbio/toolkit
    ```
1. Issue the following command:
    ```
    python atk.py serotypes
    ```
1. Once the command finishes, you will find the output CSV files in the folder you selected. Two files will be generated: a full report and a summary. The file names will contain the date when the report was generated.
