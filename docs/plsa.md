[addgene/openbio/docs](https://addgene.github.io/openbio)
# The Pooled Libraries Skew Analysis (plsa) command
The __plsa__ command helps Addgene’s scientist team assess the skew of a given pooled library during pooled library amplification. The command plots a Lorenz curve based on the counts provided from a pooled library sequencing experiment. This plot should be created for the table of counts generated from each library sequenced. Addgene scientists compare these plots for pooled libraries before and after amplification to ensure no further skewing has occurred.

The initial script was authored by Max Juchheim.

## Configuration

## Procedure
1. Make sure you have [downloaded](https://github.com/addgene/openbio/archive/master.zip) and expanded the latest code into your Home folder
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
    python atk.py plsa XXX
    ```
1. Once the command finishes, XXX
