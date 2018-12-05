[Addgene](https://addgene.org) is a nonprofit organization whose mission is to accelerate research and discovery by improving access to useful research materials and information. Since the company was founded in 2004, scientists have shared more than 70,000 unique published reagent samples via Addgene’s repository. We have also fulfilled requests for more than 1 million plasmid, viral vector, and other material types to scientists in 100 countries.

Addgene conducts a rigorous quality control process for all [plasmids](https://blog.addgene.org/plasmids-101-an-inside-look-at-ngs-plasmid-quality-control) and 
[viral vectors](https://blog.addgene.org/aav-vector-quality-control-going-the-extra-mile) we distribute. Part of this quality control involves sequencing the materials using Next Generation Sequencing (NGS). As our NGS data volume grows, we have had to find ways to automate many parts of our analysis. We do so through internally developed bioinformatics tools. Our language of choice is Python.

The aim of this code repository is to make these tools available to the broader scientific community.
Stay tuned as *the Addgene Toolkit* grows!

# Setting up your Python environment
First of all, [download](https://github.com/addgene/openbio/archive/master.zip) and expand the `addgene/openbio` GitHub repository. 
If you already have a Python environment of your liking, you simply need to navigate to the `openbio-master` root, issue the command:
```
pip install -r requirements.txt
```
and move on.

If you don't already have a Python environment, and you have a Mac, you are welcome to use the shell script that we use internally. You need to have admin privileges on your Mac. 
Follow [these instructions](https://addgene.github.io/openbio/setup).

# The Addgene Toolkit
## Basics
The entry point to the Addgene Toolkit is the command `atk`.
Open up a Terminal and navigate to the toolkit directory:
```
cd openbio-master/toolkit
```
If you used our shell script to set up, activate the virtual environment by issuing the following command:
```
workon openbio
```
Run `atk` as follows to learn more:
```
python atk.py --help
```
This lists the tools that are available as commands and their parameters. The general pattern to invoke a command is:
```
python atk.py [command]
```

## Commands
1. __[Serotypes](https://addgene.github.io/openbio/serotypes)__ - detect and report specific signatures in NGS data, useful to differentiate viral serotypes.
1. __Recombination__ (coming soon) - detect recombination in NGS data, used for Cre-Lox quality control.

