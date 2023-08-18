[Addgene](https://addgene.org) is a nonprofit organization whose mission is to accelerate research and discovery by improving access to useful research materials and information. Since the company was founded in 2004, scientists have shared more than 100,000 unique published reagent samples via Addgene’s repository. We have also fulfilled requests for more than 1.5 million plasmid, viral vector, and other material types to scientists in 100+ countries.

Addgene conducts a rigorous quality control process for all [plasmids](https://blog.addgene.org/plasmids-101-an-inside-look-at-ngs-plasmid-quality-control) and 
[viral vectors](https://blog.addgene.org/aav-vector-quality-control-going-the-extra-mile) we distribute. Part of this quality control involves sequencing the materials using Next Generation Sequencing (NGS). As our NGS data volume grows, we have had to find ways to automate many parts of our analysis. We do so through internally developed bioinformatics tools. Our language of choice is Python 3.

The aim of this code repository is to make these tools available to the broader scientific community.
Stay tuned as *the Addgene Toolkit* grows!

# Setting up your environment
First of all [download](https://github.com/addgene/openbio/archive/main.zip) the 
`addgene/openbio` zipfile, move it to your Home folder, and expand it. 
If you already have a Python 3 environment of your liking, you simply need to navigate to the `openbio-main` root,
issue the command:
```
pip install -r requirements.txt
```
and move on.

If you don't already have a Python environment, you may use our Docker container, which includes everything you need to use the toolkit. 
Please follow [these instructions](https://addgene.github.io/openbio/docker) to learn how to build and use it.
# The Addgene Toolkit
## Basics
The entry point to the Addgene Toolkit is the command `atk`.
In a terminal window, navigate to the toolkit folder (if you're using our Docker container, run the
container first):
```
cd openbio-main/toolkit
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
1. __[Recombination](https://addgene.github.io/openbio/recombination)__ - detect recombination in NGS data, used for Cre-Lox quality control.
