#! /bin/bash

set -e

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo -e "\n** Creating / updating global prerequisites...\n"
$script_dir/1-setup-usr-local.sh $1

echo -e "\n** Creating / updating Bioinformatics virtual environment...\n"
$script_dir/2-setup-prodops-virtualenv.sh $1

echo -e "\n---\n"
echo -e "All pre-requisites are setup OK!\n"
echo -e "Now close all of your terminal windows and open a fresh one (so everything is loaded cleanly)\n"
echo -e "After you have done that, navigate back to the research repo directory and type: workon bioinformatics
