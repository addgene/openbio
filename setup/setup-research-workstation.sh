#! /bin/bash

set -e

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo -e "\n** Creating / updating global prerequisites...\n"
$script_dir/1-setup-usr-local.sh $1

echo -e "\n** Creating / updating Bioinformatics virtual environment...\n"
$script_dir/2-setup-research-virtualenv.sh $1

$script_dir/3-setup-profile.sh $1

echo -e "\n---\n"
echo -e "Your environment was correctly set up!\n"
echo -e "Now close the terminal window and open a fresh one (so everything is loaded cleanly)\n"
echo -e "After you have done that, navigate back to the research directory to start executing scripts.\n"

