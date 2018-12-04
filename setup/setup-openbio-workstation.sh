#! /bin/bash

set -e

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo -e "\n** Creating / updating global prerequisites...\n"
$script_dir/1-setup-usr-local.sh $1

echo -e "\n** Creating / updating virtual environment...\n"
$script_dir/2-setup-openbio-virtualenv.sh $1

echo -e "\n---\n"
echo -e "Your environment was correctly set up!\n"
echo -e "Now close the terminal window and open a fresh one (so everything is loaded cleanly)\n"

