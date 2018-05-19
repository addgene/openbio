#! /bin/bash

ve_name="research"

export PIP_RESPECT_VIRTUALENV=true
export PIP_VIRTUALENV_BASE=$WORKON_HOME
export VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python
export VIRTUALENVWRAPPER_VIRTUALENV_ARGS='--no-site-packages'
export WORKON_HOME=$HOME/.virtualenvs
source virtualenvwrapper.sh

deactivate||echo "Not in a virtualenv"

if [ "$1" == "nuke" ]; then
    rmvirtualenv $ve_name
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $script_dir

git_dir=$(git rev-parse --show-toplevel)
pushd $git_dir

workon $ve_name

if [ $? -ne 0 ]; then
    mkvirtualenv $ve_name -p python2.7
fi

workon $ve_name

pip install pip==9.0.3

pip install -r requirements.txt
