#! /bin/bash

ve_name="openbio"

export PIP_RESPECT_VIRTUALENV=true
export PIP_VIRTUALENV_BASE=$WORKON_HOME
export VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python3
export WORKON_HOME=$HOME/.virtualenvs
unset VIRTUALENVWRAPPER_VIRTUALENV_ARGS
source virtualenvwrapper.sh

deactivate || echo "Not in a virtualenv"

if [ "$1" == "nuke" ]; then
    rmvirtualenv $ve_name
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd "$script_dir" || exit 1

git_dir=$(git rev-parse --show-toplevel)
pushd "$git_dir" || exit 1

if [ ! -e "$WORKON_HOME/$ve_name/bin/activate" ]; then
    mkvirtualenv $ve_name -p python3
else
    # always ensure latest python is used in the virtualenv
    gfind "${WORKON_HOME}/$ve_name/" -type l -xtype l -delete
    virtualenv -p python3 "${WORKON_HOME}/$ve_name"
fi

workon $ve_name

pip3 install --upgrade pip

pushd "$script_dir" || exit 1

pip3 install -r ../requirements.txt
