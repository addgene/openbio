#! /bin/bash

ve_name="openbio"

export PIP_RESPECT_VIRTUALENV=true
export PIP_VIRTUALENV_BASE=$WORKON_HOME
export VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python3
export VIRTUALENVWRAPPER_VIRTUALENV_ARGS='--no-site-packages'
export WORKON_HOME=$HOME/.virtualenvs
source virtualenvwrapper.sh

deactivate||echo "Not in a virtualenv"

if [ "$1" == "nuke" ]; then
    rmvirtualenv $ve_name
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd "$script_dir" || exit 1

git_dir=$(git rev-parse --show-toplevel)
pushd "$git_dir" || exit 1

if [ ! -e "$WORKON_HOME/$ve_name/bin/activate" ]; then
    mkvirtualenv $ve_name
fi

workon $ve_name

pip3 install -r requirements.txt

grep VIRTUALENVWRAPPER_PYTHON $HOME/.bash_profile || echo "export VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python3" >> $HOME/.bash_profile
grep "source virtualenvwrapper.sh" $HOME/.bash_profile || echo "source virtualenvwrapper.sh" >> $HOME/.bash_profile
