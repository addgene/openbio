#! /bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ -e ~/.bashrc ]; then
    config_file="$HOME/.bashrc"
elif [ -e ~/.bash_profile ]; then
    config_file="$HOME/.bash_profile"
else
    config_file="$HOME/.profile"
fi
ve_rc_file="virtualenvrc"

grep_line='[ -f ~/.virtualenvrc ] && source ~/.virtualenvrc'
grep -q -F "$grep_line" $config_file || echo $grep_line >> $config_file
cat $script_dir/$ve_rc_file > $HOME/.$ve_rc_file

echo -e "A virtualenv helper script placed at $HOME/.${ve_rc_file} has been ref'd in your ${config_file}.\n"
