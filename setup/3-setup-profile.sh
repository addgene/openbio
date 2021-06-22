#! /bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -e ~/.bashrc ]; then
    config_file="$HOME/.bashrc"
else
    config_file="$HOME/.bash_profile"
fi
ve_rc_file="virtualenvrc"

grep_line='[ -f ~/.virtualenvrc ] && source ~/.virtualenvrc'
grep -q -F "$grep_line" "$config_file" || echo "$grep_line" >> "$config_file"
cat "$script_dir/config/$ve_rc_file" > "$HOME/.$ve_rc_file"

echo -e "1) A virtualenv helper script placed at $HOME/.${ve_rc_file} has been ref'd in your ${config_file}.\n"
echo -e "The point of this script is to automatically activate virtualenvs when you 'cd' to repos with the same name\n"
echo -e "If you use a shell other than bash, you will need to add/ref the contents of $script_dir/config/$ve_rc_file to your env by your preferred means\n"
