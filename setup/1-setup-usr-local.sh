#! /bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
sudo_message="\nIf prompted for 'password' enter your mac/workstation password\n"

if [ "$1" == "nuke" ]; then
    echo -e "\n\n** You have selected the nuke option.  This will delete absolutely everything in /usr/local.  Hit ctrl-c now if this was unintended **\n\n"
    read -rp "Hit enter to continue.  Ctrl-c to exit"
    echo -e "$sudo_message"
    sudo rm -rf /usr/local/*
fi

(xcode-select --install && echo -e "\n\n" && read -rp "Installing xcode command line tools.  Hit enter to continue once it is done") || true

test -d /usr/local/Homebrew || /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

chown_cmd="chown -R $(whoami):admin /usr/local/*"

$chown_cmd || (echo "At least one file/dir may not be owned by you.  Enter your workstation password if prompted" && sudo $chown_cmd)

chmod 755 /usr/local/share

if [ ! -d /usr/local/man ]; then
    echo -e "$sudo_message"
    sudo mkdir /usr/local/man
    sudo chown "$(whoami):admin" /usr/local/man
fi

brew update
brew upgrade

while read -r homebrew; do
    brew install "${homebrew}"
done < "${script_dir}/packages/homebrew.txt"

brew cleanup

# ref /usr/local
/usr/local/bin/pip3 install --upgrade pip
while read -r pips; do
    /usr/local/bin/pip3 install "${pips}"
done < "${script_dir}/packages/pip.txt"
