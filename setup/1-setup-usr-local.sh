#! /bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
sudo_message="\nIf prompted for 'password' enter your mac/workstation password\n"

if [ "$1" == "nuke" ]; then
    echo -e "\n\n** You have selected the nuke option.  This will delete absolutely everything in /usr/local.  Hit ctrl-c now if this was unintended **\n\n"
    read -p "Hit enter to continue.  Ctrl-c to exit"
    echo -e $sudo_message
    sudo rm -rf /usr/local/*
fi

xcode-select --install || true

test -d /usr/local/Homebrew || /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
echo -e $sudo_message
sudo chown -R $(whoami):admin $(brew --prefix)/Homebrew

if [ -d /usr/local/Cellar ]; then
    sudo chown -R $(whoami):admin $(brew --prefix)/Cellar
fi

if [ -d /usr/local/Caskroom ]; then
    sudo chown -R $(whoami):admin $(brew --prefix)/Caskroom
fi

if [ -d /usr/local/var/homebrew ]; then
    sudo chown -R $(whoami):admin $(brew --prefix)/Caskroom
fi

sudo chmod 755 /usr/local/share

if [ ! -d /usr/local/man ]; then
    echo -e $sudo_message
    sudo mkdir /usr/local/man
    sudo chown $(whoami):admin /usr/local/man
fi

brew update
brew upgrade

while read homebrew; do
    brew install ${homebrew}
done < ${script_dir}/packages/homebrew.txt

brew cleanup

# ref /usr/local
/usr/local/bin/pip install --upgrade pip
while read pips; do
    /usr/local/bin/pip install ${pips}
done < ${script_dir}/packages/pip.txt
