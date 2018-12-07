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

echo -e "\nUpdating permissions, this may take a while..."
echo -e $sudo_message
if [ ! -d /usr/local/Frameworks ]; then
    sudo mkdir /usr/local/Frameworks
fi
if [ ! -d /usr/local/man ]; then
    sudo mkdir /usr/local/man
fi

sudo chown -R $(whoami):admin $(brew --prefix)/*
sudo chmod 755 /usr/local/share

brew update
brew upgrade

while read homebrew; do
    brew install ${homebrew}
done < ${script_dir}/packages/homebrew.txt

brew cleanup

ln -s /usr/local/Cellar/python/3.7.1/bin/python3 /usr/local/bin/python
ln -s /usr/local/Cellar/python/3.7.1/bin/pip3 /usr/local/bin/pip

/usr/local/bin/pip install --upgrade pip
while read pips; do
    /usr/local/bin/pip install ${pips}
done < ${script_dir}/packages/pip.txt
