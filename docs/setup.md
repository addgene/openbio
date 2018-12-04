[addgene/openbio/docs](https://addgene.github.io/openbio)
# Setting up your development environment in a Mac

You will use the Terminal application to run a script that installs everything you need. A few notes on the Terminal app:
* The Terminal application can be launched by going to Applications > Utilities and double clicking on Terminal. We suggest adding this application to your Dock.
* In the procedure below, commands that you need to type in the terminal window are in `this font`.
* While on the Terminal, you can use TAB to complete a command after typing a few letters. Hit ENTER to execute the command. Use the Up arrow to get to the previous commands.

## Procedure
1. [Download](https://github.com/addgene/openbio/archive/master.zip) the `addgene/openbio` GitHub repository and move it to your Home folder.
1. Unzip it by double clicking on it. The top-level folder is called `openbio-master`
1. Launch the Terminal application and navigate to the setup folder:
    ```
    cd openbio-master/setup
    ```
1. Type the following command:
    ```
    ./setup-openbio-workstation.sh
    ```
    This does two things in parallel:
    * Launches installation of XCode, a set of Mac tools for developers - you’ll see this in a separate window.
    * Downloads and configures everything else you need - you’ll see text scrolling in the Terminal indicating what’s going on. At some point you may be asked for your password - this is the password you use for your local Mac account.
1. When the script is done, you’ll see a success message. Close your Terminal window and open a new one so everything is loaded freshly.
1. Test that you can run our Python scripts:
    * Activate your Python environment - after this, the prompt in the Terminal changes to indicate your environment (openbio) is active:
    ```
    workon openbio
    ```
    
    * Navigate to the toolkit folder: 
    ```
    cd openbio-master/toolkit
    ```
    * Execute the `atk.py` command with the `--help` option to see a help message:
    ```
    python atk.py --help
    ```
1. In order to use the toolkit, you will need to enter your own parameters in a python file (`parameters.py`), which can be done with any text editor. Because XCode was installed, your Mac will want to open .py files with XCode, which we don’t want. The first time you need to edit a Python file, right click on the file and choose “Open with…”. To permanently associate a text editor with .py files: select the file and hit <command-i>, go to the “Open with” section and select your favorite text editor - if you don’t have a favorite, just use TextEdit.app.
