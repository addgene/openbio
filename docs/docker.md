[addgene/openbio/docs](https://addgene.github.io/openbio)
# Using the Docker container

## Prerequisites

### Docker runtime
To build and run Docker containers, you need to have Docker installed in your computer. If you are on a Mac or 
Windows, download and install [Docker Desktop](https://www.docker.com/products/docker-desktop).

### Terminal application
You will need to use a Terminal application to build and run the Docker container.

If you're on a Mac, The Terminal app can be launched by going to Applications > Utilities and double clicking on 
Terminal.
We suggest adding this application to your Dock. The commands that you need to type in the Terminal window are in `this 
font`.
* Cheat sheet for the Mac Terminal app:
    * Use TAB to complete a command or file name after typing a few letters.
    * Hit ENTER to execute the command.
    * Use Ctrl-C to interrupt a command.
    * Use the Up arrow to get to the previously typed commands.

## Building the Docker container
1. [Download](https://github.com/addgene/openbio/archive/main.zip) 
   the `addgene/openbio` repository and move it to your Home folder.
1. Unzip it by double-clicking on it. The top-level folder is called `openbio-main`
1. Launch the Terminal application and navigate to the top-level folder:
    ```
    cd openbio-main
    ```
1. Type the following command:
    ```
    docker compose build openbio
    ```

    This will take a few minutes. 

## Running the Docker container

1. Issue the following command:

    ```
    docker compose run openbio
    ```
    This will take you to a shell prompt in the container. From this shell prompt you can use the Addgene Toolkit as 
   explained [here](https://addgene.github.io/openbio) or in the individual command instructions.

1. Exit the container by typing the command `exit` at the shell prompt.
2. In order to use the toolkit, you will need to enter your own parameters in `parameters.yaml`, which 
   can be done with any text editor. 
   2. Inside the container, your Home directory is available as `/Home`. Please use this root path when entering 
      folder information in the parameter file.
