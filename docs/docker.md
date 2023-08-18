[addgene/openbio/docs](https://addgene.github.io/openbio)
# Using the Docker container

## Prerequisites

### Docker runtime
To build and run Docker containers, you need to have Docker installed in your computer. If you are on a Mac or 
Windows, download and install [Docker Desktop](https://www.docker.com/products/docker-desktop).

## Building the Docker container
1. [Download](https://github.com/addgene/openbio/archive/main.zip) 
   the `addgene/openbio` repository and move it to your Home folder. Unzip it by double-clicking on it. The top-level folder is called `openbio-main`
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
   Inside the container, your Home directory is available as `/Home`. Please use this root path when entering folder information in the parameter file.
