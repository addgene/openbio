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

To start the container, issue the following command:

 ```
 docker compose run openbio
 ```

This will take you to a shell prompt in the container. From this shell prompt you can use the Addgene Toolkit as 
explained [here](https://addgene.github.io/openbio) or in the individual command instructions.

Exit the container by typing the command `exit` at the shell prompt.

## Working directory

Inside the container, the path`/workdir` should be used to specify input and output file paths. Please use this root 
path in `parameters.yml`. 

The `/workdir` path is mapped by default to the`openbio-main/workdir` folder in your host 
computer. This means that any files you place in the `openbio-main/workdir` folder will be available inside the container in the `/workdir` 
folder. Similarly, any files you create in the `/workdir` folder inside the container will be available in the 
`openbio-main/workdir` folder in your host computer. This is a convenient way to share your work files between 
the host computer and the container.


If you want to use a different folder in your host computer to share files with the container, start the container 
with the following command:
 ```
 OPENBIO_WORKDIR=[YOUR_DIRECTORY] docker compose run openbio
 ```
For example, if you want to use the host computer folder `/Users/Harry/fastq-data`, use 
the following command:
 ```
 OPENBIO_WORKDIR=/Users/Harry/fastq-data docker compose run openbio
 ```
