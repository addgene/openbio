import logging
import os
import sys
from typing import List, Dict

import yaml
from Bio import SeqIO
from setuptools import glob

L = logging.getLogger(__name__)
PARAM_FILE_NAME = 'parameters.yml'

def get_params_for_command(command_name: str, param_file: str=PARAM_FILE_NAME) -> Dict:
    L.info('Using configuration file: {}'.format(param_file))
    with open(param_file, 'r') as ymlfile:
        params = yaml.safe_load(ymlfile)

    return params.get(command_name, {}) if params else {}

def read_fastq_file_with_reverse_complements(file_name:str) -> List[str]:
    reads = []
    with open(file_name, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            reads.append(str(record.seq))
            reads.append(str(record.seq.reverse_complement()))
    return reads


def read_fastq_file(file_name:str, to_upper:bool=True) -> List[str]:
    with open(file_name, "r") as handle:
        if to_upper:
            reads = [str.upper(str(record.seq)) for record in SeqIO.parse(handle, "fastq")]
        else:
            reads = [str(record.seq) for record in SeqIO.parse(handle, "fastq")]
    return reads


def setup_dirs(params: Dict):
    input = params.get('input_folder')
    output = params.get('input_folder')

    if not os.path.isdir(input):
        raise ValueError('The folder ' + input + ' does not exist.')

    if not os.path.isdir(output):
        L.info('\nThe folder ' + output + ' does not exist. Creating it...\n')
        os.mkdir(output)


def get_fastq_files(params: Dict) -> List[str]:
    # grab all files ending in .fastq
    input_files = [input_file for input_file in glob.glob("{}/*.fastq".format(params.get('input_folder')))]
    input_files.extend([input_file for input_file in glob.glob("{}/*.fq".format(params.get('input_folder')))])
    if not input_files:
        raise ValueError('No FASTQ files in folder: ' + params.get('input_folder'))

    return input_files


def log_to_stdout(level):
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    logging.getLogger().addHandler(handler)
    logging.getLogger().setLevel(level)
