import logging
import os
import sys

from Bio import SeqIO
from setuptools import glob

L = logging.getLogger(__name__)


def read_fastq_file(file_name):
    # type: (str) -> List[str]
    reads = []
    with open(file_name, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            reads.append(str(record.seq))
            reads.append(str(record.seq.reverse_complement()))
    return reads


def setup_dirs(params):
    input = params.input_folder
    output = params.output_folder

    if not os.path.isdir(input):
        raise ValueError('The folder ' + input + ' does not exist.')

    if not os.path.isdir(output):
        L.info('\nThe folder ' + output + ' does not exist. Creating it...\n')
        os.mkdir(output)


def get_fastq_files(params):
    # type: (object) -> List[str]

    # grab all files ending in .fastq
    input_files = [input_file for input_file in glob.glob("{}/*.fastq".format(params.input_folder))]
    if not input_files:
        raise ValueError('No FASTQ files in folder: ' + params.input_folder)

    return input_files


def log_to_stdout(level):
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    logging.getLogger().addHandler(handler)
    logging.getLogger().setLevel(level)