import argparse
import os

import click
import logging
import sys
from collections import namedtuple
from datetime import date

from parameters import Serotype_Report_Parameters as params
from utils import read_fastq_file, get_fastq_files, setup_dirs, log_to_stdout

L = logging.getLogger(__name__)

row_template = '{},{},{},{},{},{}\n'
summary_row_template = '{},{},{},{},{},{}\n'
SignatureCount = namedtuple('SignatureCount', ['name', 'count'])

class Processor(object):

    def __init__(self, input_files):
        self.input_files = input_files
        self.rows = []
        self.rows.append(row_template.format(
                'File Name', 'Read count', 'Signature Name',
                'Occurrences', 'Is top occurrence?', 'Does it match?')
        )
        self.summary_rows = []
        self.summary_rows.append(summary_row_template.format(
                'File Name', 'Read count', 'Top Signature',
                'Occurrences of Top Signature', 'Other Signatures Found', 'Does it Match?')
        )
        self.full_output_file = _get_output_filename('full')
        self.summary_output_file = _get_output_filename('summary')

    def process(self):
        for input_file in self.input_files:
            L.info('Processing file: {}...'.format(input_file))
            self._process_one_file(input_file)

        _write_to_file(self.rows, self.full_output_file)
        _write_to_file(self.summary_rows, self.summary_output_file)


    def _process_one_file(self, input_file):
        filename = os.path.basename(input_file)
        reads = read_fastq_file(input_file)
        number_of_reads = len(reads)

        signature_counts = []
        for name, signature in params.signatures.iteritems():
            count = self._count_signature(reads, signature)
            signature_counts.append(SignatureCount(name=name, count=count))

        signature_counts.sort(key=lambda x:x.count, reverse=True)

        found_signatures = []
        for index, signature_count in enumerate(signature_counts):
            if signature_count.count:
                found_signatures.append('{} ({})'.format(signature_count.name, signature_count.count))
            is_top = 'YES' if index == 0 and signature_count.count else ''
            is_match = 'YES' if (
                    index == 0 and signature_count.count and signature_count.name.lower() in input_file.lower()
            ) else ''

            self.rows.append(row_template.format(
                    filename, number_of_reads, signature_count.name, signature_count.count, is_top, is_match)
            )

        top_signature = signature_counts[0]
        other_signatures = ''
        if len(found_signatures) > 1:
            other_signatures = ', '.join(found_signatures[1:])

        is_match = 'YES' if top_signature.name.lower() in input_file.lower() else ''
        summary_row = summary_row_template.format(
                filename, number_of_reads, top_signature.name,
                top_signature.count, other_signatures, is_match)

        self.summary_rows.append(summary_row)


    def _count_signature(self, reads, signature):
        count = 0
        for read in reads:
            count += read.lower().count(signature.lower())
        return count


def _write_to_file(rows, output_file):
    with open(output_file, 'w') as file:
        for row in rows:
            file.write(row)

def _get_output_filename(suffix):
    local_name = '{}_serotype_report_{}.csv'.format(date.today().isoformat(), suffix)
    return os.path.join(params.output_folder, local_name)

##### Command #####
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
def serotype_report_command():
    """
    Reporty specific sequences in VGS data. This command reads FASTQ files from the input folder
    specified in Serotype_Report_Parameters, extracts the reads containing the sequence and counts the occurrences.
    The Research Team is using this tool to detect specific sequences in the capsid genes,
    thereby differentiating between various serotypes.
    """

    setup_dirs(params)
    input_files = get_fastq_files(params)
    L.info('Found {} FASTQ files.'.format(len(input_files)))
    processor = Processor(input_files)
    processor.process()
    L.info('Done! the following two files were created:\n')
    L.info(processor.full_output_file)
    L.info(processor.summary_output_file)

if __name__ == "__main__":
    log_to_stdout(logging.INFO)
    try:
        serotype_report_command()
    except ValueError as e:
        L.error('\nError: ' + e.message)
        sys.exit(1)
