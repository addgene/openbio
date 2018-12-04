from __future__ import division

import click
import os
import logging

from collections import namedtuple
from datetime import date
from pandas import DataFrame
from parameters import Serotypes_Parameters as params
from utils import read_fastq_file, get_fastq_files, setup_dirs

L = logging.getLogger(__name__)

#### Command ####
@click.command()
def serotypes():
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

##################

main_columns = ['File Name', 'Read count', 'Signature Name', 'Occurrences', 'Is top occurrence?', 'Does it match?']
summary_columns = ['File Name', 'Read count', 'Top Signature', 'Occurrences of Top Signature', 'Does it mmtch?',
                   'Other Signatures Found']

SignatureCount = namedtuple('SignatureCount', ['name', 'count'])


class Processor(object):

    def __init__(self, input_files):
        self.input_files = input_files
        self.main_rows = []
        self.summary_rows = []
        self.full_output_file = self._get_output_filename('full')
        self.summary_output_file = self._get_output_filename('summary')

    def process(self):
        for input_file in self.input_files:
            L.info('Processing file: {}'.format(input_file))
            self._process_one_file(input_file)

        main_data_frame = DataFrame(self.main_rows, columns=main_columns)
        main_data_frame.to_csv(self.full_output_file, index=False)

        summary_data_frame = DataFrame(self.summary_rows, columns=summary_columns)
        summary_data_frame.to_csv(self.summary_output_file, index=False)


    def _process_one_file(self, input_file):
        filename = os.path.basename(input_file)
        reads = read_fastq_file(input_file)
        # We include the reverse complement in reads, but don't need to include in count
        number_of_reads = len(reads) / 2

        signature_counts = []
        for name, signature in params.signatures.items():
            count = self._count_signature(reads, signature)
            signature_counts.append(SignatureCount(name=name, count=count))

        # sort the counts such that the highest count is reported first
        signature_counts.sort(key=lambda x:x.count, reverse=True)

        found_signatures = []
        for index, signature_count in enumerate(signature_counts):
            if signature_count.count:
                found_signatures.append('{} ({})'.format(signature_count.name, signature_count.count))
            is_top = 'YES' if index == 0 and signature_count.count else ''
            if is_top == 'YES':
                is_match = 'YES' if (
                        index == 0 and signature_count.count and self._is_match(signature_count.name, input_file)
                ) else 'NO'
            else:
                is_match = ''

            self.main_rows.append(
                (filename, number_of_reads, signature_count.name, signature_count.count, is_top, is_match)
            )

        top_signature = signature_counts[0]
        other_signatures = ''
        if len(found_signatures) > 1:
            other_signatures = '; '.join(found_signatures[1:])

        is_match = 'YES' if self._is_match(top_signature.name, input_file) else ''


        summary_row = (filename, number_of_reads, top_signature.name,
                top_signature.count, is_match, other_signatures)

        self.summary_rows.append(summary_row)

    @staticmethod
    def _count_signature(reads, signature):
        count = 0
        for read in reads:
            count += read.lower().count(signature.lower())
        return count

    @staticmethod
    def _is_match(signature_name, input_file):
        if not signature_name or not input_file:
            return False
        if '/' in input_file:
            input_file = os.path.basename(input_file)
        return signature_name.split('-')[0].lower() in input_file.lower()

    @staticmethod
    def _get_output_filename(suffix):
        local_name = '{}_serotype_report_{}.csv'.format(date.today().isoformat(), suffix)
        return os.path.join(params.output_folder, local_name)


