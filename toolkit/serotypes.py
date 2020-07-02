import click
import os
import logging

from collections import namedtuple
from datetime import date, datetime
from pandas import DataFrame
from utils import (
    read_fastq_file_with_reverse_complements,
    get_fastq_files, setup_dirs,
    get_params_for_command
)

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
    start = datetime.now()
    params = get_params_for_command('serotypes')
    setup_dirs(params)
    input_files = get_fastq_files(params)
    L.info('Found {} FASTQ files.'.format(len(input_files)))
    processor = Processor(input_files)
    processor.process()
    L.info('\nDone! the following two files were created:\n')
    L.info(processor.full_output_file)
    L.info(processor.summary_output_file)
    delta = datetime.now() - start
    L.info('\nCommand took {} seconds.'.format(delta.total_seconds()))

##################

main_columns = ['File Name', 'Read count', 'Signature Name', 'Occurrences', 'Is top occurrence?',
                'Does filename match detected serotype?']
summary_columns = ['File Name', 'Read count', 'Top Signature', 'Occurrences of Top Signature',
                   'Does filename match detected serotype?', 'Other Signatures Found']

SignatureCount = namedtuple('SignatureCount', ['name', 'count'])


class Processor(object):

    def __init__(self, input_files):
        self.params = get_params_for_command('serotypes')
        self.signatures = self.params.get('signatures')
        self.output_folder = self.params.get('output_folder', '')
        if not self.signatures:
            raise ValueError('Could not find "signatures" parameter in parameters.yml')

        self.input_files = input_files
        self.main_rows = []
        self.summary_rows = []
        self.full_output_file = self._get_output_filename('full')
        self.summary_output_file = self._get_output_filename('summary')

    def process(self):
        for input_file in self.input_files:
            L.info('Processing file: {}'.format(input_file))
            try:
                self._process_one_file(input_file)
            except ValueError as e:
                L.error('\nThe file {} is not a valid FASTQ file. Issue: {}'.format(input_file, str(e)))
                raise ValueError('Exiting now. Please fix the file issue and re-run the command.')

        main_data_frame = DataFrame(self.main_rows, columns=main_columns)
        main_data_frame.to_csv(self.full_output_file, index=False)

        summary_data_frame = DataFrame(self.summary_rows, columns=summary_columns)
        summary_data_frame.to_csv(self.summary_output_file, index=False)


    def _process_one_file(self, input_file):
        filename = os.path.basename(input_file)
        reads = read_fastq_file_with_reverse_complements(input_file)
        # We include the reverse complement in reads, but don't need to include in count
        number_of_reads = len(reads) / 2

        signature_counts = []
        aggregate_signature_counts = []
        for name, signature_list in self.signatures.items():
            total_count = 0
            for index, signature in enumerate(signature_list):
                subname = '{}-{}'.format(name, index + 1)
                count = self._count_signature(reads, signature)
                total_count += count
                signature_counts.append(SignatureCount(name=subname, count=count))
            aggregate_signature_counts.append(SignatureCount(name=name, count=total_count))

        # sort the counts such that the highest count is reported first
        signature_counts.sort(key=lambda x:x.count, reverse=True)
        aggregate_signature_counts.sort(key=lambda x:x.count, reverse=True)

        for index, signature_count in enumerate(signature_counts):
            is_top = 'YES' if index == 0 and signature_count.count else ''
            if is_top == 'YES':
                is_match = 'YES' if (
                        index == 0 and signature_count.count and self._is_match(signature_count.name.split('-')[0],
                                                                                input_file)
                ) else 'NO'
            else:
                is_match = ''

            self.main_rows.append(
                (filename, number_of_reads, signature_count.name, signature_count.count, is_top, is_match)
            )

        found_signatures = []
        for aggregate_signature_count in aggregate_signature_counts:
            if aggregate_signature_count.count:
                found_signatures.append('{} ({})'.format(aggregate_signature_count.name,
                                                         aggregate_signature_count.count))


        other_signatures = ''
        if len(found_signatures) > 1:
            other_signatures = '; '.join(found_signatures[1:])

        if not self._undetermined(aggregate_signature_counts):
            top_signature = aggregate_signature_counts[0]
            is_match = 'YES' if self._is_match(top_signature.name, input_file) else ''
            name = top_signature.name
            count = top_signature.count
        else:
            is_match = ''
            name = 'Not deteremined'
            count = ''

        summary_row = (filename, number_of_reads, name,
                count, is_match, other_signatures)

        self.summary_rows.append(summary_row)

    @staticmethod
    def _undetermined(signature_counts):
        if not signature_counts[0].count:
            return True
        if len(signature_counts) > 1:
            # Tie
            return signature_counts[0].count == signature_counts[1].count
        return False

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

    def _get_output_filename(self, suffix):
        local_name = '{}_serotype_report_{}.csv'.format(date.today().isoformat(), suffix)
        return os.path.join(self.output_folder, local_name)


