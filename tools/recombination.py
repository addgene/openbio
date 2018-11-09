from __future__ import division
from collections import Counter

import argparse
import click
import logging
import os
import sys
from parameters import Recombination_Parameters as params
from utils import get_fastq_files, read_fastq_file, setup_dirs,log_to_stdout

L = logging.getLogger(__name__)


class GlobalData(object):
    def __init__(self, seed_sequence_name):
        self.seed_sequence_name = seed_sequence_name
        self.total_reads = 0
        self.seed_occurrences = 0
        self.full_sequence_occurrences = 0
        self.recombination_occurrences = 0
        self.recombination_percentage = 0


class Pattern(object):
    def __init__(self, head_sequence, seed_sequence, tail_sequence):
        self.head_sequence = head_sequence
        self.seed_sequence = seed_sequence
        self.tail_sequence = tail_sequence
        self.full_pattern = head_sequence + seed_sequence + tail_sequence

    def __hash__(self):
        return hash((self.full_pattern))
    def __eq__(self, other):
        return (self.full_pattern) == (other.full_pattern)

    def __ne__(self, other):
        return not(self == other)


class Processor(object):
    def __init__(self, input_file, seed_sequence_name, seed_sequence):
        self.input_file = input_file
        self.seed_sequence = seed_sequence
        self.global_data = GlobalData( seed_sequence_name)

    def process(self):
        """
        Processes FASTQ input file and returns a list of CSV rows
        """
        reads = read_fastq_file(self.input_file)

        self.global_data.total_reads = len(reads) / 2
        self.global_data.seed_occurrences = sum(1 for read in reads if self.seed_sequence in read)
        extracted_patterns = [self._extract_pattern(read) for read in reads]

        patterns = Counter(pattern for pattern in extracted_patterns if pattern)

        self.global_data.full_sequence_occurrences = sum(patterns.values())
        return self._patterns_to_rows(patterns)

    def _patterns_to_rows(self, patterns):
        rows = []
        sorted_patterns = patterns.most_common()
        # The two most frequent patterns are likely to be the correct sequences
        reference_patterns = [pattern_tuple[0] for pattern_tuple in sorted_patterns[:2]]

        reference_head_0 = reference_patterns[0].head_sequence
        reference_head_1 = reference_patterns[1].head_sequence
        reference_tail_0 = reference_patterns[0].tail_sequence
        reference_tail_1 = reference_patterns[1].tail_sequence

        reference_heads = [reference_head_0, reference_head_1]
        reference_tails = [reference_tail_0, reference_tail_1]

        # Recombined patterns are when when head and tail fromr eference patterns are flipped
        recombined_patterns = [
            Pattern(reference_head_0, self.seed_sequence, reference_tail_1),
            Pattern(reference_head_1, self.seed_sequence, reference_tail_0)
        ]

        recombination_occurrences = patterns[recombined_patterns[0]] + patterns[recombined_patterns[1]]
        self.global_data.recombination_occurrences = recombination_occurrences
        self.global_data.recombination_percentage = 100 * (
                recombination_occurrences /
                (recombination_occurrences + patterns[reference_patterns[0]] + patterns[reference_patterns[1]])
        )

        for pattern_tuple in sorted_patterns:
            occurrences = pattern_tuple[1]
            pattern = pattern_tuple[0]
            note = self._get_note(pattern, reference_patterns, recombined_patterns)

            if not note:
                if pattern.head_sequence in reference_heads and pattern.tail_sequence not in reference_tails:
                    note = "Head match"
                elif pattern.head_sequence not in reference_heads and pattern.tail_sequence in reference_tails:
                    note = "Tail match"
            rows.append(
                '{},{},{},{},{},{}\n'.format(pattern.head_sequence, self.seed_sequence,
                                         pattern.tail_sequence, pattern.full_pattern, occurrences, note))

        return rows

    def _extract_pattern(self, read):
        if self.seed_sequence in read:
            seed_sequence_start = read.index(self.seed_sequence)
            # Seed sequence is too close to the beginning of the read
            if seed_sequence_start < params.HEAD:
                return None

            seed_sequence_end = seed_sequence_start + len(self.seed_sequence)
            # Seed sequence is too close to the end of the read
            if len(read) - seed_sequence_end <  params.TAIL:
                return None

            head_sequence = read[seed_sequence_start - params.HEAD: seed_sequence_start]
            assert len(head_sequence) == params.HEAD
            tail_sequence = read[seed_sequence_end: seed_sequence_end + params.TAIL]
            assert len(tail_sequence) == params.TAIL
            return Pattern(head_sequence, self.seed_sequence, tail_sequence)

    @staticmethod
    def _get_note(pattern, reference_patterns, recombined_patterns):
        if pattern in reference_patterns:
            return "Correct sequence"
        if pattern in recombined_patterns:
            return "Clean recombination"
        return ""

    def process_and_write_to_file(self, output_directory):
        rows = self.process()
        input_file_no_path = os.path.basename(self.input_file)
        output_file_no_path = self._get_output_file_name(input_file_no_path)
        output_file = os.path.join(output_directory, output_file_no_path)

        row_pattern = '{},{},{},{},{},{}\n'
        with open(output_file, 'w') as file:
            file.write(row_pattern.format(
                    'File Name',
                    input_file_no_path,
                    '', '', '', ''
            ))

            file.write(row_pattern.format(
                    'Total Reads',
                    self.global_data.total_reads,
                    '', '', '', ''
            ))

            file.write(row_pattern.format(
                    'Occurrences of {} (includes both strands)'.format(self.global_data.seed_sequence_name),
                    self.global_data.seed_occurrences,
                    '', '', '', ''
            ))

            file.write(row_pattern.format(
                    'Occurrences of a Full Sequence (includes both strands)',
                    self.global_data.full_sequence_occurrences,
                    '', '', '', ''
            ))

            file.write(row_pattern.format(
                    'Recombination occurrences',
                    self.global_data.recombination_occurrences,
                    '', '', '', ''
            ))
            file.write(row_pattern.format('% Recombination', self.global_data.recombination_percentage, '', '', '', ''))

            file.write(',,,,,\n')
            file.write(row_pattern.format(
                    'Preceding {} bases'.format(params.HEAD),
                    self.global_data.seed_sequence_name,
                    'Subsequent {} bases'.format(params.TAIL),
                    'Full Sequence',
                    'Occurrences',
                    'Notes'
            )
            )
            for row in rows:
                file.write(row)
        return output_file

    def _get_output_file_name(self, input_file):
        return '{}-{}-{}-{}.csv'.format(
                input_file[:input_file.index('.fastq')],
                params.HEAD,
                self.global_data.seed_sequence_name,
                params.TAIL
        )

##### Command #####
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
def recombination_command():
    """
    Detect recombination in VGS data. FASTQ files containing the data are read from
    an input folder and processed to return the sequences before and after a seed sequence.
    The output are CSV files, one per FASTQ file and per seed sequence.'
    """
    setup_dirs(params)
    input_files = get_fastq_files(params)
    output_files = []
    for input_file in input_files:
        for seed_sequence_name, seed_sequence in params.seed_sequences.items():
            L.info('Processing file: {} with seed sequence: {}...'.format(input_file, seed_sequence_name))
            processor = Processor(input_file, seed_sequence_name, seed_sequence.upper())
            output_file = processor.process_and_write_to_file(params.output_folder)
            output_files.append(output_file)

    L.info('\nDone! Your output is in the following {} files:'.format(str(len(output_files))))
    for output_file in output_files:
        L.info(output_file)


if __name__ == "__main__":
    log_to_stdout(logging.INFO)
    try:
        recombination_command()
    except ValueError as e:
        L.error('\nError: {}'.format(e.message))
        sys.exit(1)
