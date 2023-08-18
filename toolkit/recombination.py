from __future__ import division
from collections import Counter

import click
import logging
import os

from pandas import DataFrame

from utils import (
    get_fastq_files,
    get_params_for_command,
    read_fastq_file_with_reverse_complements,
    setup_dirs
)


L = logging.getLogger(__name__)

#### Command ####
@click.command()
def recombination():
    """
    Detect recombination in VGS data. FASTQ files containing the data are read from
    an input folder and processed to return the sequences before and after a seed sequence.
    The output are CSV files, one per FASTQ file and per seed sequence.'
    """
    params = get_params_for_command('recombination')
    if 'seed_sequences' not in params:
        raise ValueError('Could not find "seed_sequences" parameter in parameters.yml')

    setup_dirs(params)
    input_files = get_fastq_files(params)
    output_files = []
    for input_file in input_files:
        for seed_sequence_name, seed_sequence in params.get('seed_sequences').items():
            L.info('Processing file: {} with seed sequence: {}...'.format(input_file, seed_sequence_name))
            processor = Processor(input_file, seed_sequence_name, seed_sequence.upper())
            try:
                output_file = processor.process_and_write_to_file(params.get('output_folder'))
                if output_file:
                    output_files.append(output_file)
            except Exception as e:
                L.error('    >> Error: {}. Ignoring file...'.format(e))

    L.info('\nDone! Your output is in the following {} files:'.format(str(len(output_files))))
    for output_file in output_files:
        L.info(output_file)

##################


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
        self.params = get_params_for_command('recombination')
        self.head = self.params.get('HEAD', 0)
        self.tail = self.params.get('TAIL', 0)

        self.input_file = input_file
        self.seed_sequence = seed_sequence
        self.global_data = GlobalData( seed_sequence_name)

    def process(self):
        """
        Processes FASTQ input file and returns a list of CSV rows
        """
        reads = read_fastq_file_with_reverse_complements(self.input_file)

        self.global_data.total_reads = len(reads) / 2
        self.global_data.seed_occurrences = sum(1 for read in reads if self.seed_sequence in read)
        extracted_patterns = [self._extract_pattern(read) for read in reads]

        # The file did not contain the seed
        if not any(extracted_patterns):
            L.warning('    >> Error: the file does not contain the seed {}. Ignoring file...'.format(self.seed_sequence))
            return []

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

        # Recombined patterns are when head and tail from reference patterns are flipped
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
                    (pattern.head_sequence, self.seed_sequence, pattern.tail_sequence,
                     pattern.full_pattern, occurrences, note)
            )

        return rows

    def _extract_pattern(self, read):
        if self.seed_sequence in read:
            seed_sequence_start = read.index(self.seed_sequence)
            # Seed sequence is too close to the beginning of the read
            if seed_sequence_start < self.head:
                return None

            seed_sequence_end = seed_sequence_start + len(self.seed_sequence)
            # Seed sequence is too close to the end of the read
            if len(read) - seed_sequence_end <  self.tail:
                return None

            head_sequence = read[seed_sequence_start - self.head: seed_sequence_start]
            assert len(head_sequence) == self.head
            tail_sequence = read[seed_sequence_end: seed_sequence_end + self.tail]
            assert len(tail_sequence) == self.tail
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

        if not rows:
            return ''

        preamble = []
        preamble.append(('File Name', input_file_no_path))
        preamble.append(('Total Reads', self.global_data.total_reads))
        preamble.append(('Occurrences of {} (includes both strands)'.format(self.global_data.seed_sequence_name),
                    self.global_data.seed_occurrences))
        preamble.append(('Occurrences of a Full Sequence (includes both strands)',
                    self.global_data.full_sequence_occurrences))
        preamble.append(('Recombination occurrences',
                    self.global_data.recombination_occurrences))
        preamble.append(('% Recombination', self.global_data.recombination_percentage))
        preamble.append(('',''))

        preamble.append(('Preceding {} bases'.format(self.head),
                    self.global_data.seed_sequence_name,
                    'Subsequent {} bases'.format(self.tail),
                    'Full Sequence',
                    'Occurrences',
                    'Notes'))

        df = DataFrame(preamble + rows)
        df.to_csv(output_file, header=False, index=False)

        return output_file

    def _get_output_file_name(self, input_file):
        return '{}-{}-{}-{}.csv'.format(
                input_file[:input_file.index('.fastq')].replace(' ', '_'),
                self.head,
                self.global_data.seed_sequence_name,
                self.tail
        )
