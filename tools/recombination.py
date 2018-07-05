from __future__ import division
from collections import Counter, OrderedDict

import click
import glob
import os
import sys
from Bio import SeqIO
from parameters import HEAD, TAIL, seed_sequences, input_folder, output_folder


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
        reads = []
        with open(self.input_file, "rU") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                reads.append(str(record.seq))
                reads.append(str(record.seq.reverse_complement()))

        self.global_data.total_reads = len(reads) / 2
        self.global_data.seed_occurrences = sum(1 for read in reads if self.seed_sequence in read)
        patterns = Counter(
                pattern for pattern in (self._extract_pattern(read) for read in reads)
                if pattern
        )
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
                '%s,%s,%s,%s,%s,%s\n' % (pattern.head_sequence, self.seed_sequence,
                                         pattern.tail_sequence, pattern.full_pattern, occurrences, note))

        return rows

    def _extract_pattern(self, read):
        if self.seed_sequence in read:
            seed_sequence_start = read.index(self.seed_sequence)
            # Seed sequence is too close to the beginning of the read
            if seed_sequence_start < HEAD:
                return None

            seed_sequence_end = seed_sequence_start + len(self.seed_sequence)
            # Seed sequence is too close to the end of the read
            if len(read) - seed_sequence_end <  TAIL:
                return None

            head_sequence = read[seed_sequence_start - HEAD: seed_sequence_start]
            assert len(head_sequence) == HEAD
            tail_sequence = read[seed_sequence_end: seed_sequence_end + TAIL]
            assert len(tail_sequence) == TAIL
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

        with open(output_file, 'w') as file:
            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('File Name', input_file_no_path, '', '', '', ''))

            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('Total Reads', self.global_data.total_reads, '', '', '', ''))

            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('Occurrences of ' + self.global_data.seed_sequence_name + ' (includes both strands)',
                          self.global_data.seed_occurrences, '', '', '', ''))

            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('Occurrences of a Full Sequence (includes both strands)',
                          self.global_data.full_sequence_occurrences, '', '', '', ''))

            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('Recombination occurrences', self.global_data.recombination_occurrences, '', '', '', ''))
            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('% Recombination', self.global_data.recombination_percentage, '', '', '', ''))

            file.write(',,,,,\n')
            file.write('%s,%s,%s,%s,%s,%s\n' %
                         ('Preceding ' + str(HEAD) + ' bases',
                          self.global_data.seed_sequence_name,
                          'Subsequent ' + str(TAIL) + ' bases',
                          'Full Sequence',
                          'Occurrences',
                          'Notes'))
            for row in rows:
                file.write(row)
        return output_file

    def _get_output_file_name(self, input_file):
        return (
                input_file[:input_file.index('.fastq')]
                + '-' + str(HEAD) + '-' + self.global_data.seed_sequence_name + '-' + str(TAIL)
                + '.csv'
        )


#### Commmand

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--input',
              default='',
              help='Override folder with input FASTQ files defined in parameters.py.'
                   'Full path or relative to here. Example: /Users/Harry/fastq-data')
@click.option('--output',
              default='',
              help='Override folder for output CSV files defined in parameters.py. '
                   'Full path or relative to here. If the folder doesn\'t exist '
                   'it will be created. Example: /Users/Harry/csv-data')
def recombination_command(input, output):
    """
    Detect recombination in VGS data. FASTQ files containing the data are read from
    an input folder and processed to return the sequences before and after a seed sequence.
    The output are CSV files, one per FASTQ file and per seed sequence.'
    """
    if not input:
        input = input_folder

    if not output:
        output = output_folder

    if not os.path.isdir(input):
        raise ValueError('The folder ' + input + ' does not exist.')

    if not os.path.isdir(output):
        click.echo('\nThe folder ' + output + ' does not exist. Creating it...\n')
        os.mkdir(output)

    # grab all files ending in .fastq
    input_files = [input_file for input_file in glob.glob("%s/*.fastq" % input)]
    if not input_files:
        click.echo('No FASTQ files in folder: ' + input)

    output_files = []
    for input_file in input_files:
        for seed_sequence_name, seed_sequence in seed_sequences.items():
            click.echo('Processing file: ' + input_file + ' with seed sequence: ' + seed_sequence_name + '...')
            processor = Processor(input_file, seed_sequence_name, seed_sequence.upper())
            output_file = processor.process_and_write_to_file(output)
            output_files.append(output_file)

    print('\nDone! Your output is in the following ' + str(len(output_files)) + ' files:')
    for output_file in output_files:
        click.echo(output_file)


if __name__ == "__main__":
    try:
        recombination_command()
    except ValueError as e:
        click.echo('\nError: ' + e.message)
        sys.exit(1)
