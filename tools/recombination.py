from __future__ import division
from collections import Counter, OrderedDict

import click
import glob
import os
import sys
from Bio import SeqIO
from parameters import HEAD, TAIL, seed_sequences, input_folder, output_folder


def split_pattern(pattern, seed_sequence):
    """
    returns head sequence (before seed) and tail sequence (after seed) from a pattern
    """
    assert seed_sequence in pattern
    seed_sequence_start = pattern.index(seed_sequence)
    seed_sequence_end = seed_sequence_start + len(seed_sequence)
    head_sequence = pattern[: seed_sequence_start]
    tail_sequence = pattern[seed_sequence_end: len(pattern)]
    assert len(head_sequence) == HEAD
    assert len(tail_sequence) == TAIL
    return head_sequence, tail_sequence


def process_patterns(patterns, seed_sequence):
    rows = []
    sorted_patterns = patterns.most_common()
    # The two most frequent patterns are likely to be the correct sequences
    reference_patterns = [pattern_tuple[0] for pattern_tuple in sorted_patterns[:2]]
    reference_head_0, reference_tail_0 = split_pattern(reference_patterns[0], seed_sequence)
    reference_head_1, reference_tail_1 = split_pattern(reference_patterns[1], seed_sequence)

    reference_heads = [reference_head_0, reference_head_1]
    reference_tails = [reference_tail_0, reference_tail_1]

    recombined_patterns = [
        reference_head_0 + seed_sequence + reference_tail_1,
        reference_head_1 + seed_sequence + reference_tail_0
    ]

    recombination_occurrences = patterns[recombined_patterns[0]] + patterns[recombined_patterns[1]]
    recombination_percentage = 100 * (
        recombination_occurrences /
        (recombination_occurrences + patterns[reference_patterns[0]] + patterns[reference_patterns[1]])
    )

    for pattern_tuple in sorted_patterns:
        occurrences = pattern_tuple[1]
        pattern = pattern_tuple[0]
        note = get_note(pattern,
                        reference_patterns,
                        recombined_patterns)
        head_sequence, tail_sequence = split_pattern(pattern, seed_sequence)
        if not note:
            if head_sequence in reference_heads and tail_sequence not in reference_tails:
                note = "Head match"
            elif head_sequence not in reference_heads and tail_sequence in reference_tails:
                note = "Tail match"
        rows.append('%s,%s,%s,%s,%s,%s\n' % (head_sequence, seed_sequence, tail_sequence, pattern, occurrences, note))

    return recombination_percentage, rows


def get_note( pattern, reference_patterns , recombined_patterns ):
    if pattern in reference_patterns:
        return "Correct sequence"
    if pattern in recombined_patterns:
        return "Clean recombination"
    return ""

def write_to_file(input_file, output_directory, global_data, csv_rows):
    input_file_no_path = os.path.basename(input_file)
    output_file_no_path = get_output_fn(input_file_no_path, global_data.seed_sequence_name)
    output_file = os.path.join(output_directory, output_file_no_path)

    with open(output_file, 'w') as handle:
        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('File Name', input_file_no_path, '', '', '', ''))

        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('Total Reads', global_data.total_reads, '', '', '', ''))

        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('Occurrences of ' + global_data.seed_sequence_name + ' (includes both strands)',
                      global_data.seed_occurrences, '', '', '', ''))

        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('Occurrences of a Full Sequence (includes both strands)',
                      global_data.full_sequence_occurrences, '', '', '', ''))

        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('% Recombination', global_data.recombination_percentage, '', '', '', ''))


        handle.write(',,,,,\n')
        handle.write('%s,%s,%s,%s,%s,%s\n' %
                     ('Preceding ' + str(HEAD) + ' bases',
                      global_data.seed_sequence_name,
                      'Subsequent ' + str(TAIL) + ' bases',
                      'Full Sequence',
                      'Occurrences',
                      'Notes'))
        for row in csv_rows:
            handle.write(row)
    return output_file


def get_output_fn(fn, seed_sequence_name):
    return fn[:fn.index('.fastq')] + '-' + str(HEAD) + '-' + seed_sequence_name + '-' + str(TAIL) + '.csv'

class GlobalData(object):
    def __init__(self, seed_sequence_name):
        self.seed_sequence_name = seed_sequence_name
        self.total_reads = 0
        self.seed_occurrences = 0
        self.full_sequence_occurrences = 0
        self.recombination_percentage = 0


def process(fn, seed_sequence_name, seed_sequence):
    global_data = GlobalData(seed_sequence_name)

    reads = []
    with open(fn, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            reads.append(str(record.seq))
            reads.append(str(record.seq.reverse_complement()))

    global_data.total_reads = len(reads)/2
    # global_data.reads_with_seed = sum(1 for read in reads if seed_sequence in reads)
    global_data.seed_occurrences = sum(1 for read in reads if seed_sequence in read)
    patterns = Counter( pattern for pattern in (extract_pattern(read, seed_sequence) for read in reads) if pattern)
    global_data.full_sequence_occurrences = sum(patterns.values())
    global_data.recombination_percentage, rows = process_patterns(patterns, seed_sequence)
    return global_data, rows

def extract_pattern(read, seed_sequence):
    if seed_sequence in read:
        seed_sequence_start = read.index(seed_sequence)
        seed_sequence_end = seed_sequence_start + len(seed_sequence)
        pattern = read[seed_sequence_start - HEAD: seed_sequence_end + TAIL]
        expected_length = HEAD + len(seed_sequence) + TAIL
        # this will exclude matches where the seed was found too close to the beginning or end of the read,
        # such that HEAD or TAIL number of bases couldn't be fetched
        if len(pattern) == expected_length:
            return pattern
    return None


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
def process_all(input, output):
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
            global_data, rows = process(input_file, seed_sequence_name, seed_sequence.upper())
            output_file = write_to_file(input_file, output, global_data, rows)
            output_files.append(output_file)

    print('\nDone! Your output is in the following ' + str(len(output_files)) + ' files:')
    for output_file in output_files:
        click.echo(output_file)


if __name__ == "__main__":
    try:
        process_all()
    except ValueError as e:
        click.echo('\nError: ' + e.message)
        sys.exit(1)
