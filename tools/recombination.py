import click
import glob
import os
import sys
from Bio import SeqIO
from parameters import HEAD, TAIL, seed_sequences


def split_pattern(pattern, seed_sequence):
  """
  returns head sequence (before seed) and tail sequence (after seed) from a pattern
  """
  assert seed_sequence in pattern
  seed_sequence_start = pattern.index(seed_sequence)
  seed_sequence_end = seed_sequence_start + len(seed_sequence)
  head_sequence = pattern[: seed_sequence_start]
  tail_sequence = pattern[seed_sequence_end : len(pattern)]
  assert len(head_sequence) == HEAD
  assert len(tail_sequence) == TAIL
  return head_sequence, tail_sequence


def patterns_to_csv_rows(patterns, seed_sequence_name, seed_sequence):
    rows = []
    sorted_patterns = sorted(patterns.keys(), key=lambda k: -patterns[k])
    for pattern in sorted_patterns:
      occurrences = patterns[pattern]
      head_sequence, tail_sequence = split_pattern(pattern, seed_sequence)
      rows.append('%s,%s,%s,%s\n' % (head_sequence, seed_sequence, tail_sequence, occurrences))
    return rows


def write_to_file(input_file, output_directory, csv_rows, seed_sequence_name):
    input_file_no_path = os.path.basename(input_file)
    output_file_no_path = get_output_fn(input_file_no_path, seed_sequence_name)
    output_file = os.path.join(output_directory, output_file_no_path)

    with open(output_file, 'w') as handle:
        handle.write('%s,%s,%s,%s\n' % ('File Name: ' + input_file_no_path, '', '', ''))
        handle.write('%s,%s,%s,%s\n' %
                 ('Preceding ' + str(HEAD) + ' bases',
                  seed_sequence_name,
                  'Subsequent ' + str(TAIL) + ' bases',
                  'Occurrences'))
        for row in csv_rows:
            handle.write(row)
    return output_file


def get_output_fn(fn, seed_sequence_name):
    return fn[:fn.index('.fastq')] + '-' + str(HEAD) +'-' + seed_sequence_name + '-' + str(TAIL) + '.csv'


def process(fn, seed_sequence_name, seed_sequence):
  reads = []
  with open(fn, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      reads.append(str(record.seq))
      reads.append(str(record.seq.reverse_complement()))

  patterns = {}
  for read in reads:
    if seed_sequence in read:
      seed_sequence_start = read.index(seed_sequence)
      seed_sequence_end = seed_sequence_start + len(seed_sequence)
      pattern = read[seed_sequence_start - HEAD : seed_sequence_end + TAIL]
      expected_length = HEAD + len(seed_sequence) + TAIL
      # this will exclude matches where the seed was found too close to the beginning or end of the read, such that
      # HEAD or TAIL number of bases couldn't be fetched
      if len(pattern) == expected_length:
        if pattern not in patterns:
          patterns[pattern] = 0
        patterns[pattern] += 1

  return patterns_to_csv_rows(patterns, seed_sequence_name, seed_sequence)


#### Commmand

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--input',
              default='data',
              prompt='Folder with input files. Full path or relative to here. Hit ENTER for default shown in brackets',
              help='Folder with input FASTQ files. Full path or relative to here. Defaults to the relative folder "data". '
                   'Example: /Users/Harry/fastq-data')
@click.option('--output',
              default='data',
              prompt='Folder for output files. Full path or relative to here. Hit ENTER for default shown in brackets',
              help='Folder for output CSV files. Full path or relative to here. If the folder doesn\'t exist '
                   'it will be created. Defaults to the relative folder "data". '
                   'Example: /Users/Harry/csv-data')
def process_all(input, output):
    """
    Detect recombination in VGS data. FASTQ files containing the data are read from
    an input directory and processed to return the sequences before and after a seed sequence.
    The output are CSV files, one per FASTQ file and per seed sequence.'
    """
    if not os.path.isdir( input ):
        raise ValueError('The folder ' + input + ' does not exist.')

    if not os.path.isdir( output ):
        click.echo('\nThe folder ' + output + ' does not exist. Creating it....\n')
        os.mkdir( output )

    # grab all files ending in .fastq
    input_files = [input_file for input_file in glob.glob("%s/*.fastq" % input)]
    if not input_files:
        click.echo('No FASTQ files in folder: ' + input)

    output_files = []
    for input_file in input_files:
        for seed_sequence_name, seed_sequence in seed_sequences.items():
            click.echo('Processing file: ' + input_file + ' with seed sequence: ' + seed_sequence_name + '...')
            rows = process(input_file, seed_sequence_name, seed_sequence)
            output_file = write_to_file(input_file, output, rows, seed_sequence_name)
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

