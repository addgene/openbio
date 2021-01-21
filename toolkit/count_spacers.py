"""
    Analyze sequencing data for sgRNA library distribution.

    The original version of this tool was described in:
    Joung, Julia et al.
    “Genome-scale CRISPR-Cas9 knockout and transcriptional activation screening.”
    Nature protocols vol. 12,4 (2017): 828-863. doi:10.1038/nprot.2017.016
    See:
    https://github.com/fengzhanglab/Screening_Protocols_manuscript/blob/master/count_spacers.py
"""

import click
import csv
import io
import logging
import numpy
import os
import re
from collections import OrderedDict
from datetime import datetime
from logging import Logger

from utils import get_params_for_command, log_to_stdout, read_fastq_file, get_fastq_files

L: Logger = logging.getLogger(__name__)


#### Command ####

@click.command()
@click.option('-i', '--input', 'library_file', help='Name of file containing library sequences')
@click.option('-fd', '--fastq-folder', 'input_folder', help='Name of directory containing FASTQ files')
@click.option('-od', '--output-folder', 'output_folder', help='Name of directory for output')
@click.option('-g', '--guide-g', 'guide_g', flag_value=True, help='Presence of guanine before spacer', default=True)
@click.option('-u', '--unattended', 'unattended', flag_value=True, help='Skip user prompts', default=False)
@click.option('-p', '--param-file', 'param_file', help='Name of yaml file with configuration parameters',
			  default='parameters.yml')

def count_spacers(library_file, input_folder, output_folder, guide_g, unattended, param_file):
	"""
    Analyze sequencing data for sgRNA library distribution.
    Parameters can be configured in the  `parameters.yml` file, through command line options,
    or through the interactive prompts.
    Values passed in the command line override values in `parameters.yml`.
    """
	config = Config(library_file, input_folder, output_folder, guide_g, unattended, param_file)
	spacer_counter = SpacerCounter(config)
	spacer_counter.execute()
	L.info('Done!')

##################


class SpacerCounter():
	def __init__(self, config):
		self.config = config

	def execute(self):
		"""
		creates a file with guide counts from fastq_file
		dictionary: guide sequence as key, guide count as entry
		"""
		L.info('\nProcessing...')

		# Add 'G' to key sequence if included in library
		if self.config.guide_g:
			self.config.key_region_sequence += 'G'

		# open library sequences and initiate dictionary of read counts for each guide
		library_sequences = self._get_library_sequences()

		# read fastq files in the specified directory
		files_to_process = get_fastq_files({'input_folder': self.config.input_folder})

		for file in files_to_process:
			self._process_fastq_file(file, library_sequences)

	def _process_fastq_file(self, file, library_sequences):
		num_reads = 0  # total number of reads processed
		perfect_matches = 0  # guides with perfect match to library
		non_perfect_matches = 0  # number of guides without a perfect match to the library
		key_not_found = 0  # count of reads where key was not found
		base_filename = os.path.basename(file)
		base_filename = os.path.splitext(base_filename)[0]
		# open fastq file and obtain reads
		reads = read_fastq_file(file)

		# process reads in fastq file
		for read_sequence in reads:  # contains the seq and Qscore etc.
			num_reads += 1
			key_region = read_sequence[self.config.key_region_start:self.config.key_region_end]
			key_index = key_region.find(self.config.key_region_sequence)
			if key_index >= 0:
				start_index = key_index + self.config.key_region_start + len(self.config.key_region_sequence)
				guide = read_sequence[start_index:(start_index + 20)]
				if guide in library_sequences:
					library_sequences[guide] += 1
					perfect_matches += 1
				else:
					non_perfect_matches += 1
			else:
				key_not_found += 1

		# create ordered dictionary with guides and respective counts and output as a csv file
		sorted_library_sequences = OrderedDict(sorted(library_sequences.items(), key=lambda t: t[0]))
		self._write_counts(sorted_library_sequences, base_filename)

		# percentage of guides that matched perfectly
		all_matches = perfect_matches + non_perfect_matches
		# avoid division by zero
		if all_matches:
			percent_matched = round(perfect_matches / all_matches * 100, 6)
		else:
			percent_matched = 0

		# percentage of undetected guides with no read counts
		all_guides = len(library_sequences.values())
		guides_with_reads = numpy.count_nonzero(library_sequences.values())
		guides_no_reads = all_guides - guides_with_reads
		percent_no_reads = round(guides_no_reads / all_guides * 100, 6)
		# skew ratio of top 10% to bottom 10% of guide counts
		top_10 = numpy.percentile(list(library_sequences.values()), 90)
		bottom_10 = numpy.percentile(list(library_sequences.values()), 10)
		if top_10 != 0 and bottom_10 != 0:
			skew_ratio = top_10 / bottom_10
		else:
			skew_ratio = 'Not enough perfect matches to determine skew ratio'

		# Write analysis statistics to statistics.txt
		statistics = []
		statistics.append('Number of perfect guide matches: {}'.format(perfect_matches))
		statistics.append('Number of nonperfect guide matches: {}'.format(non_perfect_matches))
		statistics.append('Number of reads where key was not found: {}'.format(key_not_found))
		statistics.append('Number of reads processed: {}'.format(num_reads))
		statistics.append('Percentage of guides that matched perfectly: {}'.format(percent_matched))
		statistics.append('Percentage of undetected guides: {}'.format(percent_no_reads))
		statistics.append('Skew ratio of top 10% to bottom 10%: {}'.format(skew_ratio))
		statistics = '\n'.join(statistics)
		self._write_statistics(statistics, base_filename)
		L.info(statistics)

		# Write parameters and user input values to log.txt
		self._write_log(base_filename)

	# TODO rewrite with pandas
	def _get_library_sequences(self):
		# Open library sequences file and initiate dictionary of read counts for each guide
		with io.open(self.config.library_file, newline=None) as infile:
			reader = csv.reader(infile)
			return {rows[0]: 0 for rows in reader}

	def _get_output_file_name(self, base_filename, suffix):
		local_name = '{}{}'.format(base_filename, suffix)
		return os.path.join(self.config.output_folder, local_name)

	# TODO change to pandas
	def _write_counts(self, sorted_library_sequences, base_filename):
		filename = self._get_output_file_name(base_filename, '_library_counts.csv')
		L.info('Generating {}...'.format(filename))
		with io.open(filename, 'w') as f:
			mywriter = csv.writer(f, delimiter=',')
			for guide in sorted_library_sequences:
				count = sorted_library_sequences[guide]
				mywriter.writerow([guide, count])

	def _write_statistics(self, statistics: str, base_filename):
		filename = self._get_output_file_name(base_filename, '_statistics.txt')
		L.info('Generating {}...'.format(filename))
		with io.open(filename, 'w') as f:
			f.write(statistics)

	def _write_log(self, base_filename):
		filename = self._get_output_file_name(base_filename, '_log.txt')
		L.info('Generating {}...'.format(filename))
		with io.open(filename, 'w') as f:
			f.write('Date and time: {}\n{}'.format(datetime.now(), self.config))


class Config():
	"""
	This class contains all parameters necessary to run `count_spacers`. It handles resolving their values
	(i.e. options given in the config_file are overriden by options given in the command_line) and validating them.
	It provides an interactive mode where the user is prompted for new values.
	"""
	def __init__(self, library_file, input_folder, output_folder, guide_g, unattended, param_file):
		params = get_params_for_command('count_spacers', param_file=param_file)
		self.input_folder = input_folder or params.get('input_folder', '')
		self.library_file = library_file or params.get('library_file', '')
		self.output_folder = output_folder or params.get('output_folder', '')
		self.guide_g = guide_g or params.get('guide_g', '')
		self.sample = params.get('sample', '')
		self.key_region_start = params.get('key_region_start', '')
		self.key_region_end = params.get('key_region_end', '')
		self.key_region_sequence = params.get('key_region_sequence', '')
		# Only prompt the user if in interactive mode
		if not unattended:
			self._ask_if_changes()
		params_ok = self._validate_params()
		if not params_ok and not unattended:
			click.confirm('Try again?', default=True, abort=True)
			self._prompt()

		if not os.path.isdir(self.output_folder):
			L.info('\nThe folder ' + self.output_folder + ' does not exist. Creating it...\n')
			os.mkdir(self.output_folder)

	def _validate_params(self):
		"""
		Verifies that all parameters without explicit defaults were defined. Verifies that files exist.
		"""
		errors = []

		if not self.library_file:
			errors.append('Library file must be specified.')
		if not self.input_folder:
			errors.append('FASTQ folder must be specified.')

		elif not os.path.isdir(self.input_folder):
			errors.append('FASTQ folder specified is not a directory.')

		if self.library_file:
			file_exists = os.path.exists(self.library_file)
			if not file_exists:
				errors.append('Could not find file: {}'.format(self.library_file))

		if not self.sample:
			errors.append('Sample name must be specified.')

		if not self.key_region_start:
			errors.append('Key region start must be specified.')

		if not self.key_region_end:
			errors.append('Key region end must be specified.')

		if not self.key_region_sequence:
			errors.append('Key region sequence must be specified.')

		if errors:
			L.warning('Please correct the following errors:\n{}'.format('\n'.join(errors)))
			return False
		return True

	def _ask_if_changes(self):
		"""
		Asks the user if changes are needed after listing the values found so far
		"""
		click.secho('\nThe following parameter values were found\n', bold=True)
		self._print_config()
		if click.confirm('\nWould you like to change any of these values?', default=False):
			self._prompt()

	def _prompt(self):
		"""
		Prompts the user to enter individual parameter values or accept existing values
		"""
		click.echo(
				'Please specify the'
				' values when prompted. '
				'The current value is shown in brackets, hit enter to keep it. At any moment, use Ctrl-c to exit.'
		)
		self.library_file = click.prompt('Library file name', type=str, default=self.library_file)
		self.input_folder = click.prompt('FASTQ folder name', type=str, default=self.input_folder)
		self.output_folder = click.prompt('Output folder', type=str, default=self.output_folder)
		self.guide_g = click.prompt('Guide G?', type=bool, default='yes' if self.guide_g else 'no')
		self.sample = click.prompt('Sample name', type=str, default=self.sample)
		self.key_region_start = click.prompt(
				'Start index of key region, integer numbers only', type=int, default=self.key_region_start
		)
		self.key_region_end = click.prompt(
				'End index of key region, integer numbers only', type=int, default=self.key_region_end
		)

		# identifies sequence before guide to determine guide position while only accepting ATC
		key = None
		while not key:
			key = click.prompt(
					'Key region sequence, A,T,C and G only: ',
					value_proc=self._validate_sequence,
					default=self.key_region_sequence
			)
			if not key:
				click.echo('This only accepts characters A,T, C, and G, please try again.')
			else:
				self.key_region_sequence=key

		self._confirm()

	def _confirm(self):
		"""
		Asks the user to confirm existing values
		"""
		click.secho('\nThe following parameter values will be used:\n', bold=True)
		self._print_config()
		next = click.prompt(
				'\nContinue with these values?',
				type=click.Choice(['y', 'n', 'exit'], case_sensitive=False),
				default='y'
		)
		if next == 'y':
			return
		if next == 'n':
			self._prompt()
		elif next == 'exit':
			click.get_current_context().exit(0)

	def _validate_sequence(self, sequence):
		"""
		Checks that the sequence provided contains only ATCG, case insensitive
		"""
		if re.match(r'^([atgcATGC])*$', sequence):
			return sequence.upper()
		return None

	def _print_config(self):
		"""
		Prints the current configuration
		"""
		click.echo(str(self))

	def __str__(self):
		"""
		Serializes the current configuration to a string
		"""
		config = []
		config.append('Library file: {}'.format(self.library_file))
		config.append('FASTQ folder: {}'.format(self.input_folder))
		config.append('Output folder: {}'.format(self.output_folder))
		config.append('Guide G? {}'.format('yes' if self.guide_g else 'no'))
		config.append('Sample name: {}'.format(self.sample))
		config.append('Start index of key region: {}'.format(self.key_region_start))
		config.append('End index of key region: {}'.format(self.key_region_end))
		config.append('Key region sequence: {}'.format(self.key_region_sequence))
		return ('\n'.join(config))


if __name__ == "__main__":
	log_to_stdout(logging.INFO)
	try:
		count_spacers()
	except Exception as e:
		L.exception(e)
		raise SystemExit('\nCommand terminated with an error: ' + str(e))
