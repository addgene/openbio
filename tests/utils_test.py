import unittest
from toolkit.utils import get_fastq_files, get_params_for_command
import os


class UtilsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.resources_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'resources')

    def test_get_fastq_files(self):
        expected_files = [os.path.join(self.resources_directory, 'file1.fastq'),
                          os.path.join(self.resources_directory, 'file2.fastq'),
                          os.path.join(self.resources_directory, 'file3.fq')]
        actual_files = get_fastq_files({'input_folder': self.resources_directory})
        self.assertEqual(expected_files, actual_files)

    def test_get_params_for_command(self):
        params = get_params_for_command('count_spacers', param_file=os.path.join(self.resources_directory,
                                                                                 'parameters.yml'))
        self.assertEqual('data/Set_A.csv', params.get('library_file'))
        self.assertEqual('data/', params.get('input_folder'))
        self.assertEqual('analysis', params.get('output_folder'))
        self.assertEqual('tkov3', params.get('sample'))
        self.assertEqual(True, params.get('guide_g'))
        self.assertEqual(30, params.get('key_region_start'))
        self.assertEqual(55, params.get('key_region_end'))
        self.assertEqual('CGAAACACC', params.get('key_region_sequence'))


if __name__ == '__main__':
    unittest.main()
