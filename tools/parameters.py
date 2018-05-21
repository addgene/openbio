# recombination.py parameters

# The number of bases to retrieve after the seed sequences
TAIL = 60
# The number of bases to retrieve before the seed sequence
HEAD = 60

seed_sequences = {
  "before_loxP":    "TAAAGCGGCCGTCGACGATAT",
  "before_lox2272": "GGTCGATGGTGAAGCATTGGT",
}

# Change these to the folders you prefer - use an absolute path if necessary, e.g. /Users/Harry/fastq-data and
# /Users/Harry/csv-data
#

input_folder = "data"
output_folder = "data"