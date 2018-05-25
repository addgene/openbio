# recombination.py parameters

# The number of bases to retrieve after the seed sequences
TAIL = 10
# The number of bases to retrieve before the seed sequence
HEAD = 10

seed_sequences = {
    "loxP": "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
    "lox2272": "ATAACTTCGTATAGGATACTTTATACGAAGTTAT",
}

# Change these two values to the folders you prefer - use an absolute path e.g. /Users/Harry/fastq-data and
# /Users/Harry/csv-data or a path relative to the tools directory.
# You may use the same folder for input and output.

input_folder = "data"
output_folder = "data"