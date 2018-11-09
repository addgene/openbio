########## recombination.py parameters
class Recombination_Parameters(object):

    # Change these two values to the folders you prefer - use an absolute path e.g. /Users/Harry/fastq-data and
    # /Users/Harry/csv-data or a path relative to the tools directory.
    # You may use the same folder for input and output.
    input_folder = "/Users/Daniela/Home/Dev/Workspace/research/data/bla"
    output_folder = "/Users/Daniela/Home/Dev/Workspace/research/data/bla"

    # The number of bases to retrieve before the seed sequence
    HEAD = 10

    # The number of bases to retrieve after the seed sequences
    TAIL = 10

    seed_sequences = {
        "loxP": "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
        "lox2272": "ATAACTTCGTATAGGATACTTTATACGAAGTTAT",
    }


########## serotype_report.py parameters
class Serotype_Report_Parameters(object):

    # Change these two values to the folders you prefer - use an absolute path e.g. /Users/Harry/fastq-data and
    # /Users/Harry/csv-data or a path relative to the tools directory.
    # You may use the same folder for input and output.
    input_folder = "/Users/Daniela/Home/Dev/Workspace/research/data/bla"
    output_folder = "/Users/Daniela/Home/Dev/Workspace/research/data/bla"

    # These are the signatures that will be matched. The first part is the name, the second part is the actual signature
    signatures = {
        "AAV1": "AGTGCTTCAACGGGGGCCAG",
        "AAV2": "AACAGACAAGCAGCTACCGCA",
        "AAV5": "TCCAAGCCTTCCACCTCGTCAGACGCCGAA",
        "AAV7": "AGTGAAACTGCAGGTAGTACC",
        "AAV8": "GCAAAACACGGCTCCTCAAAT",
        "AAV9": "AGTGCCCAAGCACAGGCGCA",
        "AAVrh10": "CTACAAATCTACAAATGTGGACTTTG",
        "PHPeB": "CTTTGGCGGTGCCTTTTAAGGCACAGGCGCAGA",
        "PHPS": "AGGCGGTTAGGACGTCTTTGGCACAGGCGCAGA",
        "AAVrg": "TAGCAGACCAAGACTACACAAAAACTGCT",
    }