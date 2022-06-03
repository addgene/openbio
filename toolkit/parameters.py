########## recombination.py parameters
class Recombination_Parameters(object):

    # Change these two values to the folders you prefer - use an absolute path e.g. /Users/Harry/fastq-data and
    # /Users/Harry/csv-data or a path relative to the tools directory.
    # You may use the same folder for input and output.
    input_folder = "data"
    output_folder = "data"

    # The number of bases to retrieve before the seed sequence
    HEAD = 20

    # The number of bases to retrieve after the seed sequences
    TAIL = 20

    seed_sequences = {
        "LoxP": "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
        "Lox2272": "ATAACTTCGTATAGGATACTTTATACGAAGTTAT",
        "LoxN":"ataacttcgtataaggtatactatacgaagttat",
        "FRT":"GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC",
        "F5":"GAAGTTCCTATTCTTCAAAAGGTATAGGAACTTC",
        "Lox71_66":"GCATACAT",
    }




########## serotypes.py parameters
class Serotypes_Parameters(object):

    # Change these two values to the folders you prefer - use an absolute path e.g. /Users/Harry/fastq-data and
    # /Users/Harry/csv-data or a path relative to the tools directory.
    # You may use the same folder for input and output.
    input_folder = "data"
    output_folder = "data"

    # These are the signatures that will be matched. The first part is the name, the part in brackets contains the
    # actual signatures, separated by a comma (each serotype can have multiple signatures)
    signatures = {
        "AAV1": [
            "AGTGCTTCAACGGGGGCCAG",
            "GGGCGTGAATCCATCATCAACCCTGG",
            "CCGGAGCTTCAAACACTGCATTGGACAAT"
        ],
        "AAV2": [
            "AGGCAACAGACAAGCAGCTACC",
            "AACAGACAAGCAGCTACCGCA"
        ],
        "AAV5": [
            "TCCAAGCCTTCCACCTCGTCAGACGCCGAA",
            "CACCAACAACCAGAGCTCCACCACTG",
            "GCCCGTCAGCAGCTTCATC"
        ],
        "AAV7": [
            "AGTGAAACTGCAGGTAGTACC"
        ],
        "AAV8": [
            "GCAAAACACGGCTCCTCAAAT",
            "CAGCAAGCGCTGGAACCCCGAGATCCAGTA",
            "AAATACCATCTGAATGGAAGAAATTCATTG",
            "CGTGGCAGATAACTTGCAGC",
            "ATCCTCCGACCACCTTCAACC"
        ],
        "AAV9": [
            "AGTGCCCAAGCACAGGCGCA",
            "ATCTCTCAAAGACTATTAAC",
            "GGCGAGCAGTCTTCCAGGCA"
        ],
        "AAVrh10": [
            "CTACAAATCTACAAATGTGGACTTTG"
        ],
        "PHPeB": [
            "CTTTGGCGGTGCCTTTTAAGGCACAGGCGCAGA"
        ],
        "PHPs": [
            "AGGCGGTTAGGACGTCTTTGGCACAGGCGCAGA"
        ],
        "AAVrg": [
            "TAGCAGACCAAGACTACACAAAAACTGCT"
        ],
    }
