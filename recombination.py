from Bio import SeqIO

TAIL = 60
HEAD = 60

signatures = {
  "loxP":    "TAAAGCGGCCGTCGACGATAT",
  "lox2272": "GGTCGATGGTGAAGCATTGGT",
}


def split_pattern(pattern, signature):
  """
  returns head sequence (before signature) and tail sequence (after signature) from a pattern
  """
  assert signature in pattern
  signature_start = pattern.index(signature)
  signature_end = signature_start + len(signature)
  head_sequence = pattern[: signature_start]
  tail_sequence = pattern[signature_end : len(pattern)]
  assert len(head_sequence) == HEAD
  assert len(tail_sequence) == TAIL
  return head_sequence, tail_sequence

def process(fn, signature_name, signature_value):
  reads = []
  with open(fn, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      reads.append(str(record.seq))
      reads.append(str(record.seq.reverse_complement()))

  patterns = {}
  for read in reads:
    if signature_value in read:
      signature_start = read.index(signature_value)
      signature_end = signature_start + len(signature_value)
      pattern = read[signature_start - HEAD : signature_end + TAIL]
      expected_length = HEAD + len(signature_value) + TAIL
      # this excludes cases where the signature was found too close to the beginning or end of the read, such that
      # HEAD or TAIL number of elements couldn't be fetched
      if len(pattern) == expected_length:
        if pattern not in patterns:
          patterns[pattern] = 0
        patterns[pattern] += 1

  write_to_file(fn, patterns, signature_name, signature_value)


def write_to_file(fn, patterns, signature_name, signature_value):

  output_fn = fn[:fn.index('.fastq')] + '-' + signature_name + '.csv'
  with open(output_fn, 'w') as handle:
    handle.write('%s,%s,%s,%s\n' % ('File Name: ' + fn, '', '', ''))
    handle.write('%s,%s,%s,%s\n' %
                 ('Head Sequence of size ' + str(HEAD),
                  signature_name,
                  'Tail Sequence of size ' + str(TAIL),
                  'Number of occurrences'))
    sorted_patterns = sorted(patterns.keys(), key=lambda k: -patterns[k])
    for pattern in sorted_patterns:
      occurrences = patterns[pattern]
      head_sequence, tail_sequence = split_pattern(pattern, signature_value)
      handle.write('%s,%s,%s,%s\n' % (head_sequence, signature_value, tail_sequence, occurrences))

if __name__ == "__main__":
  import sys
  import glob
  directory = sys.argv[1]
  # grab all files ending in .fastq
  fns = [fn for fn in glob.glob("%s/*.fastq" % directory)]
  for fn in fns:
    for signature_name, signature_value in signatures.items():
      print('Processing file: ' + fn + ' and seed: ' + signature_name + '...')
      process(fn, signature_name, signature_value)

