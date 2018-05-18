from Bio import SeqIO
from Bio.Seq import Seq

TAIL = 60

# e.g. sequence after mCherry
signatures = {
  "before_loxP":    "TAAAGCGGCCGTCGACGATAT",
  "before_lox2272": "GGTCGATGGTGAAGCATTGGT",
}


def split_pattern(p, s):
  assert p.startswith(s)
  return p[0:len(s)], p[len(s):]


def process(fn, signature):
  reads = []
  with open(fn, "rU") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      reads.append(str(record.seq))
      reads.append(str(record.seq.reverse_complement()))

  patterns = {}
  for reads in reads:
    if signature in reads:
      pattern = reads[reads.index(signature):reads.index(signature)+len(signature)+TAIL]
      if len(pattern) == len(signature)+TAIL:
        if pattern not in patterns:
          patterns[pattern] = 0
        patterns[pattern] += 1

  print fn
  for p in sorted(patterns.keys(), key=lambda k: -patterns[k]):
    i = patterns[p]
    a, b = split_pattern(p, signature)
    print "%s %s: %s" % (a, b, i)


for n,s in signatures.iteritems():
  print n
  process("x/A11984_sW0158_H09_R_001.fastq", s)
