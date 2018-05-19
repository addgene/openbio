#!/usr/bin/env python

import collections
import csv
import re
import random

def get_molecules():
  with open(inputseqs, 'rb') as readinput:
    reader = csv.reader(readinput.read().splitlines())
    seqsdict = collections.OrderedDict(reader)
  return seqsdict

def chunk(s, chunklen):
  s = [s[i:i+chunklen] for i in xrange(0, len(s)) if i == 0 or i+chunklen < len(s)]
  #[a.lower() for a in s] #convert all chunks to lowercase, but not working
  return s

#gets reverse complements
def rc(s):
  _r = dict(A="T",T="A",G="C",C="G")
  return "".join([_r[c.upper()] for c in s][::-1])


def uniquefy(molecules, chunklen, samplesize=None):

  # this is a little memory intensive
  counts = {}

  chunked = collections.OrderedDict()
  for molecule_id, molecule in molecules.iteritems():
    name = molecule_id
    sequence = molecule
#get chunks from hard-coded molecules
    chunks_f = chunk(sequence, chunklen)
    chunks_r = chunk(rc(sequence), chunklen)
    chunks = chunks_f+chunks_r

    # keep counts of each substring
    for s in chunks:
      if s not in counts:
        counts[s] = 0
      counts[s] += 1

    chunked[molecule_id] = { "id": molecule_id, "name": name, "chunks": chunks }

  # find unique chunks

  all_unique = {}
  len_unique_chunks = []
#put len(unique_chunks) in list
  for molecule_id, molecule_chunks in chunked.iteritems():
    chunks = molecule_chunks["chunks"]
    unique_chunks = [s for s in chunks if counts[s] == 1]
    if len(unique_chunks) == 0:
      raise Exception("m%s has no unique sequence across this batch of molecules" % molecule_id)
    #print "m%s: %s unique chunks" % (molecule_id, len(unique_chunks))
    len_unique_chunks.append(len(unique_chunks))

    if samplesize is not None:
      unique_chunks = random.sample(unique_chunks, samplesize)

    for s in unique_chunks:
      assert s not in all_unique
      all_unique[s] = molecule_id
  return (all_unique, chunked, len_unique_chunks)

def get_fastq_reads(fn):
  f = open(fn, 'rb')
  lines = f.readlines()
  return [line.strip().upper() for i,line in enumerate(lines) if i%4 == 1]

def counted(reads):
  counts = {}
  for r in reads:
    if r not in counts:
      counts[r] = 0
    counts[r] += 1
  return counts

def detect(fn, all_unique, molecules, chunklen, samplesize=None):
  reads = get_fastq_reads(fn)
  # sample reads
  if samplesize and samplesize < len(reads):
#random.sample(population, k) returns a k length list of unique elements chosen from population sequence
    reads = random.sample(reads, samplesize)
  chunked_reads = []
#get chunks from reads in fastq files
  for s in reads:
    chunked_reads.extend(chunk(s, chunklen))
  read_counts = counted(chunked_reads)

  print "%s reads, %s chunks, %s unique" % (len(reads), len(chunked_reads), len(read_counts.keys()))

  identified = {}
  for read_chunk, count in read_counts.iteritems():
    if read_chunk in all_unique:
      molecule_id = all_unique[read_chunk]

      if molecule_id not in identified:
        name = molecules[molecule_id]["name"]
        identified[molecule_id] = { "name": name, "count": count }
        #print "  m%s: %s" % (molecule_id, name)
        #print read_chunk
      else:
        identified[molecule_id]["count"] += count

  return identified

#gets sample id from file name, where sample id is whatever comes before the first underscore
def fn_sample_id(fn):
  import os.path
  h,t = os.path.split(fn)
  tokens = re.split("_", t)
  sample_id = tokens[0]
  return sample_id

def analyze(fns, verbose=False):
  #puts unique sample ids into list and counts number of samples in that list
  sample_ids = []
  for fn in fns:
    sample_id = fn_sample_id(fn)
    sample_ids.append(sample_id)
  sample_ids = list(set(sample_ids))
  print "%s samples" % len(sample_ids)

  molecules = get_molecules()
  assert len(molecules) > 0
  print "%s molecules" % len(molecules)

  all_unique, chunked, len_unique_chunks = uniquefy(molecules, 50)

  #normalized list to hold values of normalized counts  
  normalized_list = [] 

  for i, fn in enumerate(fns):
    if verbose:
      print "File: %s" % fn
    identified = detect(fn, all_unique, chunked, 50)
    print "Number of Matches"
    total = 0
    norm_list_num = 0
    seq_place = 0
    for id, m in identified.iteritems():
      if verbose:
        name = m["name"]
        #divides # matches by # unique chunks found for each molecule
        normalized_count = (m["count"]/len_unique_chunks[seq_place])
        normalized_list.append(normalized_count)
        if normalized_count > 0:
          print "%s: %s" % (name, normalized_count)
      total += normalized_count
      seq_place += 1
    print "Percent Matches"
    if total > 0:
      for id, m in identified.iteritems():
        name = m["name"]
        norm_list_num += 1
        percent = (normalized_list[norm_list_num-1]*100.0/total)
        if percent > 0:
          print "%s: %s" % (name, percent)
      #clear normalized list for next file
    normalized_list = []

  
if __name__ == "__main__":
  import sys
  import glob
  dirfn = sys.argv[1]
  inputseqs = sys.argv[2]
#look for files ending in .fastq
  fns = [fn for fn in glob.glob("%s/*.fastq" % dirfn)]
  analyze(fns, verbose=True)
