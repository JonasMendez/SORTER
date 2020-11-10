#!/usr/bin/python

#https://gist.github.com/mkweskin/8869358

"""
Reads a FASTA file and if >1 sequence has the same description line,
it only keeps the longest sequence. It outputs all the sequencs to stdout
when complete.
"""


from Bio import SeqIO
import sys
if len(sys.argv) == 1:
  print "ERROR: Please enter the filename to read as the first argument after the program name"
  sys.exit()
else:
  file = sys.argv[1]

seqs = {}
new = 0  #For testing, count of sequences as their added to seqs
existing = 0  #For testing, count of sequences NOT added because they already exist

for seq_record in SeqIO.parse(file, "fasta"):
  if seq_record.name not in seqs:
    seqs[seq_record.name]=seq_record.seq
    new += 1
  else:
    existing += 1
    if len(seqs[seq_record.name])<=len(seq_record.seq):
       seqs[seq_record.name]=seq_record.seq
#       print "Duplicate found, ", seq_record.name

for name, seq in seqs.iteritems():
  print ">"+name
  print seq

#Uncomment below print statement to print out the sequence counts
#print new, existing, new+existing