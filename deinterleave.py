#! /bin/bash/env python
# Code to deinterleave FASTA files.
import sys
linenum = 0
infile = open(sys.argv[1], 'r')
with open('deinterleaved_' + sys.argv[1], 'w') as outfile:
    for line in infile:
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if linenum == 0:
                    outfile.write(line + '\n')
                    linenum += 1
                else:
                    outfile.write('\n' + line + '\n')
            else:
                outfile.write(line)