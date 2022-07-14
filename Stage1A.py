import os
import sys
import subprocess
import shutil
from os import path
from shutil import move
import fileinput
import os
import csv 
import sys
import subprocess
import shutil
import Bio
from os import path
from shutil import move
from shutil import copyfile
from Bio import SeqIO
import re
import itertools 
import fileinput
import argparse
import glob

#This script will make a folder for each sample, trim raw reads (trimgalore) and build contigs (SPADES), preparing samples for input into SORTER
#Reads should be prelabeled to follow the following labeling scheme:
#  @@##_taxon_R1.fastq / @@##_taxon_R2.fastq , where @@## can be any unique number and letter combination (no underscores) to differentiate individuals.
#Folders for samples will follow the labeling scheme used, proper labeling of correct progenitor taxa vs hybrids is necessary for stage 3.

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
args = parser.parse_args()

os.chdir(args.workingdir)
dst = args.workingdir + 'cleanfastq/'
direc=os.listdir(args.workingdir)
directiter = iter(sorted(direc))

#trim reads, put paired reads in folders corresponding to each sample

for file in directiter:
	if 'R1' in file:
		readst= dst + file + '/'
		os.makedirs(os.path.join(readst))
		print(readst)
		old_path = os.path.join(args.workingdir, file)
		new_path = os.path.join(readst + file)
		os.renames(old_path, new_path)
		nextread = next(directiter)
		print(nextread)
		read2path = args.workingdir + nextread
		newread2path = readst + nextread
		move(read2path, newread2path)
		os.chdir(readst)
		subprocess.call(["trim_galore --quality 20 --length 30 --paired --fastqc %s %s" % (file, nextread)], shell=True)
		os.remove(file)
		os.remove(nextread)
		os.chdir(args.workingdir)
	else:
		continue

	#Assemble Contigs

for file in os.listdir(args.workingdir +'cleanfastq/'):
	if 'R1' in file:
		print(file)
		os.chdir(args.workingdir + 'cleanfastq/'+ file)
		cdir=args.workingdir + 'cleanfastq/' + file
		print(cdir)
		for read in os.listdir(cdir):
			if 'R1_val_1.fq' in read:
				R2 = read[:-11] + 'R2_val_2.fq'
				subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (read, R2)], shell=True)
				os.chdir(args.workingdir +'cleanfastq/')
			else:
				continue
	else:
		continue


sys.exit("Trimgalore and SPADES processing has finished, exiting script")
