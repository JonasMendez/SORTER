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

# Jonas Mendez-Reneau1, Erin Sigel1
# University of Louisiana Lafayette

# Contact:jonasmrgrad@gmail.com

#BEFORE RUNNING MAKE SURE YOU HAVE DONE THE FOLLOWING:
#1. MAKE A FOLDER TO SERVE AS THE WORKING DIRECTORY FOR THE SET OF ALLOPOLYPLOIDS OR HYBRIDS TO BE PHASED AS /workingdirectory/phaseset/
#2. ADD RAW PAIR-END FASTQ READS FOR EACH SAMPLE IN PHASESET WORKING DIRECTORY, FOLLOWING THIS NAMING SCHEME FOR PAIRED END READS, RESPECTIVELY:
	# @@##_species1id_R1.fastq and @@##_species1id_R2.fastq
	#Make sure paired end reads are signified with R1/R2 as shown above
	#@@ values can be any set of two alphabetic letters and ## can be any unique set of two intergers per sample, MUST BE IDENTICAL FOR THE SAME SAMPLE
	#These Unique Identifiers are used to differentiate multiple samples of with the same species id (i.e. if multiple samples have the same 'species1id')
#4. PLACE ALL EXTERNAL PYTHON SCRIPTS IN WORKING DIRECTORY:
	#workingdirectory/getlongestcontig.py
	#workingdirectory/deinterleave.py
	#workingdirectory/seqclean.py
	#workingdirectory/annotatedupes
	#workingdirectory/keeplongest.py
#5. Know the working directory of your probe references for command line input, may be placed in /workingdirectory/
#6. Before running script make sure to load all required software and python version 

#FLAGS

#-wd WORKING DIRECTORY 
#-ref PROBE REFERENCE DIRECTORY (FASTA FILE) SAME AS PHASE1
#-c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .85 - .90 Recommended for allopolyploids/hybrids) 
#-trimgalore RUN TRIMGALORE TO TRIM RAW READS? (T/F)***
#-spades RUN SPADES ASSEMBLY? (T/F)***
#-onlyprocess (T/F) ONLY RUN TRIMGALORE AND SPADES FOR CONTIG PROCESSING; RUN AGAIN WITH -trim and -spades as F FOR PIPELINE (set as F if running processing + pipeline in one run)
#-cs TAKE SPADES CONTIGS or SCAFFOLDS? (MUST INPUT AS: scaffold or contig)
#-csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES. DEPENDING ON SUSPECTED PLOIDY, MULTIPLY BY 2 (i.e. tetraploid 4*2=8, triploid 3*2=6 ; if heterozygous variants are rertieved. 8-10 for unknown samples is recommended)
#-csl ONLY TAKE CONTIGS/SCAFFOLDS LARGER THAN N LENGTH
#-pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)
#-n PROPORTION OF MISSING DATA (N)
# ALLOWED IN PHASED SEQUENCES? (atleast 50% bp representation recommended, input as -n 50 , NOT AS DECIMAL)
#-al number of iterations for MAFFT alignments (1000 recommended)
#-indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)

#outputfiles
#_phase.fasta extention making polyploid/hybrid sample ids in the following format : >@@##_sampleid_phase(0/1)
#_diploidhit.fasta extention making polyploid/hybrid sample ids in the following format : >@@##_sampleid_diploidhit
#fully annotated as @@##_sampleid_diploidhit_phase(0/1)


#***Make 'F' if you have already RUN TRIMGALORE or SPADES for your paired end reads
#***If you have already run trimgalore, and want this script to run SPADES assembly, you must organize your working directory
#so that each set of trimmed reads has the naming scheme W@@##_species1_R1_val_1.fq and W@@##_species1_R2_val_2.fq
#Each set of trimmed .fq files is in their own folder corresponding to each sample, with each sample folder named as follows:
	#/workingdirectory/WA01_species1_R1.fastq
	#/workingdirectory/WA02_species2_R1.fastq
	#/workingdirectory/WB03_species3_R1.fastq
	#etc...
#folders have ..._R1.fastq extensions as unique identifiers for the script.


#TYPICAL COMMAND LINE
#python phase1.py -wd /workingdirectory/ -ref /workingdirectory/references.fasta  -spades T -trimgalore T -op F -c1 .85 -cs contig -csn 6 -csl 350 -pq 20 -n 50 -al 1000 -indel 0.25

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-c1", "--clust1id")
parser.add_argument("-spades", "--spadesassembly")
parser.add_argument("-trimgalore", "--trimgalore")
parser.add_argument("-op", "--onlyprocess")
parser.add_argument("-cs", "--contigscaf")
parser.add_argument("-csn", "--contigscafnum")
parser.add_argument("-csl", "--contigscaflen")
parser.add_argument("-ref", "--refdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-n", "--phasen")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
args = parser.parse_args()
phaseset=args.workingdir + 'phaseset/'
longestcontig=args.workingdir + 'getlongestcontig.py'
dintdir=args.workingdir + 'deinterleave.py'
seqclean=args.workingdir + 'seqclean.py'
annotatedupes = args.workingdir+'annotatedupes'
baitid1= ["L%d_" % x for x in range(455)]
baitid= ["L%d" % x for x in range(455)]
diploidclusters=args.workingdir + 'diploidclusters/'
keeplongest=args.workingdir + 'keeplongest.py'

# #Define command to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


#Trim Reads if Needed
if args.trimgalore is 'T':

	for file in os.listdir(phaseset):
		if 'R1' in file:
			readst= phaseset + file + '/'
			os.makedirs(readst)
			print(readst)
			old_path = os.path.join(args.workingdir, file)
			new_path = os.path.join(readst + file)
			os.renames(old_path, new_path)
			nextread = next(phaseset)
			read2path = args.workingdir + nextread
			newread2path = readst + nextread
			move(read2path, newread2path)
			os.chdir(readst)
			subprocess.call(["trim_galore --quality 20 --length 30 --paired --fastqc %s %s" % (file, nextread)], shell=True)
			os.remove(file)
			os.remove(nextread)
			os.chdir(dst)
		else:
			continue
else:

#Assemble Contigs if needed

	if args.spadesassembly is 'T':

	#Assemble allopolyploid or hybrid contigs with spades

		for file in os.listdir(phaseset):
			if 'R1' in file:
				print(file)
				os.chdir(phaseset + file)
				cdir=phaseset + file
				print(cdir)
				for read in os.listdir(cdir):
					if 'R1_val_1.fq' in read:
						R2 = read[:-11] + 'R2_val_2.fq'
						subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (read, R2)], shell=True)
						os.chdir(phaseset)
					else:
						continue
			else:
				continue

	elif args.onlyprocess is 'T':

			sys.exit("--onlyprocess = T, exiting script")
	else:

#Map Contigs to Rereferences
			os.chdir(phaseset)

			for folder in os.listdir(phaseset):
				if 'R1' in folder:
					os.chdir(phaseset + folder)
					dirpath =  phaseset + folder + "/"
					iterpath = os.listdir(dirpath)
					read = folder[:-6] + 'R1_val_1.fq'
					R2 = read[:-6] + 'R2_val_2.fq'
					for file in iterpath:
						if file.endswith("_hybrid_assembly"):
							os.chdir(phaseset + folder + "/" + file + "/")
							filepath =  phaseset + folder + "/" + file + "/"
							iterpath2 = os.listdir(filepath)
							for contig in iterpath2:
								if contig.endswith(args.contigscaf + "s.fasta"):
									conscontig= args.refdir
									subprocess.call(["bwa mem -V %s %s > %s_%smap.sam" % (conscontig, contig, folder, args.contigscaf)], shell=True)
									subprocess.call(["samtools view -S -F 4 *_%smap.sam | awk -v OFS='\t' '{print \"> \" $3 \"\\n \" $10}' > %smap.fa " % (args.contigscaf, folder + args.contigscaf)], shell=True)
									src = folder + args.contigscaf + 'map.fa'
									dst = dirpath + folder + '_' + args.contigscaf + 'map.fa'
									os.rename(src, dst)
									os.remove(folder + '_' + args.contigscaf + 'map.sam')
		os.chdir(phaseset)

		#Make fasta files for contigs which mapped to the same reference
		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				print(folder)
				os.chdir(phaseset + folder)
				readir = os.listdir(phaseset + folder)
				subprocess.call(["pwd"], shell=True)
				for file in readir:
					if file.endswith(args.contigscaf + "map.fa"):
		  				print(file)
		   				with open(file, 'r') as contigfile:
		   					for line in contigfile:
		   						for id in baitid1:
		   							if id in line:
		   								#print(line)
		   								with open(os.path.join(phaseset + folder, id), 'a') as idx:
												while True:
													try:
														idx.write(line)
														seq = next(contigfile)
														idx.write(seq)
														print(line)
														print seq
														break
													except StopIteration as e:
														print(e)
														break

		# #Rename contigs with locus ID
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				readir = os.listdir(phaseset + folder)
				subprocess.call(["pwd"], shell=True)
				for file in readir:
					if file.endswith("_"):
						with open(file, 'r') as infile:
		 					for line in infile:
		   						for id in baitid1:
		   							if id in line:
		   								print(line)
		   								name = line[:-100] + '>' + folder[:-9] + '_' + id + '\n'
		   								print(name)
		   								replaceAll(file, line, name)		
		os.chdir(phaseset)

		print('Take ' + args.contigscafnum + ' longest scaffolds for each locus per sample, then removing any scaffolds smaller than ' + args.contigscaflen + ' bp')

		#Take longest contig from contig set
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				readir = os.listdir(phaseset + folder)
				subprocess.call(["pwd"], shell=True)
				for file in readir:
					if file.endswith("_"):
		   				subprocess.call(["python %s -i %s -n %s > %slongest.fa" % (longestcontig, file, args.contigscafnum, file)], shell=True)	
		   				#Remove sequences shorter than 350bp
						subprocess.call(["seqtk seq -L %s %slongest.fa > %slongestfiltered.fa" % (args.contigscaflen, file, file)], shell=True )

		os.chdir(phaseset)

		#cluster highly similar contigs (i.e collapsing heterozygotes, or possibly homeologues depending on genetic distance)
		print('Clustering ' + args.contigscaf +'s ' + 'into Consensus Alleles at' + args.clust1id + ' Identity Threshold')

		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				print(folder)
				os.chdir(phaseset + folder)
				readir = os.listdir(phaseset + folder)
				subprocess.call(["pwd"], shell=True)
				for file in readir:
					if file.endswith("longestfiltered.fa"):
						subprocess.call(["usearch -cluster_fast %s -id %s -consout %s_cons.fa" % (file, args.clust1id, file[:-18])], shell=True)


		#annotate cluster contigs with sample ID and locus
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				readir = os.listdir(phaseset + folder)
				subprocess.call(["pwd"], shell=True)
				for file in readir:
					if file.endswith("_cons.fa"):
						with open(file, 'r') as infile:
		   					for line in infile:
		   						if '>' in line:
		   							print(line)
		   							name = line[:-100] + '>' + folder[:-9] + '_' + file[:-8] + '\n'
		   							print(name)
		   							replaceAll(file, line, name)		
		os.chdir(phaseset)

		#remove unused sequences
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				dirpath =  phaseset + folder + "/"
				iterpath = os.listdir(dirpath)
				for file in iterpath:
					if file.endswith("_"):
						print('deleting: ' + file)
						os.remove(dirpath + file)
					else:
						continue
			else:
				continue

		os.chdir(phaseset)

		#remove unused sequences
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				dirpath =  phaseset + folder + "/"
				iterpath = os.listdir(dirpath)
				for file in iterpath:
					if file.endswith("_longest.fa"):
						print('deleting: ' + file)
						os.remove(dirpath + file)
					else:
						continue
			else:
				continue

		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				dirpath =  phaseset + folder + "/"
				iterpath = os.listdir(dirpath)
				for file in iterpath:
					if file.endswith("_longestfiltered.fa"):
						print('deleting: ' + file)
						os.remove(dirpath + file)
					else:
						continue
			else:
				continue


		os.chdir(phaseset)

		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				subprocess.call(["cat *_cons.fa  > %s_allcontigs_allbaits.fasta" % (folder[:-9])], shell=True)
				print("deinterleaving " + folder[:-9] +'_allcontigs_allbaits.fasta')
				subprocess.call(["python %s %s" % (dintdir, folder[:-9] +'_allcontigs_allbaits.fasta')], shell=True)
				os.remove(folder[:-9] +'_allcontigs_allbaits.fasta')
				#annotate sequences with identical sequence names
				subprocess.call(["awk -f %s %s > %s_allcontigs_allbaits_annotated.fasta" % (annotatedupes, 'deinterleaved_' + folder[:-9] +'_allcontigs_allbaits.fasta', folder[:-9])], shell=True)
				os.remove('deinterleaved_' + folder[:-9] +'_allcontigs_allbaits.fasta')

		os.chdir(phaseset)


		#strip lines for final contigfile
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset + folder)
				for file in os.listdir(phaseset + folder):
					if file.endswith('_allcontigs_allbaits_annotated.fasta'):
						with open(file, 'r') as infile, \
							open(file[:-6] + '_final.fasta', 'w') as outfile:
							for line in infile:
								if line.strip():
									outfile.write(line)


		os.chdir(phaseset)

		# #Clean contig set, remove duplicate sequences, remove sequences smaller than 350bp

		for file in os.listdir(phaseset):
			if 'R1' in file:
				os.chdir(phaseset + file)
				path=phaseset + file
				for seqs in os.listdir(path):
					if seqs.endswith("annotated_final.fasta"):
						print(seqs)
						subprocess.call(["python %s %s %s %s" % (seqclean, seqs, args.contigscaflen, args.phasen)], shell=True)
						os.remove(seqs)
						os.remove(file[:-9] + '_allcontigs_allbaits_annotated.fasta' )

						#os.remove(file[:-9] + '_allcontigs_allbaits_annotated_final.fasta') # if you want to keep full contig dataset with dupes/seqs < 350bp do not use this line

		os.chdir(phaseset)

		# # #Map reads to locus-clustered contig-cluster references
		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				os.chdir(phaseset + folder)
				for baits in os.listdir(phaseset + folder):
					if baits.startswith('clear'):
						read = folder[:-8] + 'R1_val_1.fq'
						print(read)
						R2 = folder[:-8] + 'R2_val_2.fq'
						print(R2)
						print(baits)
						subprocess.call(["bwa index %s" % (baits)], shell=True)
						subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (baits, read, R2, folder[:-8])], shell=True)
						subprocess.call(["samtools sort %smapreads.sam -o %s" % (folder[:-8], folder[:-8] + 'mapreads.bam')], shell=True)
						subprocess.call(["samtools index  %s" % (folder[:-8] + 'mapreads.bam')], shell=True)
						subprocess.call(["samtools phase -A -Q %s -b %s %s" % (args.phasequal, folder[:-8], folder[:-8] + 'mapreads.bam')], shell=True)
						subprocess.call(["samtools sort  %s -o %s" % (folder[:-8] + '.0.bam', folder[:-8] + '.0srt.bam')], shell=True)
						subprocess.call(["samtools sort  %s -o %s" % (folder[:-8] + '.1.bam', folder[:-8] + '.1srt.bam')], shell=True)
						subprocess.call(["samtools sort  %s -o %s" % (folder[:-8] + '.chimera.bam', folder[:-8] + '.chimerasrt.bam')], shell=True)
						subprocess.call(["bcftools mpileup -Ov -d 500 -f %s %s | bcftools call -c -Ov -o %s " % (baits, folder[:-8] + '.0srt.bam', folder[:-8] + '_0.vcf' )], shell=True)
						subprocess.call(["bcftools mpileup -Ov -d 500 -f %s %s | bcftools call -c -Ov -o %s " % (baits, folder[:-8] + '.1srt.bam', folder[:-8] + '_1.vcf' )], shell=True)
						subprocess.call(["bcftools mpileup -Ov -d 500 -f %s %s | bcftools call -c -Ov -o %s " % (baits, folder[:-8] + '.chimerasrt.bam', folder[:-8] + '_chimera.vcf' )], shell=True)
						subprocess.call(["tabix %s" % (folder[:-8] + '_0.vcf')], shell=True)
						subprocess.call(["tabix %s" % (folder[:-8] + '_1.vcf')], shell=True)
						subprocess.call(["tabix %s" % (folder[:-8] + '_chimera.vcf')], shell=True)
						subprocess.call(["vcfutils.pl vcf2fq %s > %s_0.fastq" % (folder[:-8] + '_0.vcf', folder[:-8])], shell=True)
						subprocess.call(["vcfutils.pl vcf2fq %s > %s_1.fastq" % (folder[:-8] + '_1.vcf', folder[:-8])], shell=True)
						subprocess.call(["vcfutils.pl vcf2fq %s > %s_chimera.fastq" % (folder[:-8] + '_chimera.vcf', folder[:-8])], shell=True)
						subprocess.call(["seqtk seq -A %s > %s_0_Final.fasta" % (folder[:-8] + '_0.fastq', folder[:-8])], shell=True)
						subprocess.call(["seqtk seq -A %s > %s_1_Final.fasta" % (folder[:-8] + '_1.fastq', folder[:-8])], shell=True)
						subprocess.call(["seqtk seq -A %s > %s_chimera.fasta" % (folder[:-8] + '_chimera.fastq', folder[:-8])], shell=True)
						os.remove(folder[:-8] + 'mapreads.bam')
						os.remove(folder[:-8] + 'mapreads.sam')
						os.remove(folder[:-8] + '.0.bam')
						os.remove(folder[:-8] + '.0srt.bam')
						os.remove(folder[:-8] + '_0.vcf')
						os.remove(folder[:-8] + '_0.fastq')
						os.remove(folder[:-8] + '.1.bam')
						os.remove(folder[:-8] + '.1srt.bam')
						os.remove(folder[:-8] + '_1.vcf')
						os.remove(folder[:-8] + '_1.fastq')
						os.remove(folder[:-8] + '.chimera.bam')
						os.remove(folder[:-8] + '.chimerasrt.bam')
						os.remove(folder[:-8] + '_chimera.vcf')
						os.remove(folder[:-8] + '_chimera.fastq')

		os.chdir(phaseset)

# 		#concatenate cluster annotated baits for diploids and make database
 		os.chdir(diploidclusters)

		subprocess.call("cat *degap.fasta > ALLsamples_allcontigs_allbaits_clusterannotated.fasta", shell=True)
		subprocess.call("usearch -makeudb_usearch ALLsamples_allcontigs_allbaits_clusterannotated.fasta -output diploid_master.udb", shell=True)
		diploid_db = '/project/emsigel/jonasmr/rawreads/clean-fastq/diploidclusters/diploid_master.udb'

		os.chdir(phaseset)

		#Make a copy of diploid locus-clusters for each sample in phase set
		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				os.chdir(phaseset+folder)
				#os.makedirs('diploidclusters')
				dst= phaseset + folder + '/diploidclusters/'
				shutil.copytree(diploidclusters, dst)
				os.chdir(dst)
				for file in os.listdir(dst):
					if not 'deinterleaved' in file:
						os.remove(file)

		os.chdir(phaseset)


		#UBLAST Phased sequences for each sample to diploid locus-clusters to determine orthology
		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				os.chdir(phaseset+folder)
				path=phaseset+folder
				for file in os.listdir(path):
					if file.endswith('_Final.fasta'):
						subprocess.call(["awk 'BEGIN{FS=\" \"}{if(!/>/){print toupper($0)}else{print $1}}' %s > %s_cap.fasta" % (file, file[:-6])], shell=True)
						subprocess.call(["usearch -usearch_global %s -db %s -id 0.9 -top_hit_only -blast6out %s_hits.txt -strand plus" % (file[:-6]+'_cap.fasta', diploid_db, file[:-12])], shell=True)


		#Compile Polyploid into Diploid locus-cluster dataset, respectively

		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				os.chdir(phaseset+folder)
				for file in os.listdir(phaseset+folder):
					if file.endswith('0_hits.txt'):
						with open(file, 'r') as hits:
							for hitlines in hits:
								splithits=hitlines.split('	')
								splithits2=splithits[1].split('_')
										#print splithits
										# print splithits2
										# print splithits[0]
								with open(folder[:-8] + '_0_Final_cap.fasta', 'r') as phfinal:
									for line in phfinal:
										if splithits[0] in line:
											print('Match ' +splithits[0] + ' in ' + splithits[1])
											with open(phaseset + folder + '/' + 'diploidclusters/deinterleaved_' + splithits2[0] + '_' + splithits2[1] +'_degap.fasta', 'a+') as baitcluster:
												while True:
				 									try:
														splitline= line.split('_')
														print splitline
														name = '>' + splithits2[0] + '_' + splithits2[1] +'_' +splitline[0].split('>')[1] + '_'  + splitline[1] +  '_' + splithits2[3]#see if any statement for broader annotation
														print(name.rstrip('\n') + '_ph0' + '\n')
														baitcluster.write(name.rstrip('\n') + '_ph0' + '\n')
														seq = phfinal.next()
														print seq
														baitcluster.write(seq)
														break
													except StopIteration as e:
														print(e)
			 											break

		os.chdir(phaseset)

		for folder in os.listdir(phaseset):
			if 'R1' in folder:
				os.chdir(phaseset+folder)
				for file in os.listdir(phaseset+folder):
					if file.endswith('1_hits.txt'):
						with open(file, 'r') as hits:
							for hitlines in hits:
								splithits=hitlines.split('	')
								splithits2=splithits[1].split('_')
										#print splithits
										# print splithits2
										# print splithits[0]
								with open(folder[:-8] + '_1_Final_cap.fasta', 'r') as phfinal:
									for line in phfinal:
										if splithits[0] in line:
											print('Match splithits[0]= ' +splithits[0] + ' in ' + splithits[1])
											with open(phaseset + folder + '/' + 'diploidclusters/deinterleaved_' + splithits2[0] + '_' + splithits2[1]  + '_degap.fasta', 'a+') as baitcluster:
												while True:
				 									try:
														splitline= line.split('_')
														print splitline
														name = '>' + splithits2[0] + '_' + splithits2[1] +'_' +splitline[0].split('>')[1] + '_'  + splitline[1] +  '_' + splithits2[3]#see if any statement for broader annotation
														print(name.rstrip('\n') + '_ph1' + '\n')
														baitcluster.write(name.rstrip('\n') + '_ph1' + '\n')
														seq = phfinal.next()
														print seq
														baitcluster.write(seq)
														break
													except StopIteration as e:
														print(e)
			 											break

		os.chdir(phaseset)

		#Rewrite polyploid names as >WA01_polyploid_diploidhit_phase
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('degap.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								if sample in line:
									print(line)
									linspl=line.split('_')
									print(linspl)
									name = line[:-100] + '>' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '_' + linspl[5]
									print(name)
									replaceAll(file, line, name)

		os.chdir(phaseset)
		#Rename Diploid species as >WA01_diploid
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('degap.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									if not sample in line:
										print(line)
										linspl=line.split('_')
										print(linspl)
										name = name = line[:-100] + '>' + linspl[2] + '_' + linspl[3] + '\n'
										print(name)
										replaceAll(file, line, name)

		os.chdir(phaseset)

		# #Keep longest seq if two identical sequence ids are present (i.e. diploid/polyploid samples with multiple consensus alleles/cluster)
		print('Keeping longest seq if two identical sequence ids are present (i.e. diploid/polyploid samples with multiple sequences/cluster) ''\n' + 'These may represent heterozygous variants retrieved before phasing, discontinouos haplotype fragments or unclustered paralogs)')

		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				print('Processing ' + folder)
				os.chdir(phaseset+folder+'/diploidclusters/')
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('degap.fasta'):
						subprocess.call(["python %s %s > %s" % (keeplongest, file, file + '_dupremove.fa')], shell=True)
						os.remove(file)

		os.chdir(phaseset)

		#Align and trim locus-clusters
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('_dupremove.fa'):
						subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
						subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)
						os.remove(file)
						os.remove(file[:-6] + '_al.fasta')

		os.chdir(phaseset)

		#reannotate trimmed alignment to remove trimal bp annotation
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('trimmed.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									print(line)
									linspl=line.split(' ')
									name = line[:-100] + linspl[0] + '\n'
									print(name)
									replaceAll(file, line, name)

		#generate phase and diploid hit annotations files
		#Rewrite polyploid names as >WA01_polyploid_phase >WA01_polyploid_diploidhit_phase
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				print('Preparing Phase(0/1) Annotation alignments')
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('trimmed.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								with open(file[:-6]+'_phase.fasta',"a+") as phasefile:
									phasefile.write(line)

		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('phase.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								if sample in line:
									print(line)
									linspl=line.split('_')
									print(linspl)
									name = line[:-100] + linspl[0] + '_' + linspl[1] + '_' + linspl[3]
									print(name)
									replaceAll(file, line, name)

		#Rewrite polyploid names as >WA01_polyploid_diploidhit
		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				print('Preparing Diploid UBLAST Hit Annotation alignments')
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('trimmed.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								with open(file[:-6]+'_diploidhit.fasta',"a+") as phasefile:
									phasefile.write(line)

		for folder in os.listdir(phaseset):
			if 'fastq' in folder:
				os.chdir(phaseset+folder+'/diploidclusters/')
				sample=folder.split('_')[0] +'_' + folder.split('_')[1]
				for file in os.listdir(phaseset+folder+'/diploidclusters/'):
					if file.endswith('diploidhit.fasta'):
						with open(file, 'r') as infile:
							for line in infile:
								if sample in line:
									print(line)
									linspl=line.split('_')
									print(linspl)
									name = line[:-100] + linspl[0] + '_' + linspl[1] + '_' + linspl[2] + '\n' 
									print(name)
									replaceAll(file, line, name)



		sys.exit("A set of alignments with _phase.fasta and _diploidhit.fasta extensions were created with >WA01_sampleid_phase id format (_phase.fasta files) and >WA01_sampleid_diploidhit (_diploidhit.fasta)")
