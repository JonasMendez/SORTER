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

# PHASE 2 OF THE SORTED PIPELINE TAKES THE LOCUS-CLUSTER DATASET GENERATED IN PHASE 1 AND PHASES CONSENSUS ALLELES BY REMAPPING TRIMMED READS TO SAMPLE-SPECIFIC LOCUS-CLUSTERS

# -wd WORKING DIRECTORY (DIRECTORY STRING STARTING AND ENDING IN '/'. e.g: "-wd /workingdirectory/" should be the same directory used in Phase 1) 
# -pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)
# -al NUMBER OF ITERATIONS FOR MAFFT ALIGNMENTS (1000 recommended)
# -indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)
# -idformat (full/copies/onlysample/*) OUTPUTS FINAL ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:
# 	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. (Locus_ClusterID_@@##_sampleid_phase)
# 	-idformat phase = >@@##_sampleid_0 ; RECOMMENDED FOR PHASING. Keeps sample id and phase. i.e. if last annotation >1 signifies samples with multiple consensus alleles per locus-cluster that were phased.
# 	-idformat onlysample = >@@##_sampleid ; Keeps only the sample for sequence headers. Will have to decide how to manage phase or other sequence copies with identical id names.
# 	-idformat * :if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_@@##_sampleid_0 1230 bp

# -cdbonly (T/F) set to 'T' to only make blast database for phasing and blasting allopolyploid/hyrbid samples in Phase3.py. 
# 	Use this if you wish to skip phasing your diploid locus-cluster samples but still want to process polyploids/hyrbids with Phase3.py

# COMMAND LINE EXAMPLE
# python Phase2.py -wd /working/directory/ -pq 20 -al 1000 -indel .25 -idformat phase -cdbonly F

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-idformat", "--idformat")
parser.add_argument("-cdbonly", "--clustdb")

args = parser.parse_args()
diploidclusters = args.workingdir + 'diploidclusters/'
contigdir=args.workingdir + 'diploids/'
longestcontig =args.workingdir + 'getlongestcontig.py'
phasedir = args.workingdir + 'diploids_phased/' #Make sure this directory has been made
dintdir=args.workingdir + 'deinterleave.py'
direc=os.listdir(args.workingdir)
map_contigs_to_baits_dir=sorted(os.listdir(contigdir))


#Define function to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

print('SORTED Phase 2 will run with the following settings:' + '\n' + 'Working Directory = ' + args.workingdir + '\n' + 'PHASE CALLING DEPTH (SAMTOOLS PHASE -Q)= ' + args.phasequal + '\n' + 'MAFFT Alignment Iterations = ' + args.aliter + '\n' + 'Keep Indels when present in .X of Samples, X = ' + args.indelrep + '\n' + 'Alignment Sequence ID Format= ' + args.idformat +'\n'+ 'Only Format Locus-Clusters?= ' + args.clustdb +'\n')


#Degap locus-clusters
os.chdir(diploidclusters)
print('degapping locus-clusters for phasing')
for file in os.listdir(diploidclusters):
	if file.endswith('_'):
		newfilename=file+ 'degap.fasta'
		print('degapping ' + file)
		with open(newfilename, "w") as o:
			for record in SeqIO.parse(file, "fasta"):
        			record.seq = record.seq.ungap("-")
        			SeqIO.write(record, o, "fasta")

os.chdir(diploidclusters)

#deinterlieave locus-clusters
print('deinterleaving locus-clusters for phasing')
for file in os.listdir(diploidclusters):
	if file.endswith('degap.fasta'):
		print("deinterleaving " + file)
		subprocess.call(["python %s %s" % (dintdir, file)], shell=True)
		os.remove(file)

os.chdir(diploidclusters)

#add new line in cluster files
for file in os.listdir(diploidclusters):
	if 'deinterleaved' in file:
		with open(file, 'a+') as cluster:
			cluster.write('\n')

os.chdir(diploidclusters)

subprocess.call("cat *degap.fasta > ALLsamples_allcontigs_allbaits_clusterannotated.fasta", shell=True)
subprocess.call("usearch -makeudb_usearch ALLsamples_allcontigs_allbaits_clusterannotated.fasta -output diploid_master.udb", shell=True)

if args.clustdb is 'T':
	sys.exit('-cdbonly = T ; Locus-Clusters ready for use in Phase3.py blast database. Exiting.')		
else:	

	# # ##Recompile sample specific clusterbaits into fastq folders for readmapping
	for file in os.listdir(diploidclusters):
		if 'deinterleaved' in file:
			with open(file, 'r') as allsamplefile:
				for line in allsamplefile:
					if '>' in line:
						linspl=line.split('_')
						sample=linspl[2] + '_' + linspl[3] + '_allcontigs_allclusterbaits_annotated.fasta'
						print(linspl)
						with open(os.path.join(args.workingdir + linspl[2] + '_' + linspl[3] + '_R1.fastq/', sample), 'a+') as idx:
							if not line.endswith('\n'):
								idx.write('\n')
							else:
								while True:
									try:
										idx.write(line)
										seq = next(allsamplefile)
										idx.write(seq)
										print(line)
										print seq
										break
									except StopIteration as e:
										print(e)
										break

	os.chdir(args.workingdir)

	#Map reads to consensus contig references
	for folder in direc:
		if 'R1' in folder:
			os.chdir(args.workingdir + folder)
			for baits in os.listdir(args.workingdir + folder):
				if baits.endswith('allcontigs_allclusterbaits_annotated.fasta'):
					read = folder[:-8] + 'R1_val_1.fq'
					print(read)
					R2 = folder[:-8] + 'R2_val_2.fq'
					print(R2)
					print(baits)
					subprocess.call(["bwa index %s" % (baits)], shell=True)
					subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (baits, read, R2, folder[:-8])], shell=True)
					subprocess.call(["samtools sort  %smapreads.sam -o %s" % (folder[:-8], folder[:-8] + 'mapreads.bam')], shell=True)
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

	os.chdir(args.workingdir)

	# annotate phased fasta files with alternative phase

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("0_Final.fasta"):
					with open(file, 'r') as infile:
	 					for line in infile:
	 						if '>' in line:
	 							print(line)
	 							name = line.rstrip('\n') + '_ph0' + '\n'
	 							print(name)
	 							replaceAll(file, line, name)
	os.chdir(args.workingdir)

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("1_Final.fasta"):
					with open(file, 'r') as infile:
	 					for line in infile:
	 						if '>' in line:
	 							print(line)
	 							name = line.rstrip('\n') + '_ph1' + '\n'
	 							print(name)
	 							replaceAll(file, line, name)
 	os.chdir(args.workingdir)


	##Concatenate Phased seqs files
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			subprocess.call(["cat *_Final.fasta  > %s_allcontigs_allclusterbaits_contigs_phased.fasta" % (folder[:-9])], shell=True)


	#Move Phased allcontigs_allbaits files to args.contigdir 'diploids_phased/' directory
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			for file in readir:
				if file.endswith('_allcontigs_allclusterbaits_contigs_phased.fasta'):
					src = args.workingdir + folder + '/' + file
					dst = phasedir + file
					os.rename(src, dst)

	os.chdir(diploidclusters)

	DICT2= {}

	for baitcluster in os.listdir(diploidclusters):
		if baitcluster.endswith('_'):
			if baitcluster not in DICT2:
				DICT2[baitcluster]={}

	#print(DICT2)

	for folder in map_contigs_to_baits_dir:
			if folder.endswith('_'):
	 			for baitcluster in DICT2.keys(): 
	 				DICT2[baitcluster][folder + 'ph0']=[]

	for folder in map_contigs_to_baits_dir:
			if folder.endswith('_'):
	 			for baitcluster in DICT2.keys(): 
	 				DICT2[baitcluster][folder + 'ph1']=[]

	#print(DICT2)

	os.chdir(phasedir)

	#concatenated phased sequences for all samples

	subprocess.call(["cat *_allcontigs_allclusterbaits_contigs_phased.fasta  > ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta"], shell=True)

# 	# # #filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				for record in input_fasta:
					bait= record.id.split('_', 1)[0]
					baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
					print(baitcluster)
					phase = (record.id.split('_', 5)[5]).rstrip('\n')
					folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3]+ '_' + phase
					print(folder)
					seq=record.seq
					DICT2[baitcluster][folder].append(seq)
					#print(DICT2)
					print baitcluster, folder, len(DICT2[baitcluster][folder])


	# #write fasta output summary files by bait
	for baitcluster in DICT2.keys():
		if len(DICT2[baitcluster])>0:
			outfile = open(baitcluster+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT2[baitcluster].keys():
				if len(DICT2[baitcluster][folder])>0:
					seq_list = DICT2[baitcluster][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						#print '>'+baitcluster+folder+index+'\n'+seq+'\n'
						outfile.write(str('>'+baitcluster+folder+'_'+index+'\n'+seq+'\n'))

	# ##output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	ph0 = ["{}{}".format(i,'ph0') for i in columns]
	ph1 = ["{}{}".format(i,'ph1') for i in columns]
	header = ['baitcluster']+ph0+ph1
	#print header
	with open('ALLsamples_allcontigs_allbaits_SUMMARY_TABLE_phased.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		samples = sorted(DICT2.values()[0].keys())
		#samples = columns[0:]
		for baitcluster in DICT2.keys():
			writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

 	os.chdir(phasedir)

	#Align and trim clusterbaits. Regions which don't share overlap (i.e. region with unique indel) in atleast 25 percent of the samples are removed. All clusters with single sequences not aligned and compiled downstream
	#Keep nontrimmed alignment; compare?

	for file in os.listdir(phasedir):
		if file.endswith('allsamples_allcontigs.fasta'):
			subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
			subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)

	for file in os.listdir(phasedir):
		if file.endswith('_al.fasta'):
			os.remove(file)

	# # # #deinterleave
	for file in os.listdir(phasedir):
		if 'trimmed' in file:
			print("deinterleaving " + file)
			subprocess.call(["python %s %s" % (dintdir, file)], shell=True)
			os.remove(file)

	os.chdir(phasedir)

	if 'full' in args.idformat:
		#Reformat as >L100_cl0_@@##_sampleid_phase_index ; Keeps full anottation.
		for file in os.listdir(phasedir):
			if 'deinterleaved' in file:
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = line[:-1000] + linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] +'_' +linspl2[5] +'\n'
							print(name)
							replaceAll(file, line, name)
		sys.exit('Kept full sequence ID annotations; e.g. >L100_cl0_WA10_sampleid_0')
	else:
		if 'phase' in args.idformat:
			#Reformat as >@@##_sampleid_0 ; RECOMMENDED FOR PHASED DATASETS. Keeps sample id and phase. i.e. if last annotation >1 signifies samples with multiple consensus alleles per locus-cluster which were phased.
			for file in os.listdir(phasedir):
				if 'deinterleaved' in file:
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl=line.split(' ')[0]
								linspl2=linspl.split('_')
								print(linspl2)
								name = line[:-1000] + '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] + '\n'
								print(name)
								replaceAll(file, line, name)
			sys.exit('Annotated alignments as: >@@##_sampleid_0/1 (annotated with phase)')
		else:
			if 'onlysample' in args.idformat:
				#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. Will have to decide how to manage alleles/sequence copies per sample with identical id names.
					if 'deinterleaved' in file:
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									print(line)
									linspl=line.split(' ')[0]
									linspl2=linspl.split('_')
									print(linspl2)
									name = line[:-1000] + '>' + linspl2[2] + '_' + linspl2[3] + '_' + '\n'
									print(name)
									replaceAll(file, line, name)
					sys.exit('Annotated alignments as: >@@##_sampleid (no phase annotations)')
			else:
				sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")
