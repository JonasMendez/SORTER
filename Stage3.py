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

#This script uses samples processed through stages 1 - 2 and uses them as potential progenitors for hybrid samples. It does this by mapping phased
#sequences in hybrids to progenitor sequences, and determining the highest similarity to progenitor sequences as a way to determine disparate haplotypes.
#trim hybrid and generate contigs for hybrid samples using Stage1A.py to generate sample folders which can then be placed in the "phaseset" folder
#in order to run this script and compare hybrid sequences to progenitor samples.
#The script generates consensus alleles and phases them as is done for diploids in Stages 1 - 2; thus you must re-select contig/scaffold parameters
#for generating consensus-alleles and phasing. For allopolyploids, go with a number higher than 2*ploidy of the sample, as you may retrieve alleles
#for each locus, thus for a tetraploid 2*4x = 8 possible haplotypes. For allopolyploids you may want to set -csn to 20 or more if you believe paralogs are 
#prevalent in the dataset.


# -wd WORKING DIRECTORY (DIRECTORY STRING STARTING AND ENDING IN '/'. e.g: "-wd /workingdirectory/" should be the cleanfastq directory created in Phase 1) 
# -pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)
# -al NUMBER OF ITERATIONS FOR MAFFT ALIGNMENTS (1000 recommended)
# -indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)
# -idformat (full/copies/onlysample/*) OUTPUTS FINAL ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:
# 	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. (Locus_ClusterID_@@##_sampleid_phase)
# 	-idformat phase = >@@##_sampleid_0 ; RECOMMENDED FOR PHASING. Keeps sample id and phase. i.e. if last annotation >1 signifies samples with multiple consensus alleles per locus-cluster that were phased.
# 	-idformat onlysample = >@@##_sampleid ; Keeps only the sample for sequence headers. Will have to decide how to manage phase or other sequence copies with identical id names.
# 	-idformat * :if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_@@##_sampleid_0 1230 bp
# -c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .98 - .99 Recommended, INPUT AS DECIMAL) 
# -cs TAKE SPADES CONTIGS or SCAFFOLDS? (input as: scaffold or contig)
# -csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES
# -csl ONLY TAKE CONTIGS LARGER THAN N LENGTH, UNLESS ALL CONTIGS/SCAFFOLDS FOR A SCAFFOLD ARE LESS THAN USER DEFINED LENGTH, THIS IS DONE TO IMPROVE COVERAGE FOR FULL TAXON REPRESENTATION DATASETS

# COMMAND LINE EXAMPLE
# python Stage3.py -wd /working/directory/ -ref /directory/to/reference.fasta -loci 450 -cs contig -csn 10 -csl 300 -c1 .99 -pq 20 -al 1000 -indel .25

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-c1", "--clust1id")
parser.add_argument("-loci", "--locinum")
parser.add_argument("-cs", "--contigscaf")
parser.add_argument("-csn", "--contigscafnum")
parser.add_argument("-csl", "--contigscaflen")
parser.add_argument("-ref", "--refdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
args = parser.parse_args()
phaseset=args.workingdir + 'phaseset/'
longestcontig=args.workingdir + 'getlongestcontig.py'
dintdir=args.workingdir + 'deinterleave.py'
seqclean=args.workingdir + 'seqclean.py'
annotatedupes = args.workingdir+'annotatedupes'
baitid1= ["L%d_" % x for x in range(int(args.locinum))]
baitid= ["L%d" % x for x in range(int(args.locinum))]
diploidclusters=args.workingdir + 'diploids_phased/'
keeplongest=args.workingdir + 'keeplongest.py'
diploid_db = args.workingdir + 'diploids_phased/diploid_master.udb'


# #Define command to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


print('SORTED Phase 3 will run with the following settings:' + '\n' + 'Working Directory = ' + args.workingdir + '\n' + 'PHASE CALLING DEPTH (SAMTOOLS PHASE -Q)= ' + args.phasequal + '\n' + 'MAFFT Alignment Iterations = ' + args.aliter + '\n' + 'Keep Indels when present in .X of Samples, X = ' + args.indelrep)


print("Removing any previous phase files...")

os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'clusters_phased' in folder:
		shutil.rmtree(phaseset + folder)


for folder in os.listdir(phaseset):
	if 'R1.fastq' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if not 'trimming_report' in file:
				if not 'spades' in file:
					if not '_val_' in file:
						os.remove(phaseset + folder + '/'+ file)
						os.chdir(phaseset)

os.chdir(diploidclusters)

#make ublast data
subprocess.call("usearch -makeudb_usearch ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta -output diploid_master.udb", shell=True)

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
						subprocess.call(["samtools view -S -F 4 *_%smap.sam | awk -v OFS='\t' '{print \">\" $3\"_\" \"\\n \" $10}' > %smap.fa " % (args.contigscaf, folder + '_' + args.contigscaf)], shell=True)
						src = folder +'_'+ args.contigscaf + 'map.fa'
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
   								name = '>' + folder[:-9] + '_' + id + '\n'
   								print(name)
   								replaceAll(file, line, name)		
os.chdir(phaseset)

print('Take ' + args.contigscafnum + ' longest '+args.contigscaf+'s for each locus per sample, then removing any '+args.contigscaf+'s smaller than ' + args.contigscaflen + ' bp')

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
print('Clustering ' + args.contigscaf +'s ' + 'into Consensus Homeologs at' + args.clust1id + ' Identity Threshold')

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
   							name = '>' + folder[:-9] + '_' + file[:-8] + '\n'
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

for folder in os.listdir(phaseset):
	if 'fastq' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith("_longestfiltered.fa"):
				print('deleting: ' + file)
				os.remove(dirpath + file)

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
				with open(file, 'r') as infile:
					with open(file[:-6] + '_final.fasta', 'w') as outfile:
						for line in infile:
							if line.strip():
								outfile.write(line)


os.chdir(phaseset)

# # #Clean contig set (minimum 150bp, >50% N, remove identical duplicate sequences)

for file in os.listdir(phaseset):
	if 'R1' in file:
		os.chdir(phaseset + file)
		path = phaseset + file
		for seqs in os.listdir(path):
			if seqs.endswith("_final.fasta"):
				print(seqs)
				subprocess.call(["python %s %s 150 50" % (seqclean, seqs)], shell=True)
				os.remove(seqs)	
				os.remove(file[:-9] + '_allcontigs_allbaits_annotated.fasta' )
				

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
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ov -d 500 -f %s %s | bcftools call -c --ploidy 1 -Ov -p .0001 -o %s " % (baits, folder[:-8] + '.0srt.bam', folder[:-8] + '_0.vcf' )], shell=True)
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ov -d 500 -f %s %s | bcftools call -c --ploidy 1 -Ov -p .0001 -o %s " % (baits, folder[:-8] + '.1srt.bam', folder[:-8] + '_1.vcf' )], shell=True)
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ov -d 500 -f %s %s | bcftools call -c --ploidy 1 -Ov -p .0001 -o %s " % (baits, folder[:-8] + '.chimerasrt.bam', folder[:-8] + '_chimera.vcf' )], shell=True)
				subprocess.call(["tabix %s" % (folder[:-8] + '_0.vcf')], shell=True)
				subprocess.call(["tabix %s" % (folder[:-8] + '_1.vcf')], shell=True)
				subprocess.call(["tabix %s" % (folder[:-8] + '_chimera.vcf')], shell=True)
				subprocess.call(["vcfutils.pl vcf2fq %s > %s_0.fastq" % (folder[:-8] + '_0.vcf', folder[:-8])], shell=True)
				subprocess.call(["vcfutils.pl vcf2fq %s > %s_1.fastq" % (folder[:-8] + '_1.vcf', folder[:-8])], shell=True)
				subprocess.call(["vcfutils.pl vcf2fq %s > %s_chimera.fastq" % (folder[:-8] + '_chimera.vcf', folder[:-8])], shell=True)
				subprocess.call(["seqtk seq -A %s > %s_0_Final.fasta" % (folder[:-8] + '_0.fastq', folder[:-8])], shell=True)
				subprocess.call(["seqtk seq -A %s > %s_1_Final.fasta" % (folder[:-8] + '_1.fastq', folder[:-8])], shell=True)
				subprocess.call(["seqtk seq -A %s > %s_chimera.fasta" % (folder[:-8] + '_chimera.fastq', folder[:-8])], shell=True)
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
				for bam in os.listdir(phaseset + folder):
					if bam.endswith("mapreads.bam"):
						statfilename = folder[:-8] + "readstats.txt"
						with open(os.path.join(phaseset+ folder, statfilename), 'a+') as statfile:
							statfile.write( folder[:-8] + " Read Statistics" + '\n')
							statfile.write("Proportion of Reads that Mapped to Locus-Cluster Reference" + "\n")
							subprocess.call(["samtools flagstat %s >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
							statfile.write("Mean Read Depth" + "\n")
							subprocess.call(["samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
							statfile.write("Breadth of Coverage" + "\n")
							subprocess.call(["samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
							statfile.close()
							os.remove(folder[:-8] + 'mapreads.bam')


os.chdir(phaseset)

#Make a copy of diploid locus-clusters for each sample in phase set
dpdst= phaseset + 'diploidclusters_phased/'
shutil.copytree(diploidclusters, dpdst)

for folder in os.listdir(phaseset):
	if 'diploidclusters' in folder:
		dipcl=phaseset+folder
		os.chdir(phaseset+folder)
		for file in os.listdir(dipcl):
			if 'deinterleaved' in file:
				os.remove(file)
			else:
				if 'phased' in file:
					os.remove(file)
				else:
					if 'master.udb' in file:
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
									with open(phaseset + 'diploidclusters_phased/'+ splithits2[0] + '_' + splithits2[1] +'__allsamples_allcontigs.fasta', 'a+') as baitcluster:
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
									with open(phaseset + 'diploidclusters_phased/'+ splithits2[0] + '_' + splithits2[1] +'__allsamples_allcontigs.fasta', 'a+') as baitcluster:
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
		os.chdir(phaseset+'diploidclusters_phased/')
		sample=folder.split('_')[0] +'_' + folder.split('_')[1]
		for file in os.listdir(phaseset+'diploidclusters_phased/'):
			if file.endswith('allsamples_allcontigs.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if sample in line:
							print(line)
							linspl=line.split('_')
							print(linspl)
							name = '>' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '_' + linspl[5]
							print(name)
							replaceAll(file, line, name)

	# os.chdir(phaseset)
	#Rename Diploid species as >WA01_diploid_phase  >L102_cl3_WF07_pellucidumacuminatum_ph0_0
	# for folder in os.listdir(phaseset):
	# 	if 'fastq' in folder:
	# 		os.chdir(phaseset+'diploidclusters_phased/')
	# 		sample=folder.split('_')[0] +'_' + folder.split('_')[1]
os.chdir(phaseset+'diploidclusters_phased/')
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('allsamples_allcontigs.fasta'):
		print(file)
		with open(file, 'r') as infile:
			for line in infile:
				for hyb in os.listdir(phaseset):
					if "diploidclusters" not in hyb:
						hybid = hyb.split('_')[0] +'_' + hyb.split('_')[1]
						if '>' in line:
							linspl1=line.split('_')
							# print(line)
							# print(len(linspl1))
							if len(linspl1) > 5:
								print(line)
								linspl=line.split('_')
								print(linspl)
								name = name '>' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '\n'
								print(name)
								replaceAll(file, line, name)

os.chdir(phaseset+'diploidclusters_phased/')

# #Keep longest seq if two identical sequence ids are present (i.e. diploid/polyploid samples with multiple consensus alleles/cluster; allopolyploids that did not differentiat homeologs)
print('Keeping longest seq if two identical sequence ids are present ''\n')

# for folder in os.listdir(phaseset):
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('allsamples_allcontigs.fasta'):
		subprocess.call(["python %s %s > %s" % (keeplongest, file, file[:-6] + '_dupremove.fa')], shell=True)
		# os.remove(file)

os.chdir(phaseset+'diploidclusters_phased/')

#Align and trim locus-clusters
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('_dupremove.fa'):
		subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
		subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)
		os.remove(file)
		os.remove(file[:-6] + '_al.fasta')

os.chdir(phaseset+'diploidclusters_phased/')

#reannotate trimmed alignment to remove trimal bp annotation
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('trimmed.fasta'):
		with open(file, 'r') as infile:
			for line in infile:
				if '>' in line:
					print(line)
					linspl=line.split(' ')
					name = linspl[0] + '\n'
					print(name)
					replaceAll(file, line, name)
