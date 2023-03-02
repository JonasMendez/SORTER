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

# 1. WITHIN WORKING DIRECTORY (i.e. the 'cleanfastq' folder) MAKE THE FOLLOWING FOLDERS:
# 	workingdirectory/diploids (this folder outputs Phase1 alignments)
# 	workingdirectory/diploidclusters (This folder serves as a reference for Phase 3 allopolyploid/hybrid BLAST steps)
# 	workingdirectory/diploids_phased (if phasing with Phase2.py, this folder outputs Phase2 Phased alignments for diploids. Can be used for Phase 3. 
# 		***YOU WILL STILL HAVE TO RUN Phase2.phy with -cdbonly T if you wish to use Phase3.py to phase allopolyploids/hybrids without phasing in Phase2.phy)
# 	workingdirectory/phaseset (if using Phase3.py, this folder will contain hybrid/allopolyploid sample files and alignments)
# 2. PLACE ALL EXTERNAL PYTHON SCRIPTS IN WORKING DIRECTORY:
# 	workingdirectory/getlongestcontig.py
# 	workingdirectory/deinterleave.py
# 	workingdirectory/seqclean.py
# 	workingdirectory/annotatedupes
# 	workingdirectory/keeplongest.py
# 3. Know the working directory of your probe references for command line input, may be placed in /workingdirectory/
# 		REFERENCE FILES MAY HAVE ALTERNATIVE FORMATING, BUT TARGET GENES MUST BE THE BEGINNING OF THE STRING (after '>') FORMATTED AS L#_ WITH L1_ - LN_ (N=number of target reference loci) 
# 		AS FOLLOWS:
# 		>L1_*
# 		>L2_*
# 		>L3_*
# 		...
# 		>L20_*
# 		>L21_*
# 		...
# 		>L200_*
# 		etc... (*'s representing any other ID format present in your reference; The pipeline ignores this information and is just interested in identifying genes from references)
# 4. Before running script make sure to load all required software and python version



# FLAGS

# -wd WORKING DIRECTORY 
# -ref PROBE REFERENCE DIRECTORY (FASTA FILE)
# -loci NUMBER OF REFERENCE LOCI
# -c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .97 - .99 Recommended, INPUT AS DECIMAL) 
# -c2 2nd CLUSTER ID (AMONG SAMPLES FOR LOCUS-CLUSTERS, .65 - .80 Recommended***, INPUT AS DECIMAL. ***depends on taxonomic breadth of ingroup/outgroups;sensitivity analysis at different thresholds should be anlayzed for optimal clustering) 

# -reclust (T/F) IF T ONLY RERUNs LOCUS CLUSTERING (SKIPS SPADES/TRIMGALORE AND SAMPLE ANNOTATIONS BY DEFAULT). CLEARS PREVIOUS LOCUS-CLUSTER OUTPUT, BACK UP IF NEEDED TO COMPARE CLUSTERING THRESHOLDS.
# 	RUNNING THIS STEP IGNORES -c1, -cs -csl, -csn VALUES, USING VALUES FROM PREVIOUS RUN; TO CHANGE THESE FLAGS YOU MUST RERUN Phase1.py WITH -reclust F  spades F -trimgalore F -op F
# -cs TAKE SPADES CONTIGS or SCAFFOLDS? (input as: scaffold or contig)
# -csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES
# -csl ONLY TAKE CONTIGS LARGER THAN N LENGTH, UNLESS ALL CONTIGS/SCAFFOLDS FOR A SCAFFOLD ARE LESS THAN USER DEFINED LENGTH, THIS IS DONE TO IMPROVE COVERAGE FOR FULL TAXON REPRESENTATION DATASETS
# -al number of iterations for MAFFT alignments (1000 recommended)
# -indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)
# -idformat (full/copies/onlysample/*) OUTPUTS FINAL ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:
# 	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. If last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes, discontinouos haplotype fragments or unclustered paralogs, may consider higher-c2 (among samples) or -c1 (within samples) clustering values
# 	-idformat copies = >@@##_sampleid_0 ; Keeps sample id and consensus allele copy count per cluster. i.e. if last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; see above
# 	-idformat onlysample = >@@##_sampleid ; Keeps only the sample id across locus-cluster alignments, easiest for concatenation across all locus-cluster alignments
# 	-format * if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp




parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-c1", "--clust1id")
parser.add_argument("-c2", "--clust2id")
parser.add_argument("-reclust", "--recluster")
parser.add_argument("-loci", "--locinum")
parser.add_argument("-cs", "--contigscaf")
parser.add_argument("-csn", "--contigscafnum")
parser.add_argument("-csl", "--contigscaflen")
parser.add_argument("-ref", "--refdir")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-idformat", "--idformat")
args = parser.parse_args()
diploidclusters=args.workingdir + 'diploidclusters/'
longestcontig=args.workingdir + 'getlongestcontig.py'
dintdir=args.workingdir + 'deinterleave.py'
seqclean=args.workingdir + 'seqclean.py'
baitid1= ["L%d_" % x for x in range(int(args.locinum))]
baitid= ["L%d" % x for x in range(int(args.locinum))]
direc=os.listdir(args.workingdir)
contigdir=args.workingdir + 'diploids/'

#Define function to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


print('SORTED Phase 1 will run with the following settings:'+ '\n' + 'Working Directory= ' + args.workingdir + '\n' + 'Only doing reclustering step?=' + args.recluster + '\n' +'Use Spades contigs or scaffolds?=' + args.contigscaf + '\n' + 'Max Number of contigs/scaffolds per locus to use for consensus alleles = ' + args.contigscafnum + '\n' + ' Max length of contig/scaffolds retrieved = ' + args.contigscaflen + '\n' + 'Consensus Alleles (Within Samples) Clustering ID = ' + args.clust1id + '\n' +'Locus-Cluster (Among Samples) Clustering ID = ' + args.clust2id + '\n' + 'MAFFT Alignment Iterations = ' + args.aliter + '\n' + 'Keep Indels when present in .X of Samples, X = ' + args.indelrep )


if args.recluster is 'T':
	os.chdir(args.workingdir)
	#clear data from last clusetering run
	print('Preparing Locus-Cluster Directories for reclustering (i.e. using contigs from previous run)')
	for folder in os.listdir(args.workingdir):
		if folder.endswith(".fastq"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not file.endswith('map.fa'):
							if not file.endswith('cons.fa'):
								os.remove(dirpath + file)

	os.chdir(args.workingdir)

	for folder in os.listdir(contigdir):
		os.chdir(contigdir)
		os.remove(args.workingdir + 'diploids/'+ folder)

	os.chdir(args.workingdir)


	for folder in os.listdir(diploidclusters):
		os.chdir(diploidclusters)
		os.remove(args.workingdir + 'diploidclusters/'+ folder)

	os.chdir(args.workingdir)

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			subprocess.call(["cat *_cons.fa  > %s_allcontigs_allbaits_contigs.fasta" % (folder[:-9])], shell=True)

	os.chdir(args.workingdir)

	#Move allcontigs_allbaits files to contigdir 'diploids/' directory
	for folder in direc:
		if 'fastq' in folder:
			print('Moving ' + folder + ' contigs to /diploids/ directory')
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			for file in readir:
				if file.endswith('_allcontigs_allbaits_contigs.fasta'):
					src = args.workingdir + folder + '/' + file
					dst = contigdir + file[:-33]
					os.rename(src, dst)

	os.chdir(contigdir)

	#Make summary files of Consensus Alleles per Sample
	map_contigs_to_baits_dir=sorted(os.listdir(contigdir))

	for folder in map_contigs_to_baits_dir:
		print(folder)

	os.chdir(contigdir)

	DICT= {}
	for bait in baitid: 
		if bait not in DICT:
			DICT[bait]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				DICT[bait][folder]=[]


	#print(DICT)

	subprocess.call("cat *_ > ALLsamples_allcontigs_allbaits_contigs.fasta", shell=True)

	#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaits_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				for record in input_fasta:
					bait=record.id.split('_', 3)[2]
					print(bait)
					folder = record.id.split('_', 3)[0]+"_"+record.id.split('_', 3)[1]+"_"
					print(folder)
					seq=record.seq
					DICT[bait][folder].append(seq)
					#print(DICT)
					print bait, folder, len(DICT[bait][folder])
		#print DICT

	#write fasta output summary files by bait
	for bait in DICT.keys():
		if len(DICT[bait])>0:
			outfile = open(bait+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT[bait].keys():
				if len(DICT[bait][folder])>0:
					seq_list = DICT[bait][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						print '>'+bait+'_'+folder+'_'+index+'\n'+seq+'\n'
						outfile.write(str('>'+bait+'_'+folder+index+'\n'+seq+'\n'))


	#output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 

	columns =  [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['bait']+columns
	#print header
	with open('ALLsamples_consensusallele_c1'+args.clust1id+'_'+args.contigscaf+'_csl'+args.contigscaflen+'_csn'+args.contigscafnum+'_'+'SUMMARY_TABLE.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		samples = sorted(DICT.values()[0].keys())
		#samples = columns[0:]
		for bait in DICT.keys():
			writer.writerow([bait]+[len(DICT[bait][sample]) for sample in samples])


	os.chdir(contigdir)

	#Cluster contigs into orthologous sets and annotate
	print('Clustering Consensus Alleles Among Samples into Orthologous Locus-Clusters at ' + args.clust2id + ' Identity Threshold')
	for file in os.listdir(contigdir):
		if file.endswith('_allsamples_allcontigs.fasta'):
			sp=file.split('_')
			subprocess.call(["usearch -cluster_fast %s -sort length -id %s -msaout %s" % (file, args.clust2id ,sp[0] + '_cl')], shell=True)

	#Add _ to end of locus cluster files for processing
	for file in os.listdir(contigdir):
		if '_cl' in file:
			newfilename=file+'_'
			print(newfilename)
			os.rename(file,newfilename)

	os.chdir(contigdir)

	#Define Cluster IDs
	clustid= ["cl%d_" % x for x in range(1000)]

	# Annotate sequences with locus-cluster id:
	for file in os.listdir(contigdir):
		for id in baitid1:
			if id in file:
				for cid in clustid:
					if cid in file:
						with open(file, 'r') as infile:
	   						for line in infile:
	   							if '>' in line:
	   								print(line)
	   								linspl=line.split('_')
	   								print(linspl)
	   								name = linspl[0] + '_' + cid + linspl[1] + '_' + linspl[2] + '_'+ linspl[3]
	   								print(name)
	   								replaceAll(file, line, name)

	os.chdir(contigdir)


	##Align and trim regions which don't share overlap (i.e. region with unique indel, used to remove long flanking tails, can set -indel to 0.01 to keep all indels)
	#All clusters with single sequences not aligned and not compiled downstream
	for file in os.listdir(contigdir):
		if '_cl' in file:
			subprocess.call(["mafft --globalpair --maxiterate %s %s > %sal.fasta" % (args.aliter, file, file)], shell=True)
			subprocess.call(["trimal -in %s -out %s -gt %s" % (file + "al.fasta", file + 'trimmed', args.indelrep)], shell=True)
			os.remove(file +"al.fasta")

	#deinterleave trimal output
	os.chdir(contigdir)

	for file in os.listdir(contigdir):
		if 'trimmed' in file:
			print("deinterleaving " + file)
			subprocess.call(["python %s %s" % (dintdir, file)], shell=True)
			os.remove(file)

	# Move raw ualigned/untrimmed locus-clusters into diploidclusters folder.
	for file in os.listdir(contigdir):
		if '_cl' in file:
			if not 'trimmed' in file:
				src =contigdir + file
				print(src)
				dst = args.workingdir + 'diploidclusters/'+ file
				print(dst)
				os.rename(src, dst)
				print(file)

	os.chdir(args.workingdir)

	os.chdir(diploidclusters)

	#Generate summary .csv's and fastas of locus clusters
	
	DICT2= {}

	for baitcluster in os.listdir(diploidclusters):
		if baitcluster.endswith('_'):
			if baitcluster not in DICT2:
				DICT2[baitcluster]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				DICT2[baitcluster][folder]=[]


	os.chdir(diploidclusters)

	subprocess.call(["cat *_  > ALLsamples_allcontigs_allbaitclusters_contigs.fasta"], shell=True)
	# #filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				for record in input_fasta:
					bait= record.id.split('_', 1)[0]
					baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
					print(baitcluster)
					folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_'
					print(folder)
					seq=record.seq
					DICT2[baitcluster][folder].append(seq)
					#print(DICT2)
					print baitcluster, folder, len(DICT2[baitcluster][folder])
		#print DICT2

	# #write fasta output summary files by baitcluster
	for baitcluster in DICT2.keys():
		if len(DICT2[baitcluster])>0:
			outfile = open(baitcluster+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT2[baitcluster].keys():
				if len(DICT2[baitcluster][folder])>0:
					seq_list = DICT2[baitcluster][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						print '>'+baitcluster+folder+index+'\n'+seq+'\n'
						outfile.write(str('>'+baitcluster+folder+index+'\n'+seq+'\n'))


	# ##output the nested dictory to a csv file that can be exported into an excel table where rows are baitclusters, columns are samples, and cell values are number of contigs 
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['baitcluster']+columns
	#print header
	with open('ALLsamples_allclusterbaits_cid'+ args.clust2id[1:] + '_SUMMARY_TABLE.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		samples = sorted(DICT2.values()[0].keys())
		#samples = columns[0:]
		for baitcluster in DICT2.keys():
			writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

	os.chdir(contigdir)

	if 'full' in args.idformat:
		#Reformat as >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. If last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes
		for file in os.listdir(contigdir):
			if file.endswith('trimmed'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] + '\n'
							print(name)
							replaceAll(file, line, name)
	else:
		if 'copies' in args.idformat:
			#Reformat as >@@##_sampleid_0 ; Keeps sample id and consensus allele copy count per cluster. i.e. if last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes or discontinuous orthologous fragments
			for file in os.listdir(contigdir):
				if file.endswith('trimmed'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl=line.split(' ')[0]
								linspl2=linspl.split('_')
								print(linspl2)
								name = '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] + '\n'
								print(name)
								replaceAll(file, line, name)
		else:
			if 'onlysample' in args.idformat:
				#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. May have to deal with multiple allele copies per sample with identical id names.
				for file in os.listdir(contigdir):
					if file.endswith('trimmed'):
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									print(line)
									linspl=line.split(' ')[0]
									linspl2=linspl.split('_')
									print(linspl2)
									name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
									print(name)
									replaceAll(file, line, name)
			else:
				sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")

else:

	print('Preparing Directories For contig assembly and ortholog clustering...')
	for folder in os.listdir(args.workingdir):
		if folder.endswith(".fastq"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not 'trimming_report' in file:
							os.remove(dirpath + file)
						else:
							continue

	os.chdir(args.workingdir)

	#make bwa index of consensusbaits
	subprocess.call(["bwa index %s" % (args.refdir)], shell=True)						

	#map contigs to consensus baits, move contigs as scaffoldmap.fa or contigmap.fa to root sample directory (WA01_species)

	for folder in direc:
		if 'R1' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_hybrid_assembly"):
					os.chdir(args.workingdir + folder + "/" + file + "/")
					filepath =  args.workingdir + folder + "/" + file + "/"
					iterpath2 = os.listdir(filepath)
					for contig in iterpath2:
						if contig.endswith(args.contigscaf + "s.fasta"):
							print(folder)
							print(file)
							conscontig= args.refdir
							subprocess.call(["bwa mem -V %s %s > %s_%smap.sam" % (conscontig, contig, folder, args.contigscaf)], shell=True)
							subprocess.call(["samtools view -S -F 4 *_%smap.sam | awk -v OFS='\t' '{print \">\" $3\"_\" \"\\n\" $10}' > %smap.fa " % (args.contigscaf, folder + '_' + args.contigscaf)], shell=True)
							src = folder + '_' + args.contigscaf + 'map.fa'
							dst = dirpath + folder + '_' + args.contigscaf + 'map.fa'
							os.rename(src, dst)
							os.remove(folder + '_' + args.contigscaf + 'map.sam')
	os.chdir(args.workingdir)

	#make files for all contigs or scaffolds corresponding to each locus per sample

	for folder in direc:
		if 'R1' in folder:
			print(folder)
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith(args.contigscaf + "map.fa"):
	  				print(file)
	   				with open(file, 'r') as contigfile:
	   					for line in contigfile:
	   						for id in baitid1:
	   							if id in line:
	   								print(line)
	   								with open(os.path.join(args.workingdir + folder, id), 'a') as idx:
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




	#annotate contigs with locus information

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
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
	os.chdir(args.workingdir)

	print('Take ' + args.contigscafnum + ' longest '+ args.contigscaf +'s for each locus per sample, then removing any '+ args.contigscaf +'s smaller than ' + args.contigscaflen + ' bp')

	#Take longest contig from contig set
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("_"):
					subprocess.call(["python %s -i %s -n %s > %slongest.fa" % (longestcontig, file, args.contigscafnum, file)], shell=True)	
					#Remove sequences shorter than user defined length
					subprocess.call(["seqtk seq -L %s %slongest.fa > %slongestfiltered.fa" % (args.contigscaflen, file, file)], shell=True )

	os.chdir(args.workingdir)

	print('Clustering ' + args.contigscaf +'s ' + 'into Consensus Alleles at' + args.clust1id + ' Identity Threshold')

	#cluster contigs

	for folder in direc:
		if 'R1' in folder:
			print(folder)
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("longestfiltered.fa"):
					subprocess.call(["usearch -cluster_fast %s -id %s -consout %s_cons.fa" % (file, args.clust1id, file[:-18])], shell=True)


	# #annotate contig-consensus fastas with sample ID and locus
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
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
	os.chdir(args.workingdir)

	# #remove unclustered and excess Locus Sequence Files
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_"):
					print('deleting: ' + file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue

	os.chdir(args.workingdir)

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_longest.fa"):
					print('deleting: ' + file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue


	os.chdir(args.workingdir)

	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_longestfiltered.fa"):
					print(file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue

	os.chdir(args.workingdir)
	#clear data from last clusetering run
	print('Preparing Locus-Cluster Directories')
	for folder in os.listdir(args.workingdir):
		if folder.endswith(".fastq"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not file.endswith('map.fa'):
							if not file.endswith('cons.fa'):
								os.remove(dirpath + file)

	os.chdir(args.workingdir)

	for folder in os.listdir(contigdir):
		os.chdir(contigdir)
		os.remove(args.workingdir + 'diploids/'+ folder)

	os.chdir(args.workingdir)


	for folder in os.listdir(diploidclusters):
		os.chdir(diploidclusters)
		os.remove(args.workingdir + 'diploidclusters/'+ folder)

	os.chdir(args.workingdir)

	#make master file for all contigs
	for folder in direc:
		if 'fastq' in folder:
			os.chdir(args.workingdir + folder)
			subprocess.call(["cat *_cons.fa  > %s_allcontigs_allbaits_contigs.fasta" % (folder[:-9])], shell=True)

	os.chdir(args.workingdir)

	#Move allcontigs_allbaits files to contigdir 'diploids/' directory
	for folder in direc:
		if 'fastq' in folder:
			print('Moving ' + folder + ' contigs to /diploids/ directory')
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			for file in readir:
				if file.endswith('_allcontigs_allbaits_contigs.fasta'):
					src = args.workingdir + folder + '/' + file
					dst = contigdir + file[:-33]
					os.rename(src, dst)

	os.chdir(contigdir)

	# #Make summary files of Consensus Alleles per Sample
	map_contigs_to_baits_dir=sorted(os.listdir(contigdir))

	for folder in map_contigs_to_baits_dir:
		print(folder)

	os.chdir(contigdir)

	DICT= {}
	for bait in baitid: 
		if bait not in DICT:
			DICT[bait]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				DICT[bait][folder]=[]


	print(DICT)

	subprocess.call("cat *_ > ALLsamples_allcontigs_allbaits_contigs.fasta", shell=True)

	#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaits_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 		
				for record in input_fasta:
					bait=record.id.split('_', 3)[2]
					print(bait)
					folder=record.id.split('_', 3)[0]+"_"+record.id.split('_', 3)[1]+"_"
					print(folder)
					seq=record.seq
					DICT[bait][folder].append(seq)
					#print(DICT)
					print bait, folder, len(DICT[bait][folder])
		#print DICT

	#write fasta output summary files by bait
	for bait in DICT.keys():
		if len(DICT[bait])>0:
			outfile = open(bait+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT[bait].keys():
				if len(DICT[bait][folder])>0:
					seq_list = DICT[bait][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						print '>'+bait+'_'+folder+'_'+index+'\n'+seq+'\n'
						outfile.write(str('>'+bait+'_'+folder+index+'\n'+seq+'\n'))


	#output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 

	columns =  [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['bait']+columns
	#print header
	with open('ALLsamples_consensusallele_c1'+args.clust1id+'_'+args.contigscaf+'_csl'+args.contigscaflen+'_csn'+args.contigscafnum+'_'+'SUMMARY_TABLE.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		samples = sorted(DICT.values()[0].keys())
		#samples = columns[0:]
		for bait in DICT.keys():
			writer.writerow([bait]+[len(DICT[bait][sample]) for sample in samples])


	os.chdir(contigdir)

	#Cluster contigs into orthologous sets and annotate
	print('Clustering Consensus Alleles Among Samples into Orthologous Locus-Clusters at ' + args.clust2id + ' Identity Threshold')
	for file in os.listdir(contigdir):
		if file.endswith('_allsamples_allcontigs.fasta'):
			sp=file.split('_')
			subprocess.call(["usearch -cluster_fast %s -sort length -id %s -msaout %s" % (file, args.clust2id ,sp[0] + '_cl')], shell=True)

	#Add _ to end of locus cluster files for processing
	for file in os.listdir(contigdir):
		if '_cl' in file:
			newfilename=file+'_'
			print(newfilename)
			os.rename(file,newfilename)

	os.chdir(contigdir)

	#Define Cluster IDs
	clustid= ["cl%d_" % x for x in range(1000)]

	# # Annotate sequences with locus-cluster id:
	for file in os.listdir(contigdir):
		for id in baitid1:
			if id in file:
				for cid in clustid:
					if cid in file:
						with open(file, 'r') as infile:
	   						for line in infile:
	   							if '>' in line:
	   								print(line)
	   								linspl=line.split('_')
	   								print(linspl)
	   								name = linspl[0] + '_' + cid + linspl[1] + '_' + linspl[2] + '_'+ linspl[3]
	   								print(name)
	   								replaceAll(file, line, name)

	os.chdir(contigdir)


	##Align and trim regions which don't share overlap (i.e. region with unique indel, used to remove long flanking tails, can set -indel to 0.01 to keep all indels)
	#All clusters with single sequences not aligned and not compiled downstream
	for file in os.listdir(contigdir):
		if '_cl' in file:
			subprocess.call(["mafft --globalpair --maxiterate %s %s > %sal.fasta" % (args.aliter, file, file)], shell=True)
			subprocess.call(["trimal -in %s -out %s -gt %s" % (file + "al.fasta", file + 'trimmed', args.indelrep)], shell=True)
			os.remove(file +"al.fasta")

	#deinterleave trimal output
	os.chdir(contigdir)

	for file in os.listdir(contigdir):
		if 'trimmed' in file:
			print("deinterleaving " + file)
			subprocess.call(["python %s %s" % (dintdir, file)], shell=True)
			os.remove(file)

	# Move raw ualigned/untrimmed locus-clusters into diploidclusters folder.
	for file in os.listdir(contigdir):
		if '_cl' in file:
			if not 'trimmed' in file:
				src =contigdir + file
				print(src)
				dst = args.workingdir + 'diploidclusters/'+ file
				print(dst)
				os.rename(src, dst)
				print(file)

	os.chdir(args.workingdir)

	os.chdir(diploidclusters)

	#Generate summary .csv's and fastas of locus clusters

	DICT2= {}

	for baitcluster in os.listdir(diploidclusters):
		if baitcluster.endswith('_'):
			if baitcluster not in DICT2:
				DICT2[baitcluster]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				DICT2[baitcluster][folder]=[]


	os.chdir(diploidclusters)

	subprocess.call(["cat *_  > ALLsamples_allcontigs_allbaitclusters_contigs.fasta"], shell=True)
	# #filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				for record in input_fasta:
					bait= record.id.split('_', 1)[0]
					baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
					print(baitcluster)
					folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_'
					print(folder)
					seq=record.seq
					DICT2[baitcluster][folder].append(seq)
					#print(DICT2)
					print baitcluster, folder, len(DICT2[baitcluster][folder])
		#print DICT2

	# #write fasta output summary files by baitcluster
	for baitcluster in DICT2.keys():
		if len(DICT2[baitcluster])>0:
			outfile = open(baitcluster+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT2[baitcluster].keys():
				if len(DICT2[baitcluster][folder])>0:
					seq_list = DICT2[baitcluster][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						print '>'+baitcluster+folder+index+'\n'+seq+'\n'
						outfile.write(str('>'+baitcluster+folder+index+'\n'+seq+'\n'))


	# ##output the nested dictory to a csv file that can be exported into an excel table where rows are baitclusters, columns are samples, and cell values are number of contigs 
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['baitcluster']+columns
	#print header
	with open('ALLsamples_allclusterbaits_cid'+ args.clust2id[1:] + '_SUMMARY_TABLE.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		samples = sorted(DICT2.values()[0].keys())
		#samples = columns[0:]
		for baitcluster in DICT2.keys():
			writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

	os.chdir(contigdir)

	if 'full' in args.idformat:
		#Reformat as >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. If last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes
		for file in os.listdir(contigdir):
			if file.endswith('trimmed'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] + '\n'
							print(name)
							replaceAll(file, line, name)
	else:
		if 'copies' in args.idformat:
			#Reformat as >@@##_sampleid_0 ; Keeps sample id and consensus allele copy count per cluster. i.e. if last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes or discontinuous orthologous fragments
			for file in os.listdir(contigdir):
				if file.endswith('trimmed'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl=line.split(' ')[0]
								linspl2=linspl.split('_')
								print(linspl2)
								name =  '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4] + '\n'
								print(name)
								replaceAll(file, line, name)
		else:
			if 'onlysample' in args.idformat:
				#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. May have to deal with multiple allele copies per sample with identical id names.
				for file in os.listdir(contigdir):
					if file.endswith('trimmed'):
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									print(line)
									linspl=line.split(' ')[0]
									linspl2=linspl.split('_')
									print(linspl2)
									name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
									print(name)
									replaceAll(file, line, name)
			else:
				sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")



# You can concatenate trimmed clustered-loci labelled deintereleaved_...._trimmed in the 'diploids' folder for maximum likelihood phylogeny (i.e. raxml, iqtree),
#resulting MSA's may have more than one sequence per sample and could be due to:
#retention of paralogous sequences, may want to set a higher -c2 identity, atleast 75% recommended.
#heterozygous variants were retained and not collapsed into consensus alleles, check clustering id for -c1, .99 recommended
#In our assessment multiple contigs is usually the result of two locus fragments that were not contiguous and unable to be joined in the consensus allele clustering step
#The purpose of this stage is to reduce multi-copy sequences due to paralogy by identity clustering among samples for the same reference locus.
#However, some datasets based on a variety of rich to poor DNA inputs, may have significantly higher multi-copy sequences at a given locus
# due to unjoined fragmented amplicons of the same locus which is common with poor sample DNA quality.
#we leave the decision up to the user as to how to manage extra contig/scaffold sequences present in locus clusters.(e.g. choose the longest sequence)
