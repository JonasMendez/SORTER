# SORTED
Sorter of Orthologous Regions for Target Enrichment Data

Jonas Mendez-Reneau1, Erin Sigel1
University of Louisiana Lafayette

Contact:jonasmrgrad@gmail.com

#############
Dependencies:
#############
Biopython (Python 2.7.x)
Perl
Trim Galore v 0.6.4 (with cutadapt v1.1B)
bwa v 0.7.17-r1188
SAMTools v 1.9 (using htslib 1.9)
BCFTools v1.9 (using htslib 1.9)
MAFFT v 7.450
SPAdes v 3.11.1
Seqtk v 1.3-r107-dirty
Vcfutils.pl (from bcftools v 1.9)
Trimal v 1.2rev59
Usearch v 11.0.667_i86linux32

########
Scripts:
########
Phase1.py
Phase2.py
Phase3.py
annotatedupes
deinterleave.py
getlongestcontig.py
keeplongest.py
seqclean.py

########################################
General Definition of Phases and Output:
########################################

Phase1.py
This script will generate clustered locus copies at a user set threshold (.60-.65 recommended) for all of your samples and output them as alignments (per locus-cluster).
The clustering is done to separate out potential paralogs because Target-Enrichment data will often amplify them in non-model systems.
If you are interested in using this pipeline to look at possible parentage in known allopolyploids/hybrids, or identify cryptic hybrids/allopolyploids
in your current dataset, use scripts Phase2.py and Phase3.py

Phase2.py
This script will take the locus-cluster dataset from Phase1.py and phase every single sample to output a phased alignment.
This phased dataset can be used to identify potential hybrid/allopolyploid species. Phased sequences for each sample should 
be sister to each other with relatively strong support in a phylogeny, if the samples are behaving as diploids. 
If you observe phased samples that don't result as sister, and are associated with poor node support values they likely indicate
allopolyploid/hybrid species. These samples should be removed from the dataset so that you can rerun Phase1.py without them and then
phase the hybrid/allopolyploid sample(s) using Phase3.py

Phase3.py
This script will phase and assign orthology to samples suspected to be allopolyploid or hybrid in origin. It uses the locus-cluster dataset
for diploids from Phase1.py as the reference to determine which locus-cluster the phased sequences belong to, and also to assign preliminary
species parentage to each phased sequence based on sequence similarity. The script outputs two main alignment files with extensions _phase.fasta and _diploidhit.fasta
the latter annotates allopolyploids/hybrids with the top diploid blast hit and the former randomly assigns _ph0 or _ph1 annotation to phased sequences.
The diploid hit annotation can be run as is to see where separate concatenated groupings result in a phylogeny. Closely related groups can then
be collapsed into one haplotype if they result as all being associated with the same diploid species, likewise other groupings associated with
other diploid species should be collapsed. You can then run another phylogenetic analyses with the final collapsed haplotypes.
The current version of this script will output an entire set of locus-cluster alignments + the phased/annotated sequences of ONE allopolyploid/hybrid,
future version will incorporate adding sequences for multiple allopolyploids into the same alignment.
The script will accomodate multiple samples at once, but because we are re-aligning the entire Phase1 locus-cluster dataset for EVERY allopolyploid/hybrid sample,
you may consider only doing 1-4 samples at a time if you have run-time constraints <12 hours.

#####################
Running the pipeline:
#####################
Phase1.py
#########

BEFORE RUNNING MAKE SURE YOU HAVE DONE THE FOLLOWING:
1. MAKE A FOLDER TO SERVE AS THE WORKING DIRECTORY (this will be the -wd FOR ALL PIPELINE SCRIPTS)
2. ADD RAW PAIR-END FASTQ READS FOR EACH SAMPLE IN WORKING DIRECTORY, FOLLOWING THIS NAMING SCHEME FOR PAIRED END READS, RESPECTIVELY:
	@@##_species1id_R1.fastq and @@##_species1id_R2.fastq
	Make sure paired end reads are signified with R1/R2 as shown above
	@@ values can be any set of two alphabetic letters and ## can be any unique set of two intergers per sample
	These Unique Identifiers are used to differentiate multiple samples of the same species (i.e. if multiple samples have the same 'species1id')
3. WITHIN WORKING DIRECTORY MAKE THE FOLLOWING FOLDERS:
	workingdirectory/diploids (this folder outputs Phase1 alignments)
	workingdirectory/diploidclusters (This folder serves as a reference for allopolyploid/hybrid BLAST steps and will store locus-clusters)
	workingdirectory/diploids_phased (if phasing with Phase2.py, this folder outputs Phase2 Phased alignments for diploids. 
		***YOU WILL STILL HAVE TO RUN Phase2.phy with -cdbonly T if you wish to use Phase3.py to phase allopolyploids/hybrids without phasing in Phase2.phy)
	workingdirectory/phaseset (if using Phase3.py, this folder will contain hybrid/allopolyploid sample files and alignments)
4. PLACE ALL EXTERNAL PYTHON SCRIPTS IN WORKING DIRECTORY:
	workingdirectory/getlongestcontig.py
	workingdirectory/deinterleave.py
	workingdirectory/seqclean.py
	workingdirectory/annotatedupes
	workingdirectory/keeplongest.py
5. Know the working directory of your probe references for command line input, may be placed in /workingdirectory/
		REFERENCE FILES MAY HAVE ALTERNATIVE FORMATING, BUT TARGET LOCI MUST BE THE BEGINNING OF THE STRING (after '>') FORMATTED AS L#_ WITH L1_ - LN_ (N=number of target reference loci) 
		AS FOLLOWS:
		>L1_*
		>L2_*
		>L3_*
		...
		>L20_*
		>L21_*
		...
		>L200_*
		etc... (*'s representing any other ID format present in your reference; The pipeline ignores this information and is just interested in identifying genes from references)
6. Before running any script make sure to load all required software and python version


FLAGS:

-wd WORKING DIRECTORY 
-ref PROBE REFERENCE DIRECTORY (FASTA FILE)
-c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .85 - .97 Recommended) 
-c2 2nd CLUSTER ID (AMONG SAMPLES FOR LOCUS-CLUSTERS, .55 - .65 Recommended, depends on taxonomic breadth of ingroup/outgroups; general guideline would be to set a higher id for 1-4 closely related outgroups, or lower for 4+ outgroups of varying genetic distance.) 
-trim RUN TRIMGALORE TO TRIM RAW READS? (T/F)***
-spades RUN SPADES ASSEMBLY? (T/F)***
-onlyprocess (T/F) ONLY RUN TRIMGALORE AND SPADES FOR CONTIG PROCESSING; RUN AGAIN WITH -trim and -spades as F FOR PIPELINE (set as F if running processing + pipeline in one run)
-reclust (T/F) IF T ONLY RERUNs LOCUS CLUSTERING + ALIGNMENT. CLEARS PREVIOUS LOCUS-CLUSTER OUTPUT, BACK UP IF NEEDED TO COMPARE CLUSTERING THRESHOLDS.
	RUNNING -reclust T IGNORES -c1 -cs -csl -csn VALUES, USING VALUES FROM PREVIOUS RUN; TO CHANGE THESE FLAGS YOU MUST RERUN Phase1.py WITH -reclust F  spades F -trimgalore F -op F
-cs TAKE SPADES CONTIGS or SCAFFOLDS? (input as: scaffold or contig)
-csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES
-csl ONLY TAKE CONTIGS LARGER THAN N LENGTH
-al number of iterations for MAFFT alignments (1000 recommended)
-indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)
-idformat (full/copies/onlysample/*) OUTPUTS FINAL ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:
	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. If last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes, discontinouos haplotype fragments or unclustered paralogs, may consider higher-c2 (among samples) or -c1 (within samples) clustering values
	-idformat copies = >@@##_sampleid_0 ; Keeps sample id and consensus allele copy count per cluster. i.e. if last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; see above
	-idformat onlysample = >@@##_sampleid ; Keeps only the sample id across locus-cluster alignments, easiest for concatenation across all locus-cluster alignments
	-format * if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp


***Make 'F' if you have already RUN TRIMGALORE or SPADES for your paired end reads
***If you have already run trimgalore, and want this script to run SPADES assembly, you must organize your working directory
so that each set of trimmed reads has the naming scheme @@##_species1_R1_val_1.fq and @@##_species1_R2_val_2.fq
Each set of trimmed .fq files is in their own folder corresponding to each sample, with each sample folder named as follows:
	/workingdirectory/WA01_species1_R1.fastq/
	/workingdirectory/WA02_species2_R1.fastq/
	/workingdirectory/WB03_species3_R1.fastq/
	etc...
folders have ..._R1.fastq extensions as unique identifiers for the script.

COMMAND LINE EXAMPLE:
python phase1.py -wd /workingdirectory/ -ref /workingdirectory/references.fasta -op F -spades T -trimgalore T -op F -reclust F -c1 .85 -c2 .60 -spades T -trimgalore T -cs contig -csn 6 -csl 350 -al 1000 -indel 0.25 -idformat onlysample

#########
Phase2.py
#########

PHASE 2 OF THE SORTED PIPELINE TAKES THE LOCUS-CLUSTER DATASET GENERATED IN PHASE 1 AND PHASES CONSENSUS ALLELES BY REMAPPING TRIMMED READS TO SAMPLE-SPECIFIC LOCUS-CLUSTERS.
MAKE SURE YOU HAVE MADE THE /workingdirectory/diploids_phased/ DIRECTORY BEFORE RUNNING.

FLAGS:

-wd WORKING DIRECTORY (DIRECTORY STRING STARTING AND ENDING IN '/'. e.g: "-wd /workingdirectory/" should be the same directory used in Phase 1) 
-pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)
-psl MINIMUM LENGTH OF PHASED SEQUENCES (recommended 350)
-n PROPORTION OF MISSING DATA (N) ALLOWED IN PHASED SEQUENCES? (atleast 50% bp representation recommended, input as -n 50 , NOT AS DECIMAL)
-al NUMBER OF ITERATIONS FOR MAFFT ALIGNMENTS (1000 recommended)
-indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)
-idformat (full/copies/onlysample/*) OUTPUTS FINAL FASTA ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:
	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. (Locus_ClusterID_@@##_sampleid_phase)
	-idformat phase = >@@##_sampleid_0 ; RECOMMENDED FOR PHASING. Keeps sample id and phase. i.e. if last annotation >1 signifies samples with multiple consensus alleles per locus-cluster that were phased.
	-idformat onlysample = >@@##_sampleid ; Keeps only the sample for sequence headers. Will have to decide how to manage phase or other sequence copies with identical id names.
	-idformat * :if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_@@##_sampleid_0 1230 bp

-cdbonly (T/F) set to 'T' to only make blast database for phasing and blasting allopolyploid/hyrbid samples in Phase3.py. 
	Use this setting if you wish to skip phasing your diploid locus-cluster samples but still want to process polyploids/hyrbids with Phase3.py

COMMAND LINE EXAMPLE
python Phase2.py -wd /working/directory/ -pq 20 -psl 350 -n 50 -al 1000 -indel .25 -idformat phase -cdbonly F

#########
Phase3.py
#########

BEFORE RUNNING MAKE SURE YOU HAVE DONE THE FOLLOWING:
1. MAKE A FOLDER TO SERVE AS THE WORKING DIRECTORY FOR THE SET OF ALLOPOLYPLOIDS OR HYBRIDS TO BE PHASED AS /workingdirectory/phaseset/
2. ADD RAW PAIR-END FASTQ READS FOR EACH SAMPLE IN PHASESET WORKING DIRECTORY, FOLLOWING THIS NAMING SCHEME FOR PAIRED END READS, RESPECTIVELY(SAME AS Phase1.py):
	@@##_species1id_R1.fastq and @@##_species1id_R2.fastq
	Make sure paired end reads are signified with R1/R2 as shown above
	@@ values can be any set of two alphabetic letters and ## can be any unique set of two intergers per sample, MUST BE IDENTICAL FOR THE SAME SAMPLE
	These Unique Identifiers are used to differentiate multiple samples of with the same species id (i.e. if multiple samples have the same 'species1id')
4. PLACE ALL EXTERNAL PYTHON SCRIPTS IN WORKING DIRECTORY (Same as Phase1.py):
	workingdirectory/getlongestcontig.py
	workingdirectory/deinterleave.py
	workingdirectory/seqclean.py
	workingdirectory/annotatedupes
	workingdirectory/keeplongest.py
5. Know the working directory of your probe references for command line input, may be placed in /workingdirectory/
6. Before running script make sure to load all required software and python version 

FLAGS

-wd WORKING DIRECTORY 
-ref PROBE REFERENCE DIRECTORY (FASTA FILE) SAME AS Phase1.py
-c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .85 - .90 Recommended) 
-trimgalore RUN TRIMGALORE TO TRIM RAW READS? (T/F)***
-spades RUN SPADES ASSEMBLY? (T/F)***
-onlyprocess (T/F) ONLY RUN TRIMGALORE AND SPADES FOR CONTIG PROCESSING; RUN AGAIN WITH -trim and -spades as F FOR PIPELINE (set as F if running processing + pipeline in one run)
-cs TAKE SPADES CONTIGS or SCAFFOLDS? (MUST INPUT AS: scaffold or contig)
-csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES. DEPENDING ON SUSPECTED PLOIDY, MULTIPLY BY 2 TO ACCOUNT FOR HETEROZYGOUS VARIANTS (i.e. tetraploid 4*2=8, triploid 3*2=6. 8-10 for unknown samples is recommended)
-csl ONLY TAKE CONTIGS/SCAFFOLDS LARGER THAN N LENGTH
-pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)
-n PERCENT OF MISSING DATA (N) ALLOWED IN PHASED SEQUENCES? (atleast 50% bp representation recommended, input as -n 50 , NOT AS DECIMAL)
-al number of iterations for MAFFT alignments (1000 recommended)
-indel indels have to be present in atleast XX% of sequences to be kept (0.25 recommended for ~50 samples, be aware of the number of samples you are processing)


***Make 'F' if you have already RUN TRIMGALORE or SPADES for your paired end reads
***If you have already run trimgalore, and want this script to run SPADES assembly, you must organize your working directory
so that each set of trimmed reads has the naming scheme W@@##_species1_R1_val_1.fq and W@@##_species1_R2_val_2.fq
Each set of trimmed .fq files is in their own folder corresponding to each sample, with each sample folder named as follows:
	/workingdirectory/WA01_species1_R1.fastq
	/workingdirectory/WA02_species2_R1.fastq
	/workingdirectory/WB03_species3_R1.fastq
	etc...
folders have ..._R1.fastq extensions as unique identifiers for the script.

Output Files:
/workingdir/phaseset/sample/diploidclusters/*_phase.fasta extension making polyploid/hybrid sample ids in the following format : >@@##_sampleid_phase(0/1)
/workingdir/phaseset/sample/diploidclusters/*_diploidhit.fasta extension making polyploid/hybrid sample ids in the following format : >@@##_sampleid_diploidhit
the *_trimmed.fasta files are annotated as @@##_sampleid_diploidhit_phase(0/1)

COMMAND LINE EXAMPLE:
python phase1.py -wd /workingdirectory/ -ref /workingdirectory/references.fasta  -spades T -trimgalore T -op F -c1 .85 -cs contig -csn 6 -csl 350 -pq 20 -n 50 -al 1000 -indel 0.25
