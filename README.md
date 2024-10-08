# SORTER
### Sorter of Orthologous Regions for Target Enrichment Reads

Jonas Mendez-Reneau1, Erin Sigel2, Gordon Burleigh3
1. University of Louisiana Lafayette
2. University of New Hampshire, Durham
3. University of Florida, Gainesville

If you use SORTER for your study, please cite the following publication:

Jonas Mendez-Reneau, J. Gordon Burleigh, and Erin M. Sigel "Target Capture Methods Offer Insight into the Evolution of Rapidly Diverged Taxa and Resolve Allopolyploid Homeologs in the Fern Genus Polypodium s.s.," Systematic Botany 48(1), 96-109, (28 March 2023). https://doi.org/10.1600/036364423X16758873924135

As of August 2024, we have released a newer version of SORTER; The SORTER2-Toolkit has been optimized for user set-up, improved compatibility, additional summary statistics, and new tools for identifying hybrids and progenitors. We encourage users to use SORTER2 going forward, visit: https://github.com/JonasMendez/SORTER2

SORTER is a flexible user-customizable set of python scripts for building a variety of locus-alignment matrices from target-enrichment datasets. The three overarching goals of the pipeline are to:

1. Builds multiple sequence alignments of orthologous sequences for loci generated from paired-end reads in target-enrichment datasets. The initial purpose of this pipeline was to reduce multi-copy sequences due to heterozygosity or paralogy by identity clustering within and among samples for the same reference locus. The pipeline generates consensus alleles from heterozygous variants by clustering contigs associated with the same reference locus and sample to generate consensus sequences presumably representing allelic variation at IUPAC ambiguity sites. The pipeline then maps consensus alleles for all samples to the same reference and does a second round of clustering among samples to separate potential paralogs. Given appropriate clustering identity, this can filter and separate paralogs derived from the same locus into separate orthologous sets across all samples, effectively generating additional loci for analysis when paralogs are present. 

2. Phase bi-allelic variation from previously determined orthologs to build alignments from putative allelic haplotypes where each sample is represented by two sequences representing heterozygous or homozygous alleles. i.e. consensus allele sequences are phased into respective bi-allelic haplotypes for each presumably diploid sample.

3. Infer hybrid haplotypes based on similarity to potential progenitor samples. This generates a multiple sequence alignment where hybrid samples can have two or more tips, depending on the number of hybrid haplotypes present in the sample. This was originally designed and probably works best for allopolyploid hybrids where we expect lower levels of inter-homeologous recombination, but we have had success inusing the pipeline for detecting hybrid haplotypes from putative homoploid hybrids.

Contact:jonasmrgrad@gmail.com



### Dependencies:

Biopython (Python 2.7.x)

Perl

Trim Galore v 0.6.4 (with cutadapt v1.1B and FastQC v0.12.1)

bwa v 0.7.17-r1188

SAMTools v 1.9 (using htslib 1.9)

BCFTools v1.9 (using htslib 1.9)

MAFFT v 7.450

SPAdes v 3.11.1

Seqtk v 1.3-r107-dirty

Vcfutils.pl (from bcftools v 1.9)

Trimal v 1.2rev59

Usearch v 11.0.667_i86linux32


### Scripts:

Stage1A.py

Stage1B.py

Stage2.py

Stage3.py

annotatedupes

deinterleave.py

getlongestcontig.py

keeplongest.py

seqclean.py


## General Definition of Stages and Output:

### Stage1A.py

The first script trims raw reads (trimgalore) and builds contigs (SPADES) for each sample, making folders for each sample to be used downstream.

### Stage1B.py

Stage1B.py takes contigs generated by spades, first generating "consensus alleles" by collapsing highly similar contigs mapping to the same reference within a sample into consensus sequences. An identity threshold of 98-99% identity (-c1 flag) is recommended for generating consensus alleles. Consensus alleles among all samples are compiled according to reference locus for a second round of clustering (-c2 flag) to filter potentially paralogous locus sets. Sensitivity analyses of a range of identities [i.e. 70-90%] is recommended to assess the relative effect of consensus allele copy number per locus-cluster, but Edgar (2010) recommends atleast 75% identity thresholds for retrieving orthologs.

Stage1B.py outputs aligned fasta files of consensus alleles for clustered loci in the "diploids" folder with files ending in "trimmed", we leave multi-copy sequences in the alignment, leaving the choice up to the user as to how to filter samples with multiple sequences per locus (e.g. choose the longest sequence). A csv summary table for all retrieved consensus-alleles per sample per locus prior to among sample clustering can be found in the "diploids" folder, along with consensus allele files corresponding to references and samples. A secondary summary csv file found in the "diploidclusters" shows consensus-allele distributions across samples and loci after clustering. 

### Stage2.py
This script will take the putatively orthologous clustered loci for samples generated from Stage1B.py and maps sample reads to them in order to phase (samtools phase) putative bi-allelic haplotypes, assuming all samples are diploid. Output in "diploids_phased" folder results in phased alignments for each locus, where each sample has two allelic haplotypes associated with it. 

This phased dataset can potentially be used to identify hybrid/allopolyploid samples. Phased sequences for each sample should be sister to each other with relatively strong support in a phylogeny, or multiple samples of the same species should at the very least be monophyletic if samples are behaving as typical non-hybrid diploids. 
If you observe phased samples that don't result as sister to themselves or to other samples of its expected taxon and are associated with poor node support values they may indicate a hybrid individual. You need to filter potential hybrids if you want to use Stage3.py to infer hybrid haplotypes, as the inference is highly dependent on correct identification and labeling of potential progenitors. (i.e. accidentally including a hybrid as a progenitor could reduce the ability for the pipeline to accurately identify differing hybrid haplotypes)

### Stage3.py
This script will phase and assign orthology to sequences in samples suspected to be allopolyploid or hybrid in origin, outputing multiple-sequence alignments with hybrids represented by two or more pairs of phased alleles, depending on how many progenitors species contributed to the hybrid. It uses the putatively orthologous and phased datasets generated from Stage1B.py and Stage2.py as progenitor references to assign identities to different genomes or haplotypes found in a hybrid. Thus, this script assumes that the samples included through Stage1B.py and Stage2.py will serve as possible progenitor contributors to each hybrid sample. Due to this inference can be highly impacted by having miss-identified samples or cryptid hybrids included in the analysis. The script outputs alignment files, where hybrid sequences are annotated with the name of putative progenitors, allowing for analyses linking hybrid haplotypes across loci without a priori knowlegde of linkage or genomes. Disparate sequence sets may represent the same haplotype that was inferred to another species due to ILS or gene level heterogeneity, so we recommend preliminary phylogenetic analyses to assess if disparate haplotypes should be further combined (i.e. the two sets result in a shared and well supported clade)

Stage3_blastoconsensus.py uses the consensus alleles generated in Stage1.py, rather than phased sequences from Stage2.py, as the reference sequences to build a UBLAST database for inferring hybrid haplotypes and does not annotate hybrid phase pairs (i.e. haplotype pairs associated with the same homeolog). Additionally, this script generates separate datasets for each hybrid sample included, meaning each hybrid is aligned to its own the progenitor dataset (Be careful of computational wall-time if using this script to process several hybrids, as you are re-alligning the entire progenitor dataset for each hybrid sample sequentially). This approach is amenable to analyses where you are interested in only retaining the longest or most infromative phased allele for each homeolog, and only including one hybrid sample per alignment dataset. This was the method used in Mendez-Reneau et al. (2023) to faciliatate comparison of datasets with different allopolyploid species and filtering of loci with only one homeolog pair represented.


## Running the pipeline:

### Stage1A.py
1. Place all paired raw reads for all samples in the same folder.
2. Paired reads for each sample need to follow this naming scheme for the pipeline to work:

	@@##_taxon_R1.fastq and @@##_taxon_R2.fastq
	
	@@## can be any set of numbers and or letters used as unique identifiers for each sample, do not use underscores or other symbols other than numbers or letters
3. Trimmed reads and contigs for each sample will be placed in their own folder within the "cleanfastq" folder that is made. This will serve as the working directory for the rest of the scripts.

Flags:

-wd : This will be the path to the folder where you have stored the labeled paired-end reads. make sure you include '/' at the end of the path.

command line example:

python Stage1A.py -wd /path/to/paired/reads/

### Stage1B.py
BEFORE RUNNING MAKE SURE YOU HAVE DONE THE FOLLOWING:

Before running any more scripts, prepare the working directory (i.e. "cleanfastq" folder generated in Stage1A.py) with required files and folders.

1. CREATE THE FOLLOWING FOLDERS:

	workingdirectory/cleanfastq/diploids (this folder outputs Stage1B.py alignments)
	
	workingdirectory/cleanfastq/diploidclusters (un-aligned clustered loci are stored here)
	
	workingdirectory/cleanfastq/diploids_phased (if phasing with Stage2.py, this folder outputs phased alignments.)
	
	workingdirectory/cleanfastq/phaseset (if using Stage3.py, this folder will contain hybrid/allopolyploid sample files and alignments)
	

2. PLACE ALL EXTERNAL SCRIPTS IN WORKING DIRECTORY:

	workingdirectory/cleanfastq/getlongestcontig.py
	
	workingdirectory/cleanfastq/deinterleave.py
	
	workingdirectory/cleanfastq/seqclean.py
	
	workingdirectory/cleanfastq/annotatedupes
	
	workingdirectory/cleanfastq/keeplongest.py
	

3. Know the working directory of your probe references for command line input, may be placed in /workingdirectory/cleanfastq/

4. Locus references should be in FASTA format. Reference IDs may have any labels, but each reference locus must be annotated with "L#_..." at the begnning of the sequence ID (i.e. directly after the ">") to generate an index for all reference loci to be used by the pipeline. Multiple references may be used for the same targeted locus, but they must be annotated with the same locus index number.

AS FOLLOWS:
		
	>L1_*
	>L2_*
	>L3_*
	...
	>L20_*
	>L21_*
	...
	>L200_*
	etc... 
		
(*'s representing any other ID format present in your reference; The pipeline ignores this information and is just interested in identifying the index number from references)

5. Before running any script make sure to load all required software and python version


Flags:

-wd WORKING DIRECTORY : -wd /workingdirectory/cleanfastq/ (DIRECTORY PATH STRING MUST START AND END IN '/')

-ref  PATH TO REFERENCE (FASTA FILE): -ref /path/to/ref.fasta


-loci NUMBER OF REFERENCE LOCI (e.g. -loci 550)

-c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES, .97-.99 Recommended) 

-c2 2nd CLUSTER ID (AMONG SAMPLES FOR LOCUS-CLUSTERS, .70 - .90 recommended, may depend on taxonomic breadth of ingroup/outgroups.


-reclust (T/F) IF T, ONLY RERUNs LOCUS CLUSTERING + ALIGNMENT, WILL NOT REBUILD CONSENSUS ALLELES (i.e. keeps -cs, -csn, and -csl output from previous run, rerun with -reclust F if you wish to change these parameters). CLEARS PREVIOUS CLUSTERED LOCI, BACK UP OUTPUT IF NEEDED TO COMPARE CLUSTERING THRESHOLDS.


-cs TAKE SPADES CONTIGS or SCAFFOLDS? (input as: scaffold or contig)

-csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES/ORTHOLOG CLUSTERING

-csl ONLY TAKE CONTIGS/SCAFFOLDS LARGER THAN N LENGTH

-al number of iterations for MAFFT alignments, --maxiter option (e.g. -al 1000)

-indel indels have to be present in atleast XX% of sequences to be kept; trimal -gt option

-idformat (full/copies/onlysample/*) OUTPUTS FINAL ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:

	-idformat full = >L100_cl0_@@##_sampleid_0 ; Keeps full annotation. If last annotation > 0, it signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes, discontinouos haplotype fragments or unclustered paralogs, may consider higher-c2 (among samples) or -c1 (within samples) clustering values.
	
	-idformat copies = >@@##_sampleid_0 ; Keeps sample id and consensus allele copy count per cluster. i.e. if last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; see above
	
	-idformat onlysample = >@@##_sampleid ; Keeps only the sample id across locus-cluster alignments, easiest for concatenation across all locus-cluster alignments
	
	-idformat * = if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp


Command line example:

python Stage1B.py -wd /workingdirectory/cleanfastq/ -ref /workingdirectory/cleanfastq/references.fasta -loci 450 -reclust F -c1 .99 -c2 .80 -cs contig -csn 10 -csl 300 -al 1000 -indel 0.25 -idformat onlysample


### Stage2.py

Stage2.py takes the dataset generated in Stage1B.py, maps reads to consensus alleles in order to phase haplotypes using samtools phase, resulting in two haplotype bi-allelic sequences per sample, assuming all samples are diploid. The phased clustered loci are then re-aligned for analysis and output in the 'diploids_phased' folder.

*** MAKE SURE YOU HAVE MADE THE /workingdirectory/cleanfastq/diploids_phased/ DIRECTORY BEFORE RUNNING.

Flags:

-wd WORKING DIRECTORY : -wd /workingdirectory/cleanfastq/ (DIRECTORY PATH STRING MUST START AND END IN '/') should be the same directory used in Stage1B.py

-pq PHASE QUALITY; samtools phase -q flag (atleast 20 recommended)

-al number of iterations for MAFFT alignments, --maxiter option (e.g. -al 1000)

-indel indels have to be present in atleast XX% of sequences to be kept

-idformat (full/phase/onlysample/*) OUTPUTS FINAL FASTA ALIGNMENT SEQUENCE IDS IN FOLLOWING FORMATS:

	-idformat full = >L100_cl0_@@##_sampleid_ph0/ph1_0 ; Keeps full anottation.  If last annotation > 1, it signifies samples with multiple consensus alleles per locus-cluster that were phased.
	
	-idformat phase = >@@##_sampleid_ph0/ph1 ; RECOMMENDED FOR PHASING. Keeps sample id and phase annotations.
	
	-idformat onlysample = >@@##_sampleid ; Keeps only the sample for sequence headers, removing phase annotation. You will have to decide how to manage phased or other sequence copies with identical id names.
	
	-idformat * :if you mispell the above arguments or leave -id format blank, it will keep the default trimal headers; e.g. >L100_cl0_@@##_sampleid_0 1230 bp

Command line example:

python Stage2.py -wd /working/directory/ -pq 20 -al 1000 -indel .25 -idformat phase


### Stage3.py

This stage assumes samples processed through stages 1B - 2 are possible progenitor taxa contributing to putative hybrid samples. Specifically, phased allopolyploid (or unknown hybrid) sequences are UBLASTED (usearch) to phased progenitor sequences from Stage2.py in order to assess potential hybrid parentage. For best performance, a broad or informed sampling of possible progenitors as well as verified identifications of taxa and hybrids are critical. Mis-identified taxa or cryptic hybrids can confound results by giving the wrong signal associated with an erroneous taxon label, thus hybrid haplotypes might be associated with the wrong lineage and or erroneously over-represented in others causing mixed hybrid sequence sets. To improve inference, the pipeline automatically filters sequences with >50% missing data.

#### Preparing hybrid samples and running Stage3.py

If you have not already done so, run stages 1B - 2 on the progenitor samples you wish to use to compare hybrids to, excluding any other hybrid samples or potentially mis-identified samples (or you can fix the label for mis-identified samples by renaming their fastq files and sample folder name to the appropriate label, and re-run stages 1 - 2 so that the dataset has the correct taxon labels). If you have not already done so, process the hybrid samples you want to test using Stage1A.py first to trim reads and build contigs in a separate folder. You can then move the output hybrid sample folders into the "phaseset" folder associated with the "cleanfastq" directory containing progenitor samples you previously clustered and phased with stages 1B and 2.

Flags:

-wd WORKING DIRECTORY 

-ref PROBE REFERENCE DIRECTORY (FASTA FILE) SAME AS Stage1B.py

-loci NUMBER OF REFERENCE LOCI SAME AS Stage1B.py

-c1 1st CLUSTER ID (WITHIN SAMPLE FOR CONSENSUS ALLELES,  0.98 - .99 Recommended) 

-cs TAKE SPADES CONTIGS or SCAFFOLDS? (MUST INPUT AS: scaffold or contig)

-csn TAKE N CONTIGS/SCAFFOLDS PER LOCUS PER SAMPLE FOR CONSENSUS ALLELES. DEPENDING ON SUSPECTED PLOIDY, MULTIPLY BY 2 TO ACCOUNT FOR HETEROZYGOUS VARIANTS (i.e. tetraploid 4*2=8, triploid 3*2=6. A high number like 20-30 may help for allopolyploids if paralogs are prevalent in dataset)

-csl ONLY TAKE CONTIGS/SCAFFOLDS LARGER THAN N LENGTH

-pq PHASE QUALITY; SAMTOOLS -Q FLAG; MINIMUM READS TO CALL A PHASE (atleast 20 recommended)

-al number of iterations for MAFFT alignments (1000 recommended)

-indel indels have to be present in atleast XX% of sequences to be kept (INPUT AS DECIMAL)


Output alignment files stored in:

/workingdir/phaseset/diploidclusters/


Hybrids are labeled as: @@##_sampleid_diploidhit_phase(0/1)

Where diploid hit is the name of the most similar progenitor and phase matches the bi-allelic haplotype phases for each locus.

Command line example:

python Stage3.py -wd /workingdirectory/cleanfastq/ -ref /workingdirectory/cleanfastq/references.fasta -loci 450 -c1 .99 -cs contig -csn 20 -csl 300 -pq 20 -al 1000 -indel 0.20



