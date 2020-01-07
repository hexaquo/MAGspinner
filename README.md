# MAGspinner
A set of tools for decontaminating and reassembling Metagenome assembled genomes (MAGs). It is an iterative method that takes a MAG, maps reads to it, reassembles it and throws out scaffolds based on whether they differ significantly in tetranucleotide and abundance from the core MAG. It is called "MAGspinner" to illustrate that the MAG is conceptually spun in multidimensional space and loosely associated scaffolds will be ejected from the core MAG. It assumes that a MAG is relatively low in contamination to begin with, and therefore an average tetranucleotide signature can be recovered that is similar to the true signature of the MAG of interest. 

version: 0.1 (January 7, 2020) uploaded scripts that can be run in a Linux High performance computing environment using SLURM and with necessary software available as modules. R should be preloaded with necessary libraries installed.

Necessary software:\
spades\
bbmap\
samtools\
checkM

Necessary R libraries:\
Biostrings\
scales

Assumed data available for input is MAG, and metagenomic reads that have been merged and filtered through BBduk:\
BIN=Initial MAG\
SAMPLE=metagenomic sample identifier (e.g. SAMPLEB15)

Total metagenomic reads (paired, fastq format) are assumed to be in a folder with relative path ../intermediateResults/ from where Script is called and have following regular expression compatible names:\
unmerged read 1 files = ../intermediateResults/BBunmerged*${SAMPLE}_1*\
unmerged read 2 files = ../intermediateResults/BBunmerged*${SAMPLE}_2*\
merged read files = ../intermediateResults/BBmerged*$SAMPLE*\
unpaired reads = ../intermediateResults/$SAMPLE*unpaired*

