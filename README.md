# MAGspinner
A set of tools for decontaminating and reassembling Metagenome assembled genomes (MAGs). It is an iterative method that takes a MAG, maps reads to it, reassembles it and throws out scaffolds based on whether they differ significantly in tetranucleotide and abundance from the core MAG. It is called "MAGspinner" to illustrate that the MAG is conceptually spun in multidimensional space and loosely associated scaffolds will be ejected from the core MAG. It assumes that a MAG is relatively low in contamination to begin with, and therefore an average tetranucleotide signature can be recovered that is similar to the true signature of the MAG of interest. 

version: 0.1.0 (January 7, 2020) uploaded scripts that can be run in a Linux High performance computing environment using SLURM and with necessary software available as modules. R should be preloaded with necessary libraries installed.

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
merged read files = ../intermediateResults/BBmerged*${SAMPLE}*\
unpaired reads = ../intermediateResults/$SAMPLE\*unpaired\*


Craig Herbold, University of Vienna


This software is provided as-is with improvements added as I have time. The procedure was developed conceptually by Petra Pjevac and Dimitri Meier in the following manuscripts:\

Dyksma, S., Bischof, K., Fuchs, B.M., Hoffmann, K., Meier, D., Meyerdierks, A., Pjevac, P., Probandt, D., Richter, M., Stepanauskas, R. and Mußmann, M., 2016. Ubiquitous Gammaproteobacteria dominate dark carbon fixation in coastal sediments. The ISME journal, 10(8), pp.1939-1953.\
Mußmann, M., Pjevac, P., Krüger, K. and Dyksma, S., 2017. Genomic repertoire of the Woeseiaceae/JTB255, cosmopolitan and abundant core members of microbial communities in marine sediments. The ISME journal, 11(5), pp.1276-1281.\

Automation was developed into MAGspinner and described in the following manuscripts:\
<EMBARGO to be lifted Friday Feb. 7>\

Please consider citing these manuscripts if you use MAGspinner before it is developed into a stand-alone manuscript.\
