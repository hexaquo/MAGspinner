#!/bin/bash
#SBATCH --job-name binReassembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --output=binReassembly-%j.out
#SBATCH --error=binReassembly-%j.err

#sbatch ../Scripts/binReassembly.sh <BIN> <BASE>

#module load spades
#module load bbmap
#module load checkm

echo "running job on ${SLURM_JOB_NODELIST}"

# inputs are the Base of the sample read set  
BIN=$1
SAMPLE=$2

#define directory where the script is called                                                                                                                           
HOMEFOLDER=`pwd`

#define and construct workspace on /tmp/$USER                                                                                                            
THISWORKFOLDER=$TMPDIR

#sync up data
# Ensure the job fails in case the lock could not be acquired
set -e
    (
	# Wait for the lockfile for max. 1 hour (3600 s), to not block the queue forever in case of dead lock files.
	flock -w 3600 200
	# Synchronize the relevant files to the local temporary directory
	zgrep -h "^" ../intermediateResults/BBunmerged*${SAMPLE}_1* > ${THISWORKFOLDER}/${SAMPLE}.R1.fastq
	
    ) 200>~/${SAMPLE}.R1.fastq.lock
    (
	# Wait for the lockfile for max. 1 hour (3600 s), to not block the queue forever in case of dead lock files.
	flock -w 3600 200
	# Synchronize the relevant files to the local temporary directory
	zgrep -h "^" ../intermediateResults/BBunmerged*${SAMPLE}_2* > ${THISWORKFOLDER}/${SAMPLE}.R2.fastq
    ) 200>~/${SAMPLE}.R2.fastq.lock
    (
	# Wait for the lockfile for max. 1 hour (3600 s), to not block the queue forever in case of dead lock files.
	flock -w 3600 200
	# Synchronize the relevant files to the local temporary directory
	zgrep -h "^" ../intermediateResults/BBmerged*$SAMPLE* >> ${THISWORKFOLDER}/${SAMPLE}.merged.fastq
    ) 200>~/${SAMPLE}.merged.fastq.lock
    (
	# Wait for the lockfile for max. 1 hour (3600 s), to not block the queue forever in case of dead lock files.
	flock -w 3600 200
	# Synchronize the relevant files to the local temporary directory
	zgrep -h "^" ../intermediateResults/$SAMPLE*unpaired* > ${THISWORKFOLDER}/${SAMPLE}.unpaired.fastq
    ) 200>~/${SAMPLE}.unpaired.fastq.lock
    (
	# Wait for the lockfile for max. 1 hour (3600 s), to not block the queue forever in case of dead lock files.
	flock -w 3600 200
	# Synchronize the relevant files to the local temporary directory
	rsync $BIN $THISWORKFOLDER
    ) 200>${BIN}.lock
set +e

BIN=${BIN##*/}

#move to work folder
cd ${THISWORKFOLDER}

gunzip *gz
mkdir ${BIN%.gz}.reassembledBins
mkdir ${BIN%.gz}.reassembledBins.temp
echo "Bin Id,Marker lineage,# genomes,# markers,# marker sets,Completeness,Contamination,Strain heterogeneity,Genome size (bp),# ambiguous bases,# scaffolds,# contigs,N50 (scaffolds),N50 (contigs),Mean scaffold length (bp),Mean contig length (bp),Longest scaffold (bp),Longest contig (bp),GC,GC std (scaffolds > 1kbp),Coding density,Translation table,# predicted genes,0,1,2,3,4,5+" >${BIN%.gz}.reassembledBins/${BIN%.gz}.checkM.stats.txt

THISROUND=1
STOPGO="GO"
while [ "${STOPGO}" == "GO" ]
do
    for (( round=$THISROUND; round<=$((THISROUND+9)); round++))
    do
	module unload spades
	module load bbmap
	bbmap.sh ref=${BIN%.gz} threads=$SLURM_CPUS_PER_TASK -Xmx40g
	bbmap.sh in1=${SAMPLE}.R1.fastq in2=${SAMPLE}.R2.fastq outm=${BIN%.gz}.PE.98.fastq minid=0.98 idfilter=0.98  threads=$SLURM_CPUS_PER_TASK -Xmx20g
	bbmap.sh in=${SAMPLE}.unpaired.fastq outm=${BIN%.gz}.unpaired.98.fastq minid=0.98 idfilter=0.98  threads=$SLURM_CPUS_PER_TASK -Xmx20g
	bbmap.sh in=${SAMPLE}.merged.fastq outm=${BIN%.gz}.merged.98.fastq minid=0.98 idfilter=0.98  threads=$SLURM_CPUS_PER_TASK -Xmx20g
	module unload bbmap
	module load spades
	spades.py --12 ${BIN%.gz}.PE.98.fastq -s ${BIN%.gz}.unpaired.98.fastq --merged ${BIN%.gz}.merged.98.fastq -t $SLURM_CPUS_PER_TASK -k 21,33,55,77,99,127 --careful --only-assemble -o $SAMPLE.assembly --trusted-contigs ${BIN%.gz}
	module unload spades
	cp ${BIN%.gz} ${BIN%.gz}.reassembledBins.temp/${BIN%.gz}.round.$((round-1)).fa
	mv ${BIN%.gz} ${BIN%.gz}.reassembledBins/${BIN%.gz}.round.$((round-1)).fa
	cp $HOMEFOLDER/../Scripts/binFilter.R ./
	cat $SAMPLE.assembly/scaffolds.fasta > ${BIN%.gz}
	# map to new bin for cleanup
	module load bbmap
	module load samtools
	bbmap.sh ref=${BIN%.gz} threads=$SLURM_CPUS_PER_TASK -Xmx40g
	bbmap.sh in=${BIN%.gz}.PE.98.fastq out=${BIN%.gz}.pe.98.sam mappedonly=t minid=0.98 idfilter=0.98 ambiguous=all bamscript=bs.sh threads=$SLURM_CPUS_PER_TASK -Xmx20g
	sh bs.sh
	bbmap.sh in=${BIN%.gz}.unpaired.98.fastq out=${BIN%.gz}.unpaired.98.sam mappedonly=t minid=0.98 idfilter=0.98 ambiguous=all bamscript=bs.sh threads=$SLURM_CPUS_PER_TASK -Xmx20g
	sh bs.sh
	bbmap.sh in=${BIN%.gz}.merged.98.fastq out=${BIN%.gz}.merged.98.sam mappedonly=t minid=0.98 idfilter=0.98 ambiguous=all bamscript=bs.sh threads=$SLURM_CPUS_PER_TASK -Xmx20g
	sh bs.sh

	#identify contigs that are connected
	if [ -f ${BIN%.gz}.scaffoldConnections.txt ]; then rm ${BIN%.gz}.scaffoldConnections.txt; fi
	samtools view ${BIN%.gz}.pe.98_sorted.bam | awk '{if(length(readArray[$1])==0){readArray[$1]=$5}else{readArray[$1]=readArray[$1]"\t"$5}}END{for(entry in readArray){print readArray[entry]}}' | awk '{for(i=1;i<NF;i++){for(j=i+1;j<=NF;j++){if($i<$j){thisLineArray[$i"\t"$j]=1};if($j<$i){thisLineArray[$j"\t"$i]=1}}}for(entry in thisLineArray){print entry};delete thisLineArray}' >> ${BIN%.gz}.scaffoldConnections.txt
	samtools view ${BIN%.gz}.merged.98_sorted.bam | awk '{if(length(readArray[$1])==0){readArray[$1]=$5}else{readArray[$1]=readArray[$1]"\t"$5}}END{for(entry in readArray){print readArray[entry]}}' | awk '{for(i=1;i<NF;i++){for(j=i+1;j<=NF;j++){if($i<$j){thisLineArray[$i"\t"$j]=1};if($j<$i){thisLineArray[$j"\t"$i]=1}}}for(entry in thisLineArray){print entry};delete thisLineArray}' >> ${BIN%.gz}.scaffoldConnections.txt
	cat ${BIN%.gz}.scaffoldConnections.txt | awk '{if($1!=$2){if($1<$2){connectionArray[$1"\t"$2]+=1}else{connectionArray[$2"\t"$1]+=1}}}END{for(entry in connectionArray){print entry"\t"connectionArray[entry]}}' > ${BIN%.gz}.scaffoldConnectionCounts.txt

	samtools depth ${BIN%.gz}.pe.98_sorted.bam > ${BIN%.gz}.pe.98_sorted.bam.depth.txt
	samtools depth ${BIN%.gz}.unpaired.98_sorted.bam > ${BIN%.gz}.unpaired.98_sorted.bam.depth.txt
	samtools depth ${BIN%.gz}.merged.98_sorted.bam > ${BIN%.gz}.merged.98_sorted.bam.depth.txt
	module unload samtools
	module unload bbmap

	#remove sequences that have outlier tetra or coverage
	R<binFilter.R --no-save ${BIN%.gz}
	mv ${BIN%.gz}.clean ${BIN%.gz}
	
	#save connected Contigs
	grep ">" ${BIN%.gz} | tr -d ">" >cleanNames
	cat ${BIN%.gz}.scaffoldConnectionCounts.txt | grep -F -f cleanNames | awk '$3>2{print $1"\n"$2}' | sort | uniq > connectedScaffolds.txt # use connections with a depth of at least 3

	#HARD rescue cutoff only keep contigs >500 nt despite connectivity
#	HRS=500
#	grep -F -f connectedScaffolds.txt removeContigs.txt | awk '{print $NF}' >rescuedContigs.txt
#	cat ${BIN%.gz}.removed | awk '{if(substr($1,1,1)==">"){printf("\n"$0"\t")}else{printf($0)}}END{print("")}' | grep -F -f rescuedContigs.txt | awk -v HRS=$HRS '{if(length($2)>HRS){print $1"\n"$2}}'>rescuedContigs.fa
#	grep ">" rescuedContigs.fa | tr -d '>' >rescuedContigs.txt
#	cat rescuedContigs.fa >>${BIN%.gz}

	mv ${BIN%.gz}.scaffoldConnectionCounts.txt ${BIN%.gz}.reassembledBins/${BIN%.gz}.scaffoldConnectionCounts.$((round)).txt
	mv ${BIN%.gz}.removed ${BIN%.gz}.reassembledBins/${BIN%.gz}.removed.$((round)).fa
	mv removeContigs.txt ${BIN%.gz}.reassembledBins/${BIN%.gz}.removeContigs.$((round)).txt
	mv rescuedContigs.txt ${BIN%.gz}.reassembledBins/${BIN%.gz}.rescuedContigs.$((round)).txt
	mv rescuedContigs.fa ${BIN%.gz}.reassembledBins/${BIN%.gz}.rescuedContigs.$((round)).fa
        mv ${BIN%.gz}.binFilter.pdf ${BIN%.gz}.reassembledBins/${BIN%.gz}.round.$((round)).fa.binFilter.pdf
	rm -r $SAMPLE.assembly
    done
    THISROUND=$round
    cp ${BIN%.gz} ${BIN%.gz}.reassembledBins.temp/${BIN%.gz}.round.$((round-1)).fa
    cp ${BIN%.gz} ${BIN%.gz}.reassembledBins/${BIN%.gz}.round.$((round-1)).fa

    module load checkm
    echo "Bin Id,Marker lineage,# genomes,# markers,# marker sets,Completeness,Contamination,Strain heterogeneity,Genome size (bp),# ambiguous bases,# scaffolds,# contigs,N50 (scaffolds),N50 (contigs),Mean scaffold length (bp),Mean contig length (bp),Longest scaffold (bp),Longest contig (bp),GC,GC std (scaffolds > 1kbp),Coding density,Translation table,# predicted genes,0,1,2,3,4,5+" >${BIN%.gz}.reassembledBins.temp/${BIN%.gz}.checkM.stats.txt
    checkm lineage_wf -x fa -t $SLURM_CPUS_PER_TASK --pplacer_threads $SLURM_CPUS_PER_TASK ${BIN%.gz}.reassembledBins.temp tempOut
    checkm qa --tab_table -o 2 tempOut/lineage.ms tempOut | grep -v "^\[" | awk 'BEGIN{FS="\t";OFS=","};{if(NR>1){$1=$1".fa";print $0}}' | tr '\t' ',' > checkM.temp
    rm -r tempOut
    cat checkM.temp >> ${BIN%.gz}.reassembledBins.temp/${BIN%.gz}.checkM.stats.txt
    cat checkM.temp >> ${BIN%.gz}.reassembledBins/${BIN%.gz}.checkM.stats.txt
    module unload checkm

    cd ${BIN%.gz}.reassembledBins.temp/
    cp $HOMEFOLDER/../Scripts/reassemblyEvaluation.R ./
    STOPGO=`R<reassemblyEvaluation.R --no-save ${BIN%.gz}.checkM.stats.txt | tail -n 2 | head -n 1`
    cd ../
    mv ${BIN%.gz}.reassembledBins.temp/reassemblyEvaluationStats.txt ${BIN%.gz}.reassembledBins/reassemblyEvaluationStats.$((THISROUND-11)).$((THISROUND-1)).txt
    rm ${BIN%.gz}.reassembledBins.temp/*
done
mv *bam ${BIN%.gz}.reassembledBins
mv *bam*txt ${BIN%.gz}.reassembledBins
mv *98.fastq ${BIN%.gz}.reassembledBins

gzip ${BIN%.gz}.reassembledBins/*
if [ ! -d ${HOMEFOLDER}/../intermediateResults ]
then
    mkdir ${HOMEFOLDER}/../intermediateResults
fi
if [ ! -d ${HOMEFOLDER}/../intermediateResults/${BIN%.gz}.reassembledBins ]
then
    mkdir ${HOMEFOLDER}/../intermediateResults/${BIN%.gz}.reassembledBins
fi
mv ${BIN%.gz}.reassembledBins/* ${HOMEFOLDER}/../intermediateResults/${BIN%.gz}.reassembledBins
