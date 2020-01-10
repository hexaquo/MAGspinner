library(Biostrings)
library("scales")

counts<-NULL
countGenomes<-NULL
countContigs<-NULL
countsForDistribution=NULL
alpha=0.05
contigEvalLength=2000

#Pass the fasta file as a command line argument - fasta file must be folded
bin=commandArgs()[3]

foldFile=gsub(".fa",".fold",bin)
#fold data
system(paste("fold -w 1000 ",bin," > ",foldFile, sep=""))
#---Read sequences into a DNAString variable---#
sequence = readDNAStringSet(foldFile, format="fasta")
originalSequences=sequence

unpairedMappingDepth=try(read.table(file=paste(bin,".unpaired.98_sorted.bam.depth.txt",sep="")))
peMappingDepth=try(read.table(file=paste(bin,".pe.98_sorted.bam.depth.txt",sep="")))
mergedMappingDepth=try(read.table(file=paste(bin,".merged.98_sorted.bam.depth.txt",sep="")))

#---identify short contigs ---
removeContigs_length=names(sequence)[which(width(sequence)<contigEvalLength)]
if(length(removeContigs_length)>0){
	unpairedMappingDepth=unpairedMappingDepth[which(is.na(match(unpairedMappingDepth[,1],removeContigs_length))),]
	peMappingDepth=peMappingDepth[which(is.na(match(peMappingDepth[,1],removeContigs_length))),]
	mergedMappingDepth=mergedMappingDepth[which(is.na(match(mergedMappingDepth[,1],removeContigs_length))),]
	sequence=sequence[-which(width(sequence)<contigEvalLength)]
}
removeContigs=removeContigs_length

#---Compute tetranucleotide frequencies---#
tetranucleotideCounts <- (oligonucleotideFrequency(sequence,width=4,with.labels=TRUE,as.prob=FALSE) + oligonucleotideFrequency(reverseComplement(sequence),width=4,with.labels=TRUE,as.prob=FALSE))
tetranucleotide<-(1+tetranucleotideCounts)/(dim(tetranucleotideCounts)[2]+rowSums(tetranucleotideCounts))
rownames(tetranucleotide)=names(sequence)
tetranucleotidePcoA=prcomp(log10(tetranucleotide),center=TRUE,scale=TRUE)

# now create second tetra for modeled distribution
tetraDistributionSet=NULL
for(i in c(1:length(sequence))){
	for(j in sample((length(sequence[[i]])-contigEvalLength),(length(sequence[[i]])-contigEvalLength)/1000)){
		tetraDistributionSet=c(tetraDistributionSet,as.character(sequence[[i]][j:(j+contigEvalLength)]))
	}
}
tetraDistributionSetDNAstring=DNAStringSet(tetraDistributionSet)
tetraDistributionSetCounts=(oligonucleotideFrequency(tetraDistributionSetDNAstring,width=4,with.labels=TRUE,as.prob=FALSE) + oligonucleotideFrequency(reverseComplement(tetraDistributionSetDNAstring),width=4,with.labels=TRUE,as.prob=FALSE))
tetraDistribution=(1+tetraDistributionSetCounts)/(dim(tetraDistributionSetCounts)[2]+rowSums(tetraDistributionSetCounts))
tetraDistributionPcoA=prcomp(log10(tetraDistribution),center=TRUE,scale=TRUE)

# Identify abberant tetranucleotide signatures
tetraModelP_raw=pnorm(apply(tetranucleotidePcoA$x,1,sum),mean=0,sd=sqrt(sum((tetraDistributionPcoA$sdev)^2)))
tetraModelP_twotail=tetraModelP_raw
tetraModelP_twotail[tetraModelP_twotail>0.5]=1-tetraModelP_twotail[tetraModelP_twotail>0.5]
tetraModelP=p.adjust(tetraModelP_twotail,method="BH")
#tetraModelP=tetraModelP_twotail

if(sum(tetraModelP<0.25)>4){
	coefficients_tetra=lm(tetraModelP_twotail[tetraModelP<0.25]~tetraModelP[tetraModelP<0.25])$coefficients
	BH_P_cutoff_tetra=coefficients_tetra[1]+coefficients_tetra[2]*0.05
}else{
	BH_P_cutoff_tetra=min(tetraModelP_twotail[which(tetraModelP==min(tetraModelP[tetraModelP>=(alpha/2)]))])
}

removeContigs_tetra=names(which(tetraModelP<alpha/2))
removeContigs=union(removeContigs_length,removeContigs_tetra)

allContigs=NULL
if(class(peMappingDepth)!="try-error"){
	allContigs=unique(c(allContigs,as.character(unique(peMappingDepth[,1]))))
}else{
	peMappingDepth=rbind(c("PLACEHOLDER","0"),c("PLACEHOLDER","0"))
}
if(class(mergedMappingDepth)!="try-error"){
	allContigs=unique(c(allContigs,as.character(unique(mergedMappingDepth[,1]))))
}else{
	mergedMappingDepth=rbind(c("PLACEHOLDER","0"),c("PLACEHOLDER","0"))
}
if(class(unpairedMappingDepth)!="try-error"){
	allContigs=unique(c(allContigs,as.character(unique(unpairedMappingDepth[,1]))))
}else{
	unpairedMappingDepth=rbind(c("PLACEHOLDER","0"),c("PLACEHOLDER","0"))
}
medianDepthVector=rep(0,length(allContigs))
names(medianDepthVector)=allContigs
contigLengthVector=rep(0,length(allContigs))
chunkDepthVector=NULL
chunkDepthMappingVector=NULL
removeContigs_uneven=NULL
unevenList=list()

for(j in c(1:length(allContigs))){
	myContig=allContigs[j]
	thisContigDepth_pe=NULL
	thisContigDepth_merged=NULL
	thisContigDepth_unpaired=NULL
	thisChunkDepthVector=NULL
	maxLength=width(sequence)[which(names(sequence)==allContigs[j])]
	if(sum(peMappingDepth[,1]==myContig)>0){
		thisContigDepth_pe=peMappingDepth[peMappingDepth[,1]==myContig,c(2:3)]
	}
	if(sum(mergedMappingDepth[,1]==myContig)>0){
		thisContigDepth_merged=mergedMappingDepth[mergedMappingDepth[,1]==myContig,c(2:3)]
	}
	if(sum(unpairedMappingDepth[,1]==myContig)>0){
		thisContigDepth_unpaired=unpairedMappingDepth[unpairedMappingDepth[,1]==myContig,c(2:3)]
	}
	contigPositions=c(1:maxLength)
	depthVector=rep(0,maxLength)

	depthVector[match(thisContigDepth_merged[,1],contigPositions)]=depthVector[match(thisContigDepth_merged[,1],contigPositions)]+thisContigDepth_merged[,2]
	depthVector[match(thisContigDepth_pe[,1],contigPositions)]=depthVector[match(thisContigDepth_pe[,1],contigPositions)]+thisContigDepth_pe[,2]
	depthVector[match(thisContigDepth_unpaired[,1],contigPositions)]=depthVector[match(thisContigDepth_unpaired[,1],contigPositions)]+thisContigDepth_unpaired[,2]

	medianDepthVector[j]=median(depthVector)
	contigLengthVector[j]=maxLength
	
	if(maxLength>=contigEvalLength){
		startPoints=sample(c(1:(maxLength-contigEvalLength+1)),size=(maxLength-1000)/1000)
		for(m in startPoints){
			if(m+contigEvalLength<=maxLength){
				thisChunkDepth=(depthVector[m:(m+contigEvalLength-1)])
				thisChunkDepthVector<-c(thisChunkDepthVector,median(thisChunkDepth))
			}		
		}
	}
	chunkDepthVector<-c(chunkDepthVector,thisChunkDepthVector)
	if(abs(log10((median(depthVector[1:contigEvalLength])+1)/(median(rev(depthVector)[1:contigEvalLength])+1)))>1){
		removeContigs_uneven<-c(removeContigs_uneven,myContig)
		unevenList[[myContig]]=depthVector
	}
}
removeContigs=union(removeContigs,removeContigs_uneven)

# Identify abberant coverage signatures
medianDepthModelP_raw=pnorm(medianDepthVector,mean=mean(chunkDepthVector),sd=sd(chunkDepthVector))
medianDepthModelP_twotail=medianDepthModelP_raw
medianDepthModelP_twotail[medianDepthModelP_twotail>0.5]=1-medianDepthModelP_twotail[medianDepthModelP_twotail>0.5]
medianDepthModelP=p.adjust(medianDepthModelP_twotail,method="BH")

if(sum(medianDepthModelP<0.25)>4){
	coefficients_medianDepth=lm(medianDepthModelP_twotail[medianDepthModelP<0.25]~medianDepthModelP[medianDepthModelP<0.25])$coefficients
	BH_P_cutoff_depth=coefficients_medianDepth[1]+coefficients_medianDepth[2]*0.05
}else{
	BH_P_cutoff_depth=min(medianDepthModelP_twotail[which(medianDepthModelP==min(medianDepthModelP[medianDepthModelP>=(alpha/2)]))])	
}
#medianDepthModelP=medianDepthModelP_twotail
removeContigs_coverage=names(which(medianDepthModelP<(alpha/2)))
removeContigs=union(removeContigs,removeContigs_coverage)


if(length(removeContigs)>0){
	writeXStringSet(originalSequences[-match(removeContigs,names(originalSequences))],paste(bin,".clean",sep=""))
	writeXStringSet(originalSequences[match(removeContigs,names(originalSequences))],paste(bin,".removed",sep=""))
	removeContigsOutput=NULL
	if(length(removeContigs_uneven)>0){
		write(paste("evenness-based",removeContigs_uneven,sep="\t"),file="removeContigs.txt",append=TRUE)
	}	
	if(length(removeContigs_length)>0){
		write(paste("length-based",removeContigs_length,sep="\t"),file="removeContigs.txt",append=TRUE)
	}	
	if(length(removeContigs_coverage)>0){
		write(paste("coverage-based",removeContigs_coverage,sep="\t"),file="removeContigs.txt",append=TRUE)
	}
	if(length(removeContigs_tetra)>0){
		write(paste("composition-based",removeContigs_tetra,sep="\t"),file="removeContigs.txt",append=TRUE)
	}
}else{
	writeXStringSet(sequence,paste(bin,".clean",sep=""))
	write("NO CONTIGS REMOVED",file="removeContigs.txt")
}

#PLOTTING
pdf(file=gsub(".fa",".fa.binFilter.pdf",bin))
#plotting tetranucleotides
heatmapContigColors=rep("NA",dim(tetranucleotide)[1])
heatmapContigColors[match(removeContigs_tetra,rownames(tetranucleotide))]="red"
heatmapContigColors[match(removeContigs_coverage,rownames(tetranucleotide))]="blue"
heatmapContigColors[match(intersect(removeContigs_tetra,removeContigs_coverage),rownames(tetranucleotide))]="purple"

par(cex.main=0.7)
heatmap(log10(tetranucleotide),main=paste(bin,"tetra signature"),margins=c(5,15),scale="none",RowSideColors=heatmapContigColors)
legend("topleft",c("tetra flagged","cov flagged","both"),col=c("red","blue","purple"),pch=15,cex=0.5)

#plot tetra model
tetraHist=hist(apply(tetranucleotidePcoA$x,1,sum),breaks=c(floor(min(apply(tetranucleotidePcoA$x,1,sum))):ceiling(max(apply(tetranucleotidePcoA$x,1,sum)))),col="lightblue",main="tetranucleotide distribution")
tetraSeq=seq(min(tetraHist$mids),max(tetraHist$mids),length=10000)
tetraNullModel=dnorm(seq(min(tetraHist$mids),max(tetraHist$mids),length=10000),mean=0,sd=sqrt(sum((tetraDistributionPcoA$sdev)^2)))
lowerBoundP=BH_P_cutoff_tetra
upperBoundP=1-lowerBoundP
lowerBoundTetra=qnorm(lowerBoundP,mean=0,sd=sqrt(sum((tetraDistributionPcoA$sdev)^2)))
upperBoundTetra=qnorm(upperBoundP,mean=0,sd=sqrt(sum((tetraDistributionPcoA$sdev)^2)))
polygonX=c(min(tetraSeq[intersect(which(tetraSeq>lowerBoundTetra),which(tetraSeq<upperBoundTetra))]),tetraSeq[intersect(which(tetraSeq>lowerBoundTetra),which(tetraSeq<upperBoundTetra))],max(tetraSeq[intersect(which(tetraSeq>lowerBoundTetra),which(tetraSeq<upperBoundTetra))]))
polygonY=c(0,(length(apply(tetranucleotidePcoA$x,1,sum))*tetraNullModel)[intersect(which(tetraSeq>lowerBoundTetra),which(tetraSeq<upperBoundTetra))],0)
polygon(x=polygonX,polygonY,col=alpha("green", alpha = 0.2))
points(tetraSeq,length(apply(tetranucleotidePcoA$x,1,sum))*tetraNullModel,type="l",lwd=2)


#plot depth model
depthHist=hist(medianDepthVector,breaks=c(min(c(medianDepthVector,chunkDepthVector)):max(c(medianDepthVector,chunkDepthVector))),col="lightblue",main="coverage distribution")
depthSeq=seq(min(depthHist$mids),max(depthHist$mids),length=10000)
depthNullModel=dnorm(seq(min(depthHist$mids),max(depthHist$mids),length=10000),mean=mean(chunkDepthVector),sd=sd(chunkDepthVector))
lowerBoundP=BH_P_cutoff_depth
upperBoundP=1-lowerBoundP
lowerBoundCov=qnorm(lowerBoundP,mean=mean(chunkDepthVector),sd=sd(chunkDepthVector))
upperBoundCov=qnorm(upperBoundP,mean=mean(chunkDepthVector),sd=sd(chunkDepthVector))
polygonX=c(min(depthSeq[intersect(which(depthSeq>lowerBoundCov),which(depthSeq<upperBoundCov))]),depthSeq[intersect(which(depthSeq>lowerBoundCov),which(depthSeq<upperBoundCov))],max(depthSeq[intersect(which(depthSeq>lowerBoundCov),which(depthSeq<upperBoundCov))]))
polygonY=c(0,(length(medianDepthVector)*depthNullModel)[intersect(which(depthSeq>lowerBoundCov),which(depthSeq<upperBoundCov))],0)
polygon(x=polygonX,polygonY,col=alpha("green", alpha = 0.2))
points(depthSeq,length(medianDepthVector)*depthNullModel,type="l",lwd=2)

# plot uneven contigs
if(length(unevenList)>0){
	par(mfrow=c(3,1))
	for(i in c(1:length(unevenList))){
		plot(unevenList[[names(unevenList)[i]]],type="l",main=paste("UNEVEN CONTIG ",i,"\n",names(unevenList)[i],sep=""),xlab="position (bp)",ylab="coverage")
	}
}

dev.off()
