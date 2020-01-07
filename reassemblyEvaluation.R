myData=read.table(commandArgs()[3],sep=",",header=TRUE,comment.char="")
rounds=as.numeric(gsub(".*[.]","",gsub(".fa$","",myData[,1])))

correlationValues=NULL
for(i in c(3:dim(myData)[2])){
      thisCorrelation=cor.test(rounds,myData[,i])
      correlationValues=rbind(correlationValues,c(thisCorrelation$statistic,thisCorrelation$p.value))
}
correlationValues<-cbind(colnames(myData)[3:dim(myData)[2]],correlationValues)

#evaluate whether to continue
CONCLUSION="GO"
contaminationCorrelation=cor.test(rounds,myData$Contamination)
if(!is.na(contaminationCorrelation$statistic)){
	if((contaminationCorrelation$statistic>0)&&(contaminationCorrelation$p.value<0.05)){
		if(max(myData$Contamination)>10){
			CONCLUSION="STOP"
		}
	}
}
if(min(sort(as.numeric(correlationValues[,3])))>0.05){
	CONCLUSION="STOP"
}
correlationValues<-rbind(c("Value","Correlation","P-value"),correlationValues)
write.table(correlationValues,file="reassemblyEvaluationStats.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


cat(paste(CONCLUSION,"\n",sep=""))
