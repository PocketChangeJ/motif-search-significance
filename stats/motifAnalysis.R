source("motifAnalysisLib.R")

########################################################################
#################################### MAIN
#setwd("")
sample="" # <- enter the prefix of your sample

##### Analysis on ALL VERTEBRATE motifs
### Co-Occurence analysis
controls = loadMatricesCoOccurenceControl(path="Controls", motif="vertebrate")
main = loadMatricesCoOccurenceMain(path="Main", motif="vertebrate", name=sample)
zscore = zScore(interest=main, control=controls, title=paste(sample, "-vertebrate",sep="") )

### Occurence analysis
controls.Oc = loadMatricesOccurenceControl(path="Controls", motif="vertebrate")
main.Oc = loadMatricesOccurenceMain(path="Main", motif="vertebrate", name=sample)
zscore.Oc = zScoreOc(interest=main.Oc, control=controls.Oc, title=paste(sample, "-vertebrate",sep="") )

dev.off()
dev.off()

########################### Occurence - motif analysis
BH.data = pvalueBH(zscore.Oc = zscore.Oc, main.Oc= main.Oc, control.Oc= controls.Oc, error=0.05)
BH.selec = BH.data[which(BH.data$sig==1),]
BH.selec = BH.selec[order(BH.selec$pvalue.sorted, -BH.selec$zscore),]
BH.selec$indice = seq(1,dim(BH.selec)[1],1)

write.table(BH.data, file="occurenceTable.tsv", sep="\t", col.names=T, row.names=T, quote=F)
write.table(BH.selec, file="occurenceTable_sig.tsv", sep="\t", col.names=T, row.names=T, quote=F)

pdf(paste("Occurrence_",sample,".pdf", sep=""), width=12, height=7)
plotDataMotif(data=BH.selec[1:10,], title="")
dev.off()

########################### Co-Occurence - motif analysis
#### Plot selected motifs

name=rownames(BH.selec)
selec = zscore[name,]
selec = selec[, which( abs( apply(selec,2,mean) ) >= 0.5 ) ]

colors=c(seq(min(selec), -1, length=50), seq(-1, 1, length=50), seq(1, max(selec), length=50) )
my_palette = colorRampPalette(c("blue","grey","red"))(149)

pdf(paste("Zoom_",sample,".pdf", sep=""),width=11, height=11)
heatmap.2(as.matrix(t(selec)), trace="none",scale="none", dendrogram="none", density.info="none",
	col=my_palette, breaks=colors,cexRow=.7, cexCol=.6, srtCol=45, na.rm=TRUE, margins=c(6,0), keysize=1,
	lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(0.3, 4, 1 ), lhei=c(0.8, 4))
dev.off()

