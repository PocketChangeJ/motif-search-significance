library("abind")
library(gplots)

########################################################################
#################################### FUNCTIONS
loadMatrix = function(file=file){
	data=read.table(file, sep="\t", header=T)

	### Some of the name of the motif are not unique
	### As they are IDs (as row names), we have to make them unique
	rownames(data) = make.names(data[,1], unique=TRUE)
	data=data[,2:dim(data)[2]]
	return(data)
}

######### Co-occurence analysis

loadMatrices = function(path=path, motif=motif, pattern=pattern, prefix=prefix){

	### Searching for FIMO directories of a given pattern
	dirs = dir(path = path, pattern = pattern)

	matrices=data.frame()

	### For all FIMO directories found...
	for (dir in dirs){
		### Get the names of the files of a given pattern
		files = list.files(path = paste( path,"/", dir, sep="" ), pattern = prefix)
		for (file in files){

			### Reading the files
			test = loadMatrix( file=paste( path,"/", dir, "/", file, sep="" ) )

			### Putting the data into a table
			if(dim(matrices)[1] == 0){
				matrices = test
			}else{
				matrices = abind(test, matrices)
			}
		}
	}

	### Transforming the resulting table into a multi-dimentional array
	### There is one table per file
	if(length(dirs) > 1){
		matrices.array = array(matrices, dim=c(dim(matrices)[1], dim(matrices)[1], length(dirs) ) )
		dimnames(matrices.array) <- list(rownames(matrices), rownames(matrices) )
		matrices=matrices.array
	}
	return(matrices)
} 


loadMatricesCoOccurenceControl=function(path=path, motif=motif){
	pattern= paste("FIMO-", motif, ".*_[01213456789]+$", sep="")
	loadMatrices(path=path, motif=motif, prefix="CountCoOccurence", pattern=pattern)
}

loadMatricesCoOccurenceMain=function(path=path, motif=motif, name=name){
	pattern= paste("FIMO-", motif, "-", name, "$", sep="")
	loadMatrices(path=path, motif=motif, prefix="CountCoOccurence", pattern=pattern)
}

### Function computing Z-score
zScore = function(interest=interest, control=control, title=title){

	### zscore = (mean of data - expected mean)/standard deviation

	### control is a multi-dimentional matrix
	### Computing mean
	control.mean = apply(control,c(1,2),mean)
	### Computing sd
	control.sd = apply(control,c(1,2),sd)
	### Computing Z-score
	zscore = (interest-control.mean)/control.sd

	### removing columns with only NaN
	zscore.2 = zscore[ , ! apply( zscore , 2 , function(x) all(is.na(x)) ) ]
	### removing rows with only NaN
	zscore.3 = zscore.2[ ! apply( zscore.2 , 1 , function(x) all(is.na(x)) ) , ]

	### Changing Nan into 0
	zscore.3[!is.finite(as.matrix(zscore.3))] = 0

	#print(zscore.3)

	colors=c(seq(min(zscore.3), -1, length=50), seq(-1, 1, length=50), seq(1, max(zscore.3), length=50) )
	my_palette = colorRampPalette(c("blue","grey","red"))(149)

	#x11()
	dev.new(width=10, height=9)
	heatmap.2(as.matrix(zscore.3), col=my_palette, breaks=colors, symbreaks=F, symkey=T, main=title,
		scale="none",cex.axis=0.4, margins=c(8,8), trace="none" )

	if(title != ""){
		pdf(paste(title,"-CoOc.pdf", sep="") , pointsize=12)
		par(cex.main=0.7)
		heatmap.2(as.matrix(zscore.3), col=my_palette, breaks=colors, symbreaks=F, symkey=T, main=title,
			scale="none", cexRow=0.2, cexCol=0.2, margins=c(8,8), trace="none" )
		dev.off()
	}

	return(zscore.3)
}

######### Occurence analysis

loadMatricesOc = function(path=path, motif=motif, pattern=pattern, prefix=prefix){

	### Searching for FIMO directories of a given pattern
	dirs = dir(path = path, pattern = pattern)

	matrices=data.frame()

	### For all FIMO directories found...
	for (dir in dirs){
		### Get the names of the files of a given pattern
		files = list.files(path = paste( path,"/", dir, sep="" ), pattern = prefix)

		### Reading the files
		for (file in files){
			test = loadMatrix( file=paste( path,"/", dir, "/", file, sep="" ) )

			### Putting the data into a table
			if(dim(matrices)[1] == 0){
				matrices = data.frame( apply(test, 2, sum) )
				colnames(matrices)[1] = file
			}else{
				test = apply(test, 2, sum)
				matrices = abind(matrices, test)
				colnames(matrices)[dim(matrices)[2]] = file
			}
		}
	}
	return(t(matrices))
} 


loadMatricesOccurenceControl=function(path=path, motif=motif){
	pattern= paste("FIMO-", motif, ".*_[01213456789]+$", sep="")
	loadMatricesOc(path=path, motif=motif, prefix="CountOccurence", pattern=pattern)
}

loadMatricesOccurenceMain=function(path=path, motif=motif, name=name){
	pattern= paste("FIMO-", motif, "-", name, "$", sep="")
	loadMatricesOc(path=path, motif=motif, prefix="CountOccurence", pattern=pattern)
}

zScoreOc = function(interest=interest, control=control, title=""){

	### zscore = (mean of data - expected mean)/standard deviation

	### control is a multi-dimentional matrix
	### Computing mean
	control.mean = apply(control,2,mean)
	### Computing sd
	control.sd = apply(control,2,sd)
	### Computing Z-score
	zscore = (interest-control.mean)/control.sd

	### Changing Nan into 0
	#zscore[!is.finite(as.matrix(zscore))] = 0

	zscore = zscore[,order(zscore)]

	#x11()
	dev.new(width=10, height=9)
	par(mar=c(5.1, 10 ,4.1 ,2.1))
	barplot(zscore[! is.na(zscore)], horiz=TRUE, col="blue", main=title,
		names.arg=colnames(zscore)[! is.na(zscore)], las=1, cex.names=.5)

	if(title != ""){
		pdf(paste(title,"-Oc.pdf", sep=""))
		par(mar=c(5.1, 10 ,4.1 ,2.1))
		barplot(zscore[! is.na(zscore)], horiz=TRUE, col="blue", main=title,
			names.arg=colnames(zscore)[! is.na(zscore)], las=1, cex.names=.2)
		dev.off()
	}

	return(zscore)
}

#### Computing pvalue with Bonferroni & hochberg correction

pvalueBH = function(zscore.Oc = zscore, main.Oc= main, control.Oc= control, error=0.05){
	#######################################################################
	#### Computing the p-value from the Z-score from the Occurrence matrix
	f <- function(x){ pnorm(x, mean=0, sd=1) }
	f.Rev <- function(x){ 1-pnorm(x, mean=0, sd=1) }

	#### Computing p-values
	pvalue = unlist( lapply(zscore.Oc[which(zscore.Oc<0)], f) )
	pvalue = c(pvalue, unlist( lapply(zscore.Oc[which(zscore.Oc>0)], f.Rev) ) )
	pvalue = 2*pvalue

	#### Computing Benjamini & Hochberg adjusted pvalues
	## error fixed to 0.05 (bu default)
	BH.data = data.frame( pvalue.sorted=sort(pvalue ), indice=seq(1,length(pvalue),1) )
	BH = BH.data$indice*error/length(pvalue)

	BH.data = cbind(BH.data, BH)

	#### If pvalue < Bonferoni & Hochberg then p-value is significant
	sig = rep(0,length(pvalue))
	sig[which(BH.data$pvalue.sorted < BH.data$BH)] = 1

	BH.data = cbind(BH.data, sig)
	name=rownames(BH.data)

	#### Adding zscore (occurrences) and motif occurrence counts (sample + controls)
	zscore.Oc.sig = t( as.data.frame(t(zscore.Oc))[, name ] )
	main.Oc.sig = t( as.data.frame(main.Oc)[, name ] )
	colnames(main.Oc.sig) = 'occurrence'
	control.Oc.sig = t( as.data.frame(control.Oc)[, name ] )
	colnames(control.Oc.sig) = paste( "C", seq(1, dim(control.Oc.sig)[2]), sep="" )

	BH.data = cbind(BH.data, zscore=zscore.Oc.sig )
	BH.data = cbind(BH.data, occurrence=main.Oc.sig )
	BH.data = cbind(BH.data, control.Oc.sig )

	return(BH.data)

}

#### Plot data 
plotDataCluster = function(data=data, title=title){
	expected = apply(data[,7:dim(data)[2]], 1, median)
	observed = data$occurrence
	names(observed) = rownames(data)

	standard.error = function(x){ 1.96*(sd(x)/sqrt(length(x)) ) }
	se = apply(BH.selec[,7:dim(BH.selec)[2]], 1, standard.error)

	data = data.frame(
		cluster=c( names(observed), names(expected) ),
		occurrence=c(observed, expected), 
		type=c( rep("observed",length(expected)) , rep("expected",length(observed)) ),
		se = c(rep(NA,length(se)), se  )
		)

	data$type = factor(data$type, levels = rev(levels(data$type)))

	library(ggplot2)

	g <- ggplot(data, aes(x=cluster, y=occurrence, fill=type ) ) + 
		geom_bar(position=position_dodge(), stat="identity") +
	    geom_bar(position=position_dodge(), stat="identity", colour="black", show_guide=FALSE) +
	    scale_fill_manual(values=c("#000000", "#999999")) + 
	    geom_errorbar( aes(ymin=occurrence-se, ymax=occurrence+se), position=position_dodge(width=0.9), width=0.25 ) +
	    ggtitle(title)

	return(g)

}

#### Plot data motif 
plotDataMotif = function(data=data, title=title){
	expected = apply(data[,7:dim(data)[2]], 1, mean)
	observed = data$occurrence
	pvalue = data$pvalue.sorted
	names(observed) = rownames(data)

	standard.error = function(x){ 1.96*(sd(x)/sqrt(length(x)) ) }
	se = apply(data[,7:dim(data)[2]], 1, standard.error)

	data2 = data.frame(
		motifs=c( names(observed), names(expected) ),
		occurrence=c(observed, expected), 
		type=c( rep("observed",length(expected)) , rep("expected",length(observed)) ),
		se = c(rep(NA,length(se)), se  ),
		p = c(format(data$pvalue.sorted, scientific = TRUE, digits = 3), rep(NA, length(data$pvalue.sorted)) ),
		x = c(data$indice, rep(NA,length(data$indice)) ),
		y = c(apply(data.frame(expected, observed), 1, max), rep(NA,length(data$indice)) )
		)

	#data2$motifs = factor(data2$motifs, levels = c("GATA", "REST", "NR1H2_RXRA", "HNF4A") )
	#data2$motifs = factor(data2$motifs, levels = c("Zfx", "NFYA", "SP1", "Klf4", "CTCF", "REST","RREB1","Myf") )
	#data2$motifs = factor(data2$motifs, levels = c("SP1", "Klf4", "Zfx", "CTCF","Pax4") )
	data2$type = factor(data2$type, levels = rev(levels(data2$type)))
	data2$motifs = factor(data2$motifs, levels=data2[order(pvalue),"motifs"] )

	library(ggplot2)

	col=c("#000000", "#999999")

	g <- ggplot(data=data2, aes(x=motifs, y=occurrence, fill=type ) ) + 
		geom_bar(position=position_dodge(), stat="identity") +
	    geom_bar(position=position_dodge(), stat="identity", colour="black", show_guide=FALSE) +
	    scale_fill_manual(values=col) + 
	    geom_errorbar( aes(ymin=occurrence-se, ymax=occurrence+se), position=position_dodge(width=0.9), width=0.25 ) +
	    geom_text(aes(y=y, x=x, label=p ), position= position_dodge(width=0.9), vjust=-.5, color="black") +
	    ggtitle(title)

	return(g)

}


