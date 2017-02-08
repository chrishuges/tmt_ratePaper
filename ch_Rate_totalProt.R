# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_totalProt/Routput/")
#grab the PSM files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_totalProt", pattern="PSMs\\.txt", full.names=TRUE)
#read in the files into a list
psmSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(psmSet) = c(sub(".*?6-2\\_(.*?)(\\_PSMs\\.txt.*|$)", "\\1", infiles))
#grab the protein files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_totalProt", pattern="Proteins\\.txt", full.names=TRUE)
#read in the files into a list
proSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(proSet) = c(sub(".*?6-2\\_(.*?)(\\_Proteins\\.txt.*|$)", "\\1", infiles))
#grab the MS2 files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_totalProt", pattern="MSMSSpectrumInfo\\.txt", full.names=TRUE)
#read in the files into a list
msmsSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(msmsSet) = c(sub(".*?6-2\\_(.*?)(\\_MSMSSpectrumInfo\\.txt.*|$)", "\\1", infiles))

#######calculate the mean number of MS2 spectra
#create a data holder
ms2Data = data.frame()
#loop over the files
for (i in 1:length(msmsSet)){
	spec = as.data.frame(msmsSet[[i]])
	ms2Data[i,1] = names(msmsSet)[i]
	ms2Data[i,2] = nrow(spec)	
}
#write out the data
write.table(ms2Data,'ch_RawMS2Numbers.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


#######calculate the number of PSMs at 'High' confidence
#create a data holder
psmData = data.frame()
#loop over the files
for (i in 1:length(psmSet)){
	pepCols = c('Annotated.Sequence','Master.Protein.Accessions','X126','X127','X128','X129','X130','X131')
	pep = as.data.frame(psmSet[[i]])[,pepCols]
	#subset out PSMs to keep only those that provide a full set of reporter ions
	#pep = subset(pep, rowSums(is.na(pep[,3:8]))<1)
	pep = filter(pep, Master.Protein.Accessions != 'sp')
	pep = unique(pep)
	#output the data
	psmData[i,1] = names(psmSet)[i]
	psmData[i,2] = nrow(pep)	
	
}
#write out the data
write.table(psmData,'ch_RawPSMNumbers.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

#######calculate the number of Proteins at 'High' confidence
#create a data holder
proData = data.frame()
#loop over the files
for (i in 1:length(proSet)){
	proCols = c('Accession','Master','Number.of.Unique.Peptides')
	pro = as.data.frame(proSet[[i]])[,proCols]
	#subset out PSMs to keep only those that provide a full set of reporter ions
	#pep = subset(pep, rowSums(is.na(pep[,3:8]))<1)
	pro = filter(pro, Accession != 'sp')
	pro = unique(pro)
	#output the data
	proData[i,1] = names(proSet)[i]
	proData[i,2] = nrow(pro)	
	
}
#write out the data
write.table(proData,'ch_RawProNumbers.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)



#######calculate the variance of the spiked peptides
#create a data holder
spikeData = list()
#loop over the files
for (i in 1:length(psmSet)){
	pepCols = c('Annotated.Sequence','Master.Protein.Accessions','X126','X127','X128','X129','X130','X131')
	pep = as.data.frame(psmSet[[i]])[,pepCols]
	#subset out PSMs to keep only those that provide a majority set of reporter ions
	pep = subset(pep, rowSums(is.na(pep[,3:8]))<3)
	pep = filter(pep, Master.Protein.Accessions == 'SPIKE')
	pep$Annotated.Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep$Annotated.Sequence))
	pep = unique(pep)
	pep$meanHigh = rowMeans(pep[,3:5],na.rm=TRUE)
	pep$meanLow = rowMeans(pep[,6:8],na.rm=TRUE)
	pep$logFC = log2(pep$meanHigh/pep$meanLow)
	#data output
	spikeData[[i]] = pep
}

#plot the two data sets
ot = as.data.frame(spikeData[[1]])[,c(1,11)]
it = as.data.frame(spikeData[[2]])[,c(1,11)]
ot = aggregate(logFC~Annotated.Sequence,data=ot,na.action=na.pass,FUN=mean,na.rm=TRUE)
it = aggregate(logFC~Annotated.Sequence,data=it,na.action=na.pass,FUN=mean,na.rm=TRUE)
#combine the data
ot.it = merge(ot,it,by='Annotated.Sequence')


#xCol = col2rgb(brewer.pal(9,'Greens')[7])
xCol = col2rgb('grey40')
#make the plot
pdf('ch_Rate_SPIKE-FC.pdf')
plot(ot.it$logFC.x, 
		ot.it$logFC.y, 
		xlab = 'log10(TMT126 Intensity)',
		#col = rgb(xCol[1,],xCol[2,],xCol[3,],75,maxColorValue=255),
		col = rgb(xCol[1,],xCol[2,],xCol[3,],75,maxColorValue=255),
		pch = 20,
		cex = 1.5,
		xlim = c(-1,4),
		ylim = c(-1,4))
dev.off()










