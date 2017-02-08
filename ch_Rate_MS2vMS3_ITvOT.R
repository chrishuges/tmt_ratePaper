# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS2vMS3/Routput/")
#grab the PSM files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS2vMS3", pattern="PSMs\\.txt", full.names=TRUE)
#read in the files into a list
psmSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(psmSet) = c(sub(".*?6-2\\_(.*?)(\\_PSMs\\.txt.*|$)", "\\1", infiles))
#grab the protein files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS2vMS3", pattern="Proteins\\.txt", full.names=TRUE)
#read in the files into a list
proSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(proSet) = c(sub(".*?6-2\\_(.*?)(\\_Proteins\\.txt.*|$)", "\\1", infiles))
#grab the MS2 files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS2vMS3", pattern="MSMSSpectrumInfo\\.txt", full.names=TRUE)
#read in the files into a list
msmsSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(msmsSet) = c(sub(".*?6-2\\_(.*?)(\\_MSMSSpectrumInfo\\.txt.*|$)", "\\1", infiles))
#grab the quan files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS2vMS3", pattern="QuanSpectra\\.txt", full.names=TRUE)
#read in the files into a list
quanSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(quanSet) = c(sub(".*?6-2\\_(.*?)(\\_QuanSpectra\\.txt.*|$)", "\\1", infiles))




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

#######calculate the numbers of spectra with a full set of quan values
#create a data holder
qData = data.frame()
#loop over the files
for (i in 1:length(quanSet)){
	quanCols = c('Mass.Analyzer','Ion.Inject.Time.in.ms','First.Scan','X126','X127','X128','X129','X130','X131')
	quan = as.data.frame(quanSet[[i]])[,quanCols]
	quan$missing = rowSums(is.na(quan[,4:9]))
	qData[i,1] = names(quanSet)[i]
	qData[i,2] = mean(quan$missing,na.rm=TRUE)
}
#write out the data
write.table(qData,'ch_RawQuanNumbers.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

#######calculate the fill time metrics of the quan spectra
#loop over the files
for (x in 4:6){
	quanCols = c('Mass.Analyzer','Ion.Inject.Time.in.ms','First.Scan','X126','X127','X128','X129','X130','X131')
	quan = as.data.frame(quanSet[[x]])[,quanCols]
	message(nrow(quan))
	message(nrow(quan[quan$Ion.Inject.Time.in.ms==250,]))
	message(mean(quan$Ion.Inject.Time.in.ms,na.rm=TRUE))
}


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
	#subset out PSMs to keep only those that provide a full set of reporter ions
	#pep = subset(pep, rowSums(is.na(pep[,3:8]))<1)
	pep = filter(pep, Master.Protein.Accessions == 'SPIKE')
	pep$totSignal = rowSums(pep[,3:8],na.rm=TRUE)
	#keep only spectra with a majority of reporter ions
	pep = subset(pep, rowSums(is.na(pep[,3:8]))<2)
	#get the reporter ion ratios
	pep[,3:8] = apply(pep[,3:8],2,function(y) (y/pep$totSignal)*12)
	#get the percent errors
	expRatios = c(3,3,3,1,1,1)
	for (l in 3:8){
		pep[,l] = abs((pep[,l] - expRatios[l-2])/expRatios[l-2]) * 100
	}
	#data output
	spikeData[[i]] = pep
}
	

#bind all of the data together into a single frame
lsetout = data.frame()
for (x in 1:3){
	lset = as.data.frame(spikeData[[x]])
	lsetout = rbind(lsetout,lset)
}

colMeans(lsetout[,3:8],na.rm=TRUE)
mean(colMeans(lsetout[,3:8],na.rm=TRUE))

#####plot the data
#open the plot
pdf('ch_Rate_MS2vMS3_IT-MS3_SpikeQuant.pdf')
#choose your data set to plot
#set the colors
cols1 = brewer.pal(9,'RdBu')[7]
lcols = brewer.pal(9,'RdBu')[1]
#plot the data
boxplot(lsetout[,c(3:8)],
		border='black',
		col = cols1,
		outline = FALSE,
		boxlwd = 3,
		ylim = c(0,6),
		boxwex=0.5,
		staplelwd=2,
		whisklwd=2)
abline(h=3,col=lcols,lty=2)
abline(h=1,col=lcols,lty=2)
dev.off()
