# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS3nce/Routput/")
#grab the MS2 files from Proteome Discoverer
infiles = dir("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_MS3nce", pattern="PSMs\\.txt", full.names=TRUE)
#read in the files into a list
psmSet <- lapply(infiles, read.table, header=TRUE, sep='\t')
names(psmSet) = c(sub(".*?6-2\\_(.*?)(\\_PSMs\\.txt.*|$)", "\\1", infiles))

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
	pep$totSignal = rowMeans(pep[,3:8],na.rm=TRUE)
	#output the data
	psmData[i,1] = names(psmSet)[i]
	psmData[i,2] = nrow(pep)
	psmData[i,3] = mean(pep$totSignal,na.rm=TRUE)
	
}
#write out the data
write.table(psmData,'ch_RawPSMNumbers.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
















