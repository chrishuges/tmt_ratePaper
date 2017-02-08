# TODO: Add comment
# 
# Author: cshughes
###############################################################################

#read in the chromatogram files
setwd('/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_chromatograms/Routput')
chr<-read.table("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_chromatograms/ch_SM-liu_chromatogram.txt", header=TRUE, sep="\t")
#chr<-read.table("/Users/cshughes/Documents/projects/tmtPlexing/ratePaper_chromatograms/ch_19Jan2017_SM-1-1_Sp-6-2_OT-MS3_CID35-HCD60_1_chromatogram.txt", header=TRUE, sep="\t")


#make the plot
pdf('ch_Rate_Liu_chromatogram.pdf')
plot(chr[chr$Time>5 & chr$Time<120,], 
		type='l',
		col = brewer.pal(9,'RdBu')[9]
)
dev.off()


