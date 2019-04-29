setwd("/beegfs/group_dv/home/RCui/killifish_genomes/human_alfred/1kG_GRCH37/allelefreq_above0.05/consurf");
library(ggplot2);
sWet <- "YRI.nonsyn.AF.Consurf.txt"
# sWet <- "CEU.nonsyn.AF.Consurf.txt"
sDry <- "JPT.nonsyn.AF.Consurf.txt"
# sWet <- "CHB.nonsyn.AF.Consurf.txt"
# 
# sDry <- "CEU.nonsyn.AF.Consurf.txt"
#sDry <- "../../../chimp/vcfs/hg19/allelefreq_above0.05/consurf/ptt.nonsyn.AF.Consurf.txt"
#sDry <- "../../../chimp/vcfs/hg19/allelefreq_above0.05/consurf/ptv.nonsyn.AF.Consurf.txt"
#sDry <- "../../../chimp/vcfs/hg19/allelefreq_above0.05/consurf/pte.nonsyn.AF.Consurf.txt"

arrYLims <- c( 4.0 , 6)


sOutStem <- paste(basename(sWet), basename(sDry), sep=".");
sink(paste(sOutStem,".txt",sep=""));
pdf(paste(sOutStem,".pdf",sep=""), width=5, height=3);

datNonSynAllFreq <- read.table(sDry, header=T, sep="\t");
datNonSynAllFreq$Pop <- 'dry'
arrCoords1 <- datNonSynAllFreq$coord;

datNonSynAllFreq2 <- read.table(sWet, header=T, sep="\t");
datNonSynAllFreq2$Pop <- 'wet'
arrCoords2 <- datNonSynAllFreq2$coord;

#filter out shared polymorphisms or divergences:
#datNonSynAllFreq <- datNonSynAllFreq[ ( (! (datNonSynAllFreq$coord %in% arrCoords2)) & datNonSynAllFreq$derived_freq >= 0.95) | datNonSynAllFreq$derived_freq < 0.95  , ]
#datNonSynAllFreq2 <- datNonSynAllFreq2[ ( (! (datNonSynAllFreq2$coord %in% arrCoords1)) & datNonSynAllFreq2$derived_freq >= 0.95) | datNonSynAllFreq2$derived_freq < 0.95  , ]


datNonSynAllFreq <- rbind(datNonSynAllFreq, datNonSynAllFreq2);

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]][!is.na(x[[col]])] )))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


arrBreaks <- seq(0.0, 1, 0.2);
datNonSynAllFreq$AFBin <- cut(datNonSynAllFreq$derived_freq, breaks = arrBreaks);

#do stats:
for(nAF in levels(datNonSynAllFreq$AFBin) ) {
  arrDry <- datNonSynAllFreq[datNonSynAllFreq$AFBin == nAF & datNonSynAllFreq$Pop == 'dry' , 'ConsurfScore'  ]
  arrWet <- datNonSynAllFreq[datNonSynAllFreq$AFBin == nAF & datNonSynAllFreq$Pop == 'wet' , 'ConsurfScore'  ]
  print(wilcox.test(arrDry , arrWet))
}

dsum <- data_summary(datNonSynAllFreq, 'ConsurfScore'  , c('AFBin', 'Pop') );
dsum
summary(lm(datNonSynAllFreq$ConsurfScore ~ datNonSynAllFreq$derived_freq * datNonSynAllFreq$Pop));


p <- ggplot(dsum, aes(x=AFBin, y=ConsurfScore, fill=Pop)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=ConsurfScore-sem, ymax=ConsurfScore+sem), width=.2,
                position=position_dodge(.9)) + coord_cartesian(ylim = arrYLims )

p + scale_fill_manual(values=c("#F21b3f", "#08bdbd")) + theme_classic()


# boxplot(datNonSynAllFreq$ConsurfScore ~ datNonSynAllFreq$AFBin * datNonSynAllFreq$Pop );

# ggplot2.boxplot(data = datNonSynAllFreq,  addMean=TRUE, meanPointShape=23, meanPointSize=3,
#                    meanPointColor="black", meanPointFill="blue", groupName = 'AFBin' , xName = 'Pop', yName = 'ConsurfScore')
# means <- tapply(datNonSynAllFreq$ConsurfScore,datNonSynAllFreq$AFBin,mean)
# points(means,col="red",pch=18)
sink();
dev.off();