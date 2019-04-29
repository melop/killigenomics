setwd("/beegfs/group_dv/home/RCui/killifish_genomes/human_alfred/1kG_GRCH37/allelefreq_above0.05/consurf/plotallpops/");
library(ggplot2);
library(lme4)
library(lmerTest)


datDatDesc <- read.table("popdesc.txt", header = T);
arrYLims <- c( 4.0 , 6)
bExcludeX <- T;
bExcludeY <- T;
#sOutStem <- paste( "allpops", sep=".");
sOutStem <- paste( "allpops_randeffect", sep=".");
bOnlyX <- F;

#bOnlyX <- T;

#sOutStem <- paste( "allpops_chrX", sep=".");

sink(paste(sOutStem,".txt",sep=""));
pdf(paste(sOutStem,".pdf",sep=""), width=10, height=3);

datNonSynAllFreq <- NULL;

for(nRow in 1:nrow(datDatDesc)) {
  sF <- paste(datDatDesc$dir[nRow], "/", datDatDesc$pop[nRow],".nonsyn.AF.Consurf.txt", sep="");
  datNonSynAllFreq2 <- read.table(sF, header=T, sep="\t", stringsAsFactors = F);
  datNonSynAllFreq2$Pop <- datDatDesc$pop[nRow]
  datNonSynAllFreq2$Genus <- datDatDesc$genus[nRow]
  
  datNonSynAllFreq <- rbind(datNonSynAllFreq, datNonSynAllFreq2);
}

if (bExcludeX & (!bOnlyX) ) {
	datNonSynAllFreq <- datNonSynAllFreq[!grepl("X:", datNonSynAllFreq$coord , fixed=T), ];
}
if (bExcludeY) {
	datNonSynAllFreq <- datNonSynAllFreq[!grepl("Y:", datNonSynAllFreq$coord , fixed=T), ];
}

if (bOnlyX) {
	datNonSynAllFreq <- datNonSynAllFreq[grepl("X:", datNonSynAllFreq$coord , fixed=T), ];
}

datNonSynAllFreq$Pop <- factor( datNonSynAllFreq$Pop , levels = datDatDesc$pop)
datNonSynAllFreq$Genus <- factor( datNonSynAllFreq$Genus )

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
 # arrDry <- datNonSynAllFreq[datNonSynAllFreq$AFBin == nAF & datNonSynAllFreq$Genus == 'pan' , 'ConsurfScore'  ]
 # arrWet <- datNonSynAllFreq[datNonSynAllFreq$AFBin == nAF & datNonSynAllFreq$Genus == 'homo' , 'ConsurfScore'  ]
 # print(wilcox.test(arrDry , arrWet))
  cat("====================== Allele Freq Bin: ", nAF, "=================\n");
  datThis <- datNonSynAllFreq[datNonSynAllFreq$AFBin == nAF, ];
  #print(summary(lm(ConsurfScore ~ Genus/Pop , data=datThis) ));
	model = lmer(ConsurfScore ~ Genus + (1|Pop),
            data=datThis,
            REML=TRUE)
	
	print(anova(model));
	print(summary(model));
}

dsum <- data_summary(datNonSynAllFreq, 'ConsurfScore'  , c('AFBin', 'Pop') );
dsum$Pop <- factor(dsum$Pop , levels = as.character(datDatDesc$pop)  );
dsum
#summary(lm(ConsurfScore ~ derived_freq * (Genus / Pop) , data=datNonSynAllFreq));
#summary(lm(datNonSynAllFreq$ConsurfScore ~ datNonSynAllFreq$derived_freq * datNonSynAllFreq$Genus  ));
model = lmer(ConsurfScore ~ derived_freq * Genus + (1|Pop),
            data=datNonSynAllFreq,
            REML=TRUE)
print(anova(model));
print(summary(model));

fnCol1 <- colorRampPalette(c("#F21b3f", "#a00a1a"));
fnCol2 <- colorRampPalette(c("#FE9A2E", "#8A4B08"));
fnCol3 <- colorRampPalette(c("#FACC2E", "#B18904"));
fnCol4 <- colorRampPalette(c("#C8FE2E", "#688A08"));
fnCol5 <- colorRampPalette(c("#D358F7", "#6A0888"));
fnCol6 <- colorRampPalette(c("#08bdbd", "#006666"));

arrTable <- table(datDatDesc$superpop)
arrColCounts <- c(arrTable[['chimp']], arrTable[['EUR']], arrTable[['EAS']], arrTable[['AMR']], arrTable[['SAS']], arrTable[['AFR']]);

p <- ggplot(dsum, aes(x=AFBin, y=ConsurfScore, fill=Pop)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=ConsurfScore-sem, ymax=ConsurfScore+sem), width=.1, size=0.4, alpha=0.6,color="#666666", 
                position=position_dodge(.9)) + coord_cartesian(ylim = arrYLims )

p + scale_fill_manual(values=c(fnCol1(arrColCounts[1]) , fnCol2(arrColCounts[2]) , fnCol3(arrColCounts[3]),  fnCol4(arrColCounts[4]), fnCol5(arrColCounts[5]) , fnCol6(arrColCounts[6]) )) + theme_classic() + theme(legend.text=element_text(size=5), legend.key.size = unit(0.5,"line"));


# boxplot(datNonSynAllFreq$ConsurfScore ~ datNonSynAllFreq$AFBin * datNonSynAllFreq$Pop );

# ggplot2.boxplot(data = datNonSynAllFreq,  addMean=TRUE, meanPointShape=23, meanPointSize=3,
#                    meanPointColor="black", meanPointFill="blue", groupName = 'AFBin' , xName = 'Pop', yName = 'ConsurfScore')
# means <- tapply(datNonSynAllFreq$ConsurfScore,datNonSynAllFreq$AFBin,mean)
# points(means,col="red",pch=18)
sink();
dev.off();
