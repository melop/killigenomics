setwd("/beegfs/group_dv/home/RCui/killifish_genomes/codonbias/46spp/codon_usage");
library(lme4);
require(lmerTest)


sCAI <- "CAI_min_100AA.txt";
#sCAI <- "CAI_min_100AA_frameshift_1_highexpgene_ref.txt";
#sCAI <- "CAI_min_100AA_frameshift_2_highexpgene_ref.txt";
#sCAI <- "CAI_min_100AA_frameshift_1_revcomp_highexpgene_ref.txt";
#sCAI <- "CAI_min_100AA_frameshift_2_revcomp_highexpgene_ref.txt";
#sCAI <- "CAI_min_100AA_revcomp_highexpgene_ref.txt";


sOut <- "CAI_correlate_exp";
#sOut <- "CAI_correlate_exp_frameshift_1";
#sOut <- "CAI_correlate_exp_frameshift_2";
#sOut <- "CAI_correlate_exp_frameshift_1_revcomp";
#sOut <- "CAI_correlate_exp_frameshift_2_revcomp";
#sOut <- "CAI_correlate_exp_revcomp";


sExp <- "~/killifish_genomes/codonbias/gene_expression_levels/all.expression.level.txt";
sWhichTissue <- "AvgRPKM";
#sWhichTissue <- "NOR";

datCAI <- read.table(sCAI , header=T, stringsAsFactors = F, na.strings="NA");
datExp <- read.table(sExp , header= T, sep="\t", quote = '', stringsAsFactors = F, na.strings="NA");
datExp <- datExp[datExp$SdOverAvg < 0.5, ];

#datCAI <- datCAI[complete.cases(datCAI),];

arrColnames <- colnames(datCAI);
arrSpecies <- arrColnames[seq(2, ncol(datCAI), 2)];
arrLengths <- arrColnames[seq(3, ncol(datCAI), 2)];


datMerged <- merge(datCAI, datExp,  by.x='GeneName', by.y='OrthoID');

datFilter <- datMerged[datMerged[,sWhichTissue] > 1, ];

sink(file=paste(sOut , ".summary.txt", sep=""))

pdf(paste(sOut , ".pdf", sep=""));
par(mar=c(5.1,4.1,4.1,2.1));

for(sSp in arrSpecies) {
cat("doing " , sSp , "\n");
fitLm <- lm(datFilter[, sSp] ~ log(datFilter[, sWhichTissue])) 
print(summary(fitLm));
plot( log(datFilter[, sWhichTissue]) , datFilter[, sSp] , pch=19, cex=0.5, xlab="Log(Gene expression level)", ylab=sSp, col=rgb(1/2,0,1,1/4));
abline(a=fitLm$coefficients[1], b=fitLm$coefficients[2], col="red", lwd=3);
}


arrExtremes <- quantile(datFilter[, sWhichTissue] , c(0.05,0.95));
#prepare dataframe for boxplot:

datTopQuantile <- datFilter[datFilter[, sWhichTissue]>=arrExtremes[2], ];
datBottomQuantile <- datFilter[datFilter[, sWhichTissue]>=arrExtremes[1], ];
datBoxplot <-data.frame();
datAnalysis <- data.frame();

for(nSp in 1:length(arrSpecies)) {
  sSp <- arrSpecies[nSp];
  sLen <- arrLengths[nSp];
  datBoxplot <- rbind.data.frame( datBoxplot, cbind( as.numeric( datTopQuantile[, sSp]) , "highexp", sSp ));
  datBoxplot <- rbind.data.frame( datBoxplot, cbind( as.numeric(datBottomQuantile[, sSp]) , "lowexp", sSp ));
  datAnalysis <- rbind.data.frame( datAnalysis, cbind( datFilter$GeneName, datFilter[, sSp] , datFilter[, sWhichTissue], sSp , datFilter[, sLen] ));
}
colnames(datBoxplot) <- c("CAI", "Expression" , "Species");
colnames(datAnalysis) <- c("GeneName", "CAI", "Expression" , "Species", "CodonLen");
datBoxplot$CAI <- as.numeric(as.character( datBoxplot$CAI) );
datBoxplot$Expression <- as.factor(datBoxplot$Expression);
datBoxplot$Species <- as.factor(datBoxplot$Species);

datAnalysis$CAI <- as.numeric(as.character(datAnalysis$CAI ));
datAnalysis$Expression <- as.numeric(as.character(datAnalysis$Expression ));
datAnalysis$Species <- as.factor(datAnalysis$Species);
datAnalysis$GeneName <- as.factor(datAnalysis$GeneName);
datAnalysis$CodonLen <- as.numeric(as.character(datAnalysis$CodonLen));

par(mar=c(5.1,12,4.1,2.1));
boxplot( CAI  ~ Expression *Species, data=datBoxplot, notch=TRUE, col=(c( rgb(1/2,0,1,1/3), rgb(0,1/2,1,1/3))), horizontal = T, las=1, ylim=c(0.6,0.73));
nRefLow <- median(datBoxplot[datBoxplot$Species==arrSpecies[1] & datBoxplot$Expression == "lowexp", "CAI"], na.rm = T);
nRefHigh <- median(datBoxplot[datBoxplot$Species==arrSpecies[1] & datBoxplot$Expression == "highexp", "CAI"], na.rm = T);
abline(v=nRefLow , col="red", lwd=2.3);
abline(v=nRefHigh , col="red", lwd=2.3);

dev.off();

write.table(aggregate(CAI ~ Expression + Species, data=datBoxplot, FUN=mean), file = paste(sOut , ".meanCAI.txt", sep=""), sep="\t", quote=F, row.names=F );

#fm1 <- glm( CAI~ log(Expression)/GeneName + Species + CodonLen, data=datAnalysis);

datAnalysis <- within(datAnalysis, Species <- relevel(Species, ref = arrSpecies[1]))
fm2 <- lmer( CAI~ log(Expression) + Species + (1|CodonLen)  , data=datAnalysis, REML=F);
summary(fm2);
sink();
