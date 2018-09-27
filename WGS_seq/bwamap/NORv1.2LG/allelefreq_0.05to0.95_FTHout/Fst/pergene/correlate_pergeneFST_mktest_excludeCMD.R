#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_lowfreq_excl/Fst/pergene");

sSuffix <- "_exclCMD";
sPop <- "ORT";
sFileName <- "ORTWET";
#sFileName <- "ORTDRY.excludehybrid";

#sPop <- "RAC";
#sFileName <- "RACWET.excludehybrid";
#sFileName <- "RACDRY.excludehybrid";

sMKFile <- paste("../../../allelefreq_above0.2_NVGout/MKTEST_syn",sSuffix,"/",sFileName, ".pergene.alpha.txt", sep="");


sPopName <- unlist(strsplit(basename(sMKFile) ,'\\.'))[1];

pdf(file=paste(sPopName, sSuffix, ".fst.vs.mk.pdf", sep=""), width=5, height=5);

nMinSites <- 5;

datMK <- read.table( sMKFile, header=T, quote="", sep="\t");
#datMK <- datMK[datMK$DoS < 0, ];

sOut <- paste("zcat -f ", sPop, ".pergene.", ".minsite", nMinSites , ".txt.gz" , sep="");
pOut <- pipe(sOut)

datPerGeneFST <- read.table(pOut , header=T, fill = T);

#datPerGeneFST$fstabsdev <- abs(datPerGeneFST$fst - mean(datPerGeneFST$fst));

datMerge <- merge(datPerGeneFST , datMK[, c(1, 3, 4, 14) ], by.x="OrthoId", by.y = "OrthoID");

oFit <- summary(lm( datMerge$DoS ~ datMerge$fst))

plot(datMerge$fst , datMerge$DoS, cex=0.5, col=rgb(0,0,1,1/10), pch=16, xlab="Fst", ylab=paste(sPopName, "DoS"), main=oFit$coefficients[2,4]);
sColor <- 'grey';
if (oFit$coefficients[2,4] < 0.05) {
  sColor <- 'red';
  
}
abline(a=oFit$coefficients[1,1], b=oFit$coefficients[2,1], col=sColor, lwd=2)
dev.off();
