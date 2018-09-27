#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_above0.2_NVGout/MKTEST_syn_exclCMD");
#compare the DoS between the Dry and Wet populations.
#find genes that have opposing DoS

arrQuantiles <- c(0.25, 0.75);
#arrPops <- c("RACWET.excludehybrid.pergene.alpha.txt", "RACDRY.excludehybrid.pergene.alpha.txt");
arrPops <- c( "ORTWET.pergene.alpha.txt", "ORTDRY.excludehybrid.pergene.alpha.txt");
lsSplit <- strsplit(arrPops, "\\.");
sPop1 <- lsSplit[[1]][[1]]
sPop2 <- lsSplit[[2]][[1]]

sOut <- paste("CompareDoS.", sPop1 , ".", sPop2 ,sep="");
nPopCount <- length(arrPops);

#Get the list of genes that agree in the sign of alpha (>0 or <0) in all populations, and output the lowest p and a multiplied p in all pops:
datCombined <- NULL;
for(sPop in arrPops) {
  datMK <- read.table(sPop, header=T, quote='', sep="\t");
  if (is.null(datCombined) ) {
    datCombined <- datMK[ , c(1,2,3,4,5,6,7,9,10:14)];
    next;
  }
  datCombined <- merge(x=datCombined, y=datMK[, c(1,2,3,4,5,6,7,9,10:14)], by=1:7, all=T);
}

datCombined <- datCombined[complete.cases(datCombined),];
datCombined$deltaDoSP <- apply(datCombined[, c('dN.x','dS.x','pN.x','pS.x','dN.y','dS.y','pN.y','pS.y') ], 1, FUN = function(arr) {chisq.test(matrix( arr, nrow = 2, byrow = T ))$p.value} )
datCombined$deltaDoSQ <- p.adjust(datCombined$deltaDoSP, method = "fdr")
datCombined$dNdSP <- apply(datCombined[, c('dN.x','dS.x','dN.y','dS.y') ], 1, FUN = function(arr) {chisq.test(matrix( arr, nrow = 2, byrow = T ))$p.value} )
datCombined$pNpSP <- apply(datCombined[, c('pN.x','pS.x','pN.y','pS.y') ], 1, FUN = function(arr) {chisq.test(matrix( arr, nrow = 2, byrow = T ))$p.value} )

datCombined$deltaDoS <- datCombined$DoS.x - datCombined$DoS.y ;
fnQuantile <- ecdf(datCombined$deltaDoS);
datCombined$deltaDoSQuantile <- unlist(lapply(datCombined$deltaDoS, fnQuantile));
#datSig <- datCombined[datCombined$deltaDoSP <0.05 & ( datCombined$deltaDoSQuantile < 0.025 | datCombined$deltaDoSQuantile > 0.975 ), ];
datSig <- datCombined[datCombined$deltaDoSP <0.05 & ( datCombined$deltaDoSQuantile < arrQuantiles[1] | datCombined$deltaDoSQuantile > arrQuantiles[2] ) , ];
datSig <- datSig[complete.cases(datSig),]
datSig <- datSig[order(-datSig$deltaDoS), ];
pdf(file=paste(sOut, '.pdf',sep=""), width=5, height=5)
plot(datCombined$deltaDoS, -log10(datCombined$deltaDoSP ), pch=16, col=rgb(0,0,1,1/4) , xlab=paste("DoS(",sPop1,") - DoS(",sPop2,")" ), ylab="-log10(P)");
points(datSig$deltaDoS, -log10(datSig$deltaDoSP ), pch=16, col=rgb(1,0,0,1) )
dev.off();
write.table(datCombined, file=paste(sOut, '.txt',sep=""), quote=F, row.names = F, col.names=T, sep="\t");
write.table(datSig, file=paste(sOut, '.sig.txt',sep=""), quote=F, row.names = F, col.names=T, sep="\t");

# 
# datSig2 <- datCombined[datCombined$dNdSP <0.05 , ];
# plot(datCombined$deltaDoS, -log10(datCombined$dNdSP ), pch=16, col=rgb(0,0,1,1/4));
# points(datSig2$deltaDoS, -log10(datSig2$dNdSP ), pch=16, col=rgb(1,0,0,1) )
# 
# datSig3 <- datCombined[datCombined$pNpSP <0.05 , ];
# plot(datCombined$deltaDoS, -log10(datCombined$pNpSP ), pch=16, col=rgb(0,0,1,1/4));
# points(datSig3$deltaDoS, -log10(datSig3$pNpSP ), pch=16, col=rgb(1,0,0,1) )

