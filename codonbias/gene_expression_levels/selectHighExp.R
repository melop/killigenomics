arrSpp <- c("NOR", 'CTO', 'PLP', 'AAU');

datAll <- NULL;
for(sSp in arrSpp) {
  sIn <- paste("counts_.", sSp ,".txt", sep='');
  dat <- read.table(sIn, header=T, sep="\t", stringsAsFactors = F, quote='');
  dat$RPKM <- dat[, 11] / (dat[, 10]/1000) / ( sum(dat[, 11]) / 1e6 );
  if (is.null(datAll)) {
    datAll <- dat[, c(1,2,3,12)];
    colnames(datAll) <- c(colnames(dat)[1:3], sSp);
  } else {
    colnames(dat)[12] <- sSp;
    datAll <- merge(datAll, dat[, c(1,12)], by='OrthoID');
  }
}

datAll$AvgRPKM <- apply(datAll[, seq(4, ncol(datAll))] , 1, mean);
datAll$SdRPKM <- apply(datAll[, seq(4, ncol(datAll))] , 1, sd);
datAll$SdOverAvg <- datAll$SdRPKM / datAll$AvgRPKM;

nSdCutoff <- quantile(datAll$SdOverAvg , probs = 0.05, na.rm=T);
nRPKMCutoff <- quantile(datAll$AvgRPKM , probs = 0.95, na.rm=T);

datHighExpGenes <- datAll[datAll$AvgRPKM >= nRPKMCutoff & datAll$SdOverAvg <= nSdCutoff, ];
write.table(datAll, file="all.expression.level.txt" , col.names = T, row.names = F, quote = F, sep="\t");
write.table(datHighExpGenes, file="high.expression.level.txt" , col.names = T, row.names = F, quote = F, sep="\t");
