setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/Fst");

args <- commandArgs(trailingOnly=TRUE)
sPop <- args[1];

oConn <- pipe(paste("zcat -f ", sPop, "/fst.txt.gz", sep="") );
datFST <- read.table(oConn , header=F);
nWindowSize <- as.numeric(args[2]); #50000;
nStepSize <- as.numeric(args[3]); #10000;
nMinSites <- as.numeric(args[4]); #10;

sOut <- paste("gzip > ", sPop, ".slideFST.win",nWindowSize, ".step", nStepSize, ".minsite", nMinSites , ".txt.gz" , sep="");
pOut <- pipe(sOut, 'w')

#datOut <- NULL;
bWriteHeader <- T;

for(sChr in unique(datFST$V1)) {
  cat("Doing ", sChr,"\n");
  datChr <- datFST[datFST$V1 == sChr, ];
  nMin <- 1;#min(datFST$V2);
  nMax <- max(datFST$V2);
  for(nStart in seq(nMin, nMax-nWindowSize, nStepSize) ) {
    datWin <- datChr[datChr$V2 >=nStart & datChr$V2 < (nStart+nWindowSize), ];
    if (nrow(datWin) < nMinSites) {
      next;
    }
    nWinFST <- sum(datWin$V4) / sum(datWin$V4 + datWin$V5 + datWin$V6); 
    #datOut <- rbind(datOut, data.frame('chr' = sChr,'left' = nStart, 'right'=nStart+nWindowSize,  'midpoint' = (nStart + 0.5*nWindowSize), 'fst' = nWinFST ));
    datOut <- data.frame('chr' = sChr,'left' = nStart, 'right'=nStart+nWindowSize,  'midpoint' = (nStart + 0.5*nWindowSize), 'fst' = nWinFST );
    write.table(datOut, file=pOut , col.names = bWriteHeader, append = (!bWriteHeader), row.names = F, quote = F, sep="\t");
    bWriteHeader <- F;
  }
}

#write.table(datOut, file=sOut , col.names = T, row.names = F, quote = F, sep="\t");