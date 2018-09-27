
args <- commandArgs(trailingOnly=TRUE)
sPop <- args[1];
#sPop <- 'ORT';

oConn <- pipe(paste("zcat -f ../", sPop, "/fst.txt.gz", sep="") );
datFST <- read.table(oConn , header=F);
oConn2 <- pipe(paste("zcat -f genecoord.txt.gz", sep="") );
datGeneCoord <- read.table(oConn2, header=F, sep="\t");
nMinSites <- as.numeric(args[2]); #10;
#nMinSites <- 5; #10;


sOut <- paste("gzip > ", sPop, ".pergene.", ".minsite", nMinSites , ".txt.gz" , sep="");
pOut <- pipe(sOut, 'w')

#datOut <- NULL;
bWriteHeader <- T;

for(sGroup in levels(datGeneCoord$V1)) {
  #sGroup<-"Group_0_0"
  cat("Doing ", sGroup,"\n");
  datGene <- datGeneCoord[datGeneCoord$V1 == sGroup, ];
  sChr <- as.character(datGene$V5[1]);
  
  nMin <- min( c(datGene$V6 , datGene$V7) );
  nMax <- max( c(datGene$V6 , datGene$V7) );
    datWin <- datFST[datFST$V1 == sChr & datFST$V2 >=nMin & datFST$V2 <= nMax, ];
    if (nrow(datWin) < nMinSites) {
      next;
    }
    nWinFST <- sum(datWin$V4) / sum(datWin$V4 + datWin$V5 + datWin$V6); 
    #datOut <- rbind(datOut, data.frame('chr' = sChr,'left' = nStart, 'right'=nStart+nWindowSize,  'midpoint' = (nStart + 0.5*nWindowSize), 'fst' = nWinFST ));
    datOut <- data.frame('OrthoId'=sGroup ,  'chr' = sChr,'left' = nMin, 'right'=nMax,  'midpoint' = mean(nMin, nMax), 'fst' = nWinFST );
    write.table(datOut, file=pOut , col.names = bWriteHeader, append = (!bWriteHeader), row.names = F, quote = F, sep="\t");
    bWriteHeader <- F;
  
}

close(pOut)

#write.table(datOut, file=sOut , col.names = T, row.names = F, quote = F, sep="\t");
