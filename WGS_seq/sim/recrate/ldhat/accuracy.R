setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/sim/recrate/ldhat");

sRealRate <- "../recmap.txt";
arrPops <- c('ORTWET', 'ORTDRY');

arrCols <- c('blue', 'purple');

bPlottedTrue <- F;

nSmoothWindow <- 300;

fnSmooth <- function(datS, nWindowSize) {
  nMaxPos <- max(datS$Loci);
  datRet <- data.frame();
  for(nStart in seq(1,nMaxPos-nWindowSize, nWindowSize) ) {
    nEnd <- nStart + nWindowSize;
    datWin <- datS[datS$Loci >= nStart & datS$Loci <= nEnd, ];
    if (nrow(datWin) >= 1) {
      datRet <- rbind(datRet , list( Loci = (nEnd-nWindowSize/2), Mean_rho=mean(datWin$Mean_rho),Median= median(datWin$Median) ) )
    }
  }
  return(datRet);
}

datRealRate <- read.table(sRealRate, header=F);
pdf(file=paste("ldhat_accuracy_2.3xcov.smooth",nSmoothWindow,"kb.pdf", sep="" ), width=5, height=5);

for(nPop in 1:length(arrPops) ) {

  sPop <- arrPops[nPop];
  sCol <- arrCols[nPop];
  
  sIn <- paste('ldhat_pen50', '/', sPop, '/chr1/res.txt', sep='');
  
  datRate <- read.table(sIn, header = T);
  
  datRate <- datRate[datRate$Loci > 0, ];
  datS <- fnSmooth(datRate , nSmoothWindow);
  datS$normalize_rho <- datS$Mean_rho / mean(datS$Mean_rho);
  
  if (!bPlottedTrue) {
    
    datRealRate$Loci <- (datRealRate$V1 + datRealRate$V2)* max(datRate$Loci) / 2;
    datRealRate$normalize_rho <- datRealRate$V3 /mean(datRealRate$V3);
    
    plot(datRealRate$Loci, log10(datRealRate$normalize_rho),  col='red', pch=16, xlab = "Position (Kb)", ylab="Log10(normalized rho)");
    bPlottedTrue <- T;
  }
  lines(datS$Loci, log10(datS$normalize_rho), type = "l" , lwd=2, col=sCol);
}

dev.off();