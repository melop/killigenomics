sRootFolder <- "../pseudogenomes/cov_blacklist/";

arrRefFolders <- c('AAU' , 'NOR' , 'PLP' , 'CTO');
for(sRefFolder in arrRefFolders) {
  arrSubFolders <- Sys.glob( paste(sRootFolder,"/",sRefFolder, "/*", sep="") );
  for(sSpFolder in arrSubFolders) {
    dat <- read.table(paste(sSpFolder, "/out.txt", sep='') , header=F, colClasses = c('character', 'character', 'numeric', 'numeric', 'character','numeric', 'numeric', 'numeric' ) );
    
    arrQuantiles <- quantile(dat$V8, c(0.5, 0.75));
    nDiff <- 3 * (arrQuantiles[2] - arrQuantiles[1]);
    arrXRange <- c(arrQuantiles[1]-nDiff, arrQuantiles[1]+nDiff);
    
    sCondition <- ( dat$V8 <= max(arrXRange) & dat$V8 >= min(arrXRange) );
    arrVals <- dat[sCondition , 'V8' ];
    names(arrVals) <- dat[sCondition , 'V1' ];
    
    nMean <- mean(arrVals) ;
    nSd <- sd(arrVals);
    pdf(file=paste(sSpFolder, "/cov.pdf", sep='') );
    hist(dat$V8, breaks=3000, xlim=c(0, max(arrXRange)+20), ylim=c(0,0.3), freq = F, xlab = paste("Coverage (X)\n", "mean: ", nMean , " sd:",  nSd) )
    lines( seq(0,100,0.1) , dnorm( seq(0,100,0.1) , mean = nMean, sd=nSd) , col='red');
    dev.off();
    dat$P <- pnorm(dat$V8, mean = nMean, sd=nSd );
    dat[dat$V8 <= nMean, 'P'] <- 1; #if it's lower than mean always set to 1
    dat[dat$V8 > nMean, 'P'] <-  2*( 1-dat[dat$V8 > nMean, 'P']); 
    dat$Excluded <- "Keep";
    dat[dat$P < 0.001, 'Excluded'] <- "Exclude";
    write.table(dat, file=paste(sSpFolder, "/out_Pvalues.txt", sep=''), col.names = F, row.names = F, quote = F, sep="\t");
  }

}

#now read in all P values
for(sRefFolder in arrRefFolders) {
  arrSubFolders <- Sys.glob( paste(sRootFolder,"/",sRefFolder, "/*", sep="") );
  datAllSp <- NULL;
  for(sSpFolder in arrSubFolders) {
      dat2 <- read.table(paste(sSpFolder, "/out_Pvalues.txt", sep='') , header=F, colClasses = c('character', 'character', 'numeric', 'numeric', 'character','numeric', 'numeric', 'numeric', 'numeric', 'character' ) );
      if (is.null(datAllSp)) {
        datAllSp <- dat2[, c(1,2,3,4,5,7)];
      }
      
      sSp <- basename(sSpFolder);
      #CHECK TO SEE IF ORDER IS CONSISTENT
      if (!identical(dat2$V1, datAllSp$V1)) {
        cat("ERROR, gene name order differ!\n");
        break;
      }
      datAllSp[ , sSp] <- dat2$V9;
    
  }
  
  datAllSp$Exclude <- F;
  for(nRow in 1:nrow(datAllSp)) {
    if (any( datAllSp[nRow, 7: (ncol(datAllSp)-1)] < 0.001) ) {
      datAllSp[nRow, 'Exclude'] <- T;
    }
  }
  #View(datAllSp[datAllSp$Exclude,]);
  write.table(datAllSp, file=paste(getwd(), "/", sRefFolder, "_blacklist.txt", sep=''), col.names = F, row.names = F, quote = F, sep="\t");
}








