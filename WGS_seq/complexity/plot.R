#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/finalrun/complexity_all");
arrFiles <- Sys.glob("*_*_extrapolate.txt");
arrSamples <- gsub('_extrapolate\\.txt', '', arrFiles  );
arrSpecies <- gsub('_.+', "", arrSamples);
nSamples <- length(arrSamples);
nPerRow <-1;
nRows <- ceiling(nSamples / nPerRow);
arrCuts <- seq(0.6e7, 5e7, 0.6e7);

#par(mfrow=c(nRows , nPerRow))

fnPlot <- function(sSample) {
  sIn <- paste( sSample, "_extrapolate.txt" , sep="");
  sCountFile <- paste( sSample, "_sampledtotalread.txt" , sep="");
  
  datCounts <- read.table(sCountFile , header = F, skip = 1);
  
  tryCatch( {
    dat <- read.table(sIn, header=T)

    nSampledReads <- sum(datCounts$V1);
    
    plot(dat$TOTAL_READS, dat$UPPER_0.95CI, type="l" , xlab="Total Reads", ylab="Unique Reads", xlim=c(0 , 5e+7) , ylim=c(0 , 5e+7), main=gsub("[0-9A-z]+_[0-9A-z]+_([0-9A-z]+_[0-9A-z]+)" , "\\1", sSample, perl=T));
    polygon(c(dat$TOTAL_READS,rev(dat$TOTAL_READS)),c( dat$UPPER_0.95CI ,rev(dat$LOWER_0.95CI )),col=rgb(1/3,1/3,1,1/4))
    lines(dat$TOTAL_READS, dat$EXPECTED_DISTINCT)
    
    abline(v=nSampledReads , col="red")
    abline(0,1, col="blue");
  } ,
    error = function(e) {
      plot(0,1,main=sSample);
    }
  
  )
}

fnGetUnique <- function(sSample) {
  sIn <- paste( sSample, "_extrapolate.txt" , sep="");
  sCountFile <- paste( sSample, "_sampledtotalread.txt" , sep="");
  
  datCounts <- read.table(sCountFile , header = F, skip = 1);
  
  tryCatch( {
    dat <- read.table(sIn, header=T)
    
    nSampledReads <- sum(datCounts$V1);
    

    return( approx(dat$TOTAL_READS , dat$EXPECTED_DISTINCT , arrCuts )  );

  } ,
  error = function(e) {
    return(NULL);
  }
  
  )
  
}



datRet <- data.frame();
datRet$Sp <- NULL;
datRet$TotalReads <- numeric();
datRet$UniqueRatio <- numeric();

for(i in 1:length(arrSamples) ) {
  lsRet <- fnGetUnique(arrSamples[i] );
  if (is.null(lsRet)) {next;}
  datRet <- rbind( datRet, cbind(Sp = rep(arrSpecies[i] , length(lsRet$x) ) , TotalReads = lsRet$x, UniqueRatio=lsRet$y /lsRet$x ) );
}
datRet$Sp <- as.factor(datRet$Sp);
datRet$TotalReads <- as.numeric( as.character(datRet$TotalReads));
datRet$UniqueRatio <- as.numeric( as.character(datRet$UniqueRatio)) ;

pdf(file = "UniqueRatios.pdf" , width = 7, height=7)
par(mar=c(8.1,4.1,2.1,2.1));
boxplot( UniqueRatio ~  TotalReads * Sp  , dat=datRet, ylim=c(0.8,1), las=2, cex=0.2, col=c("red", "blue"))
abline(h = 0.975 , col="red")
dev.off();
