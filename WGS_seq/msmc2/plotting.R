setwd("../WGS_seq/bwamap/NORv1.2LG/msmc2/");
pdf("2pops.pdf", width=7, height=6);
sMainFolder <- "msmc2ret2/";
sBSFolder <- "bootstrapped/"
nBSReps <- 30;

arrYLim <- c(4e4, 2e6);
arrXLim <- c(1e4, 3e6);
nMu <- 2.6321e-9;


fnPlotPopPair <- function(sSp, sPop1, sPop2, arrCol, bBSPlot = F, nSplitTime=0, sThisBSFile1='',sThisBSFile2='', sRealFile="") {
 
  sInFile1 <- paste(sMainFolder,"/",sPop1, ".final.txt", sep="");
  sInFile2 <- paste(sMainFolder,"/",sPop2, ".final.txt", sep="");
  
  if (bBSPlot) {
    cat("bs: ", sThisBSFile1, sThisBSFile2,"\n");
    sInFile1 <- sThisBSFile1;
    sInFile2 <- sThisBSFile2;
    arrCol <- add.alpha(arrCol, 0.2);
  }
  
  cat("Open ", sInFile1,"\n");
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  
  cat("Open ", sInFile2,"\n");
  datMSMC2 <- read.table(sInFile2, header=T, sep="\t" );
  datMSMC2$gen <- datMSMC2$left_time_boundary/nMu;
  datMSMC2$gen[1] <- 0.01;
  datMSMC2$popsize <- (1/datMSMC2$lambda)/(2*nMu);
  
  
  if (bBSPlot == F) {
    plot( datMSMC1$gen ,datMSMC1$popsize , main=sSp,  type = 's' , log='xy', xlim=arrXLim, ylim = arrYLim, xlab="Generations" , ylab = "Pop size", lwd=3, col=arrCol[2])
    # legend("topright", legend = c(sPop1, sPop2),
    #        text.width = strwidth("1,000,000"),
    #        lty = 1, lwd=3, col = arrCol[2:3],  xjust = 1, yjust = 1,
    # )
    # 
    abline(v= nSplitTime , col=arrCol[1], lwd=2);
    
  } else {
    lines( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , lwd=3, col=arrCol[2])
    
  }
  
  lines( datMSMC2$gen ,datMSMC2$popsize ,  type = 's' , lwd=3, col=arrCol[3])
  
  

  if (bBSPlot==TRUE) {
    #plot variation
    return();
  } else {
    cat("Plot bs...\n");
    for (nRep in 1:nBSReps) {
      sBSF1 <- paste(sBSFolder,'/',sPop1, '/_', nRep, '/out.final.txt' , sep="");
      sBSF2 <- paste(sBSFolder,'/',sPop2, '/_', nRep, '/out.final.txt' , sep="");
      #cat(sBSF1);
      if (file.exists(sBSF1) && file.exists(sBSF2)) {
        fnPlotPopPair(sSp, sPop1, sPop2, arrCol, bBSPlot = T, sThisBSFile1=sBSF1,sThisBSFile2=sBSF2)
      }
    }
    
  }
  
  #construct output table:
  return();
  
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

fnPlotPopPair('RAC' , 'RACWET', 'RACDRY', c('purple', 'blue', 'red') , sRealFile="RAC_pop_sizes.txt", nSplitTime=258.2e3);
fnPlotPopPair('ORT' , 'ORTWET', 'ORTDRY', c('purple', 'blue', 'red') , sRealFile="ORT_pop_sizes.txt" , nSplitTime=205.35e3);
dev.off();

