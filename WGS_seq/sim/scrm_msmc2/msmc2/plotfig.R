setwd(".");

pdf("2pops_sim.pdf", width=7, height=6);
sMainFolder <- "msmc2ret/";

sRealFolder <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/msmc2/"; #the folder with form_RAC.txt and form_ORT.txt
arrYLim <- c(5e4, 1e6);
arrXLim <- c(1e4, 1.5e6);
nMu <- 2.6321e-9;


fnPlotPopPair <- function(sSp, sPop1, sPop2, arrCol, bBSPlot = F, sThisBSFile='', sRealFile, nSplitTime=0) {
 
  datMSMC1 <- read.table(paste(sMainFolder,"/",sPop1, ".final.txt", sep=""), header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  
  datMSMC2 <- read.table(paste(sMainFolder,"/",sPop2, ".final.txt", sep=""), header=T, sep="\t" );
  datMSMC2$gen <- datMSMC2$left_time_boundary/nMu;
  datMSMC2$gen[1] <- 0.01;
  datMSMC2$popsize <- (1/datMSMC2$lambda)/(2*nMu);
  
  
  if (bBSPlot == F) {
    plot( datMSMC1$gen ,datMSMC1$popsize , main=sSp,  type = 's' , log='xy', xlim=arrXLim, ylim = arrYLim, xlab="Generations" , ylab = "Pop size", lwd=3, col=arrCol[2])
    abline(v= nSplitTime , col=arrCol[1], lwd=2);
    
    # legend("topright", legend = c(sPop1, sPop2),
    #        text.width = strwidth("1,000,000"),
    #        lty = 1, lwd=3, col = arrCol[2:3],  xjust = 1, yjust = 1,
    # )
    
  } else {
    lines( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , lwd=3, col=arrCol[2])
    
  }
  
  lines( datMSMC2$gen ,datMSMC2$popsize ,  type = 's' , lwd=3, col=arrCol[3])
  
  
    #fnPlotVariation(sPop1, arrCol[2]);
  #fnPlotVariation(sPop2, arrCol[3]);
  if (bBSPlot==TRUE) {
    return();
  } else {
    datReal <- read.table( paste(sRealFolder, '/',sRealFile ,sep=""), header=T, sep="\t");
    points(datReal$Pop1_Gen, datReal$Pop1_Popsize, col=arrCol[2])
    points(datReal$Pop2_Gen, datReal$Pop2_Popsize, col=arrCol[3])
    
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
fnPlotPopPair('ORT' , 'ORTWET', 'ORTDRY', c('purple', 'blue', 'red') , sRealFile="ORT_pop_sizes.txt", nSplitTime=205.35e3);
dev.off();

