setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/PLPv1.2LG/msmc2");
#pdf("allpops.pdf", width=7, height=6);
#pdf("praslin.pdf", width=7, height=6);
#pdf("mahesouth.pdf", width=7, height=6);
pdf("allwild.pdf", width=7, height=6);

sMainFolder <- "msmc2ret/";
sBSFolder <- "bootstrapped/"
nBSReps <- 30;

sRealFolder <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/smcpp_raw_hicov_realigned_varqualmask/"; #the folder with form_RAC.txt and form_ORT.txt
arrYLim <- c(1e3, 5e5);
arrXLim <- c(5e3, 0.5e6);
#nMu <- 2.6321e-9;
nMu <- 3.72221404879E-9; # estimated from divergence between PML and PLP, this is per year

fnPlotPop <- function(sPop, arrCol, bBSPlot = F, bAppendPlot=T, sThisBSFile1 = "") {
 
  sInFile1 <- paste(sMainFolder,"/",sPop, ".final.txt", sep="");
  arrCol <- add.alpha(arrCol, 0.8);
  
  if (bBSPlot) {
    cat("bs: ", sThisBSFile1,"\n");
    sInFile1 <- sThisBSFile1;
    arrCol <- add.alpha(arrCol, 0.2);
  }
  
  cat("Open ", sInFile1,"\n");
  datMSMC1 <- read.table(sInFile1, header=T, sep="\t" );
  datMSMC1$gen <- datMSMC1$left_time_boundary/nMu;
  datMSMC1$gen[1] <- 0.01;
  datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*nMu);
  
  
  if (bBSPlot == F && (!bAppendPlot) ) {
    plot( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , log='xy', xlim=arrXLim, ylim = arrYLim, xlab="Generations" , ylab = "Pop size", lwd=3, col=arrCol[1])
  } else {
    lines( datMSMC1$gen ,datMSMC1$popsize ,  type = 's' , lwd=3, col=arrCol[1])
    
  }
  

  if (bBSPlot==TRUE) {
    #plot variation
    return();
  } else {
    cat("Plot bs...\n");
    for (nRep in 1:nBSReps) {
      sBSF1 <- paste(sBSFolder,'/',sPop, '/_', nRep, '/out.final.txt' , sep="");
      #cat(sBSF1);
      if (file.exists(sBSF1)) {
        fnPlotPopPair(sPop, arrCol, bBSPlot = T, sThisBSFile1=sBSF1 )
      }
    }
    
    # datReal <- read.table( paste(sRealFolder, '/',sRealFile ,sep=""), header=T, sep="\t");
    # points(datReal$Pop1_Gen, datReal$Pop1_Popsize, col=arrCol[2])
    # points(datReal$Pop2_Gen, datReal$Pop2_Popsize, col=arrCol[3])
    # 
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

#Wild Praslin:
fnPlotPop('PLP_A_PLP-A3_D501_D701', 'blue' , bAppendPlot = F);
fnPlotPop('PLP_A_PLP-A4_D502_D702', 'blue' , bAppendPlot = T);
fnPlotPop('PLP_B_PLP-B3_D504_D704', 'red' , bAppendPlot = T);
fnPlotPop('PLP_C_PLP-C1_D506_D706', 'orange' , bAppendPlot = T);
fnPlotPop('PLP_C_PLP-C2_D507_D707', 'orange' , bAppendPlot = T);
fnPlotPop('PLP_C_PLP-C3_D508_D708', 'orange' , bAppendPlot = T);
fnPlotPop('PLP_C_PLP-C4_D502_D701', 'orange' , bAppendPlot = T);
fnPlotPop('PLP_F_PLP-F1_D503_D701', 'purple' , bAppendPlot = T);
fnPlotPop('PLP_F_PLP-F2_D504_D702', 'purple' , bAppendPlot = T);
fnPlotPop('PLP_F_PLP-F3_D505_D703', 'purple' , bAppendPlot = T);

#Domesticated Praslin (La Digue)
#fnPlotPop('PLP_D_PLP-D0_ref', 'darkgreen' , bAppendPlot = T);
#fnPlotPop('PLP_D_PLP-D1_D503_D702', 'darkgreen' , bAppendPlot = T);

#Domesticated Mahe' North.
# fnPlotPop('PLP_G_PLP-G1_D506_D704', 'yellow' , bAppendPlot = T);
# fnPlotPop('PLP_G_PLP-G3_D508_D706', 'yellow' , bAppendPlot = T);
# fnPlotPop('PLP_G_PLP-G5_D502_D708', 'yellow' , bAppendPlot = T);
# Wild Mahe South
fnPlotPop('PLP_E_PLP-E1_D505_D704', 'turquoise' , bAppendPlot = T);
fnPlotPop('PLP_E_PLP-E2_D506_D705', 'turquoise' , bAppendPlot = T);
fnPlotPop('PLP_E_PLP-E3_D507_D706', 'turquoise' , bAppendPlot = T);
fnPlotPop('PLP_E_PLP-E4_D508_D707', 'turquoise' , bAppendPlot = T);
fnPlotPop('PLP_E_PLP-E5_D501_D708', 'turquoise' , bAppendPlot = T);

# fnPlotPop('PLP_E_PLP-E5_D501_D708', 'black' , bAppendPlot = T);


dev.off();

