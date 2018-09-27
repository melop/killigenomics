#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_0.01to0.98_NVGout/anavar/runanavar");
#sSrcDir <- "RACDRY.excludehybrid";
#sSrcDir <- "RACWET.excludehybrid";
sSrcDir <- "ORTDRY.excludehybrid";
#sSrcDir <- "ORTWET";
arrRunTypes <- c('full', 'sametheta_', 'sametheta_nopolerr_');
arrParamCounts <- c(59, 58, 54);
nBestTolerance <- 0.05; #tolerate this much difference in likelihood
lsDat <- list();
lsBest <- list();

options(digits=16)

for( i in 1:length(arrRunTypes) ) {
  sRunType <- arrRunTypes[i];
  if (sRunType == 'full') {
    sRunType <- '';
  }
  nParamCount <- arrParamCounts[i];
  arrFiles <- Sys.glob(paste(sSrcDir,'/',sRunType, 'out.rep*.txt', sep=""));
  datAnavar <- NULL;
  for(sFile in arrFiles) {
    datAnavar <- rbind(datAnavar, read.table(sFile, header=T, skip = 6, sep="\t"));
    
  }
  
  #adjust the order of gamma categories, such that sel_gamma_1 < sel_gamma_2 < sel_gamma3
  
  arrColC1 <- grep('sel_.*_1' , colnames(datAnavar));
  arrColC2 <- grep('sel_.*_2' , colnames(datAnavar));
  arrColC3 <- grep('sel_.*_3' , colnames(datAnavar));
  datCols <- rbind(arrColC1, arrColC2, arrColC3);
  
  for(nRow in 1:nrow(datAnavar)) {
    datVals <- data.frame(); 
    for(nR2 in 1:nrow(datCols) ) {
      datVals <- rbind(datVals, as.numeric(datAnavar[nRow ,datCols[nR2,]]) );
    }
    colnames(datVals) <- c('sel_theta' , 'sel_gamma' , 'sel_e' );
    datVals <- datVals[order(datVals$sel_gamma) , ];
    
    for(nR2 in 1:nrow(datVals) ) {
      datAnavar[nRow, datCols[nR2,] ] <- datVals[nR2, ];
    }
    
  }
  
  datAnavar <- datAnavar[order(-datAnavar$lnL) , ];
  datBest <- datAnavar[1,];
  lsBest[[i]] <- datBest;
  lsDat[[i]] <- datAnavar;
}

sOutStem <- paste(basename(sSrcDir), ".out", sep="");
sink(file=paste(sOutStem, ".txt", sep=""));
pdf(file=paste(sOutStem, ".pdf", sep=""), width=7, height=5);

#perform tests LRT
for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  nLnL1 <- lsBest[[i]]$lnL;
  nAIC <- 2*nParamCount1 - 2*nLnL1;
  cat(sRunType1, "AIC=",nAIC, "\n");
}

for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  nLnL1 <- lsBest[[i]]$lnL;

  for( j in 1:length(arrRunTypes) ) {
    if (i == j) {
      next;
    }
    
    sRunType2 <- arrRunTypes[j];
    nParamCount2 <- arrParamCounts[j];
    nLnL2 <- lsBest[[j]]$lnL;
    
    if (nParamCount2 > nParamCount1) {
      next;
    }
    
    #model 1 is the more complex model.
    nStat <- 2 * (nLnL1 - nLnL2);
    nDF <- nParamCount1 - nParamCount2;
    nP <- 1 - pchisq(nStat,nDF);
    cat(sRunType1, "vs.", sRunType2, "chistat=", nStat , "LRT=", nP, "d.f.=",nDF, "\n");
  }
  
}


for( i in 1:length(arrRunTypes) ) {
  sRunType1 <- arrRunTypes[i];
  nParamCount1 <- arrParamCounts[i];
  datAnavar <- lsDat[[i]][ (lsBest[[i]]$lnL - lsDat[[i]]$lnL) < nBestTolerance , ];
  datAnavar[, c('sel_theta_prop_1','sel_theta_prop_2', 'sel_theta_prop_3')] <- datAnavar[, grep('sel_theta_', colnames(datAnavar)) ] / rowSums(datAnavar[, grep('sel_theta_', colnames(datAnavar)) ]);
  nTotal <- 10000;
  arrG1 <- rep(datAnavar$sel_gamma_1[1], datAnavar[1,'sel_theta_prop_1'] * nTotal );
  arrG2 <- rep(datAnavar$sel_gamma_2[1], datAnavar[1,'sel_theta_prop_2'] * nTotal );
  arrG3 <- rep(datAnavar$sel_gamma_3[1], datAnavar[1,'sel_theta_prop_3'] * nTotal );
  
  hist(arrG1, col='blue', lty="blank", xlim = c(min(datAnavar$sel_gamma_1)-10, max(datAnavar$sel_gamma_3) )+10, main=sRunType1 , freq = T);
  hist(arrG2, col='darkgreen' , lty="blank", add=T, freq = T);
  hist(arrG3, col='red' , lty="blank", add=T, freq = T);
  text(arrG1[1]+(-arrG1[1])/4, 3500, paste(round(datAnavar[1,'sel_theta_prop_1'] * 100) ,"% : ", "Gamma 1 = " , round(datAnavar$sel_gamma_1[1], digits=2),"\n",
                        round(datAnavar[1,'sel_theta_prop_2'] * 100) ,"% : ", "Gamma 2 = " ,  round(datAnavar$sel_gamma_2[1], digits=2),"\n",
                        round(datAnavar[1,'sel_theta_prop_3'] * 100) ,"% : ", "Gamma 3 = " ,  round(datAnavar$sel_gamma_3[1], digits=2),"\n",
                        sep=""), cex=1.4, adj=c(0,0.5));
 # boxplot(datAnavar$sel_gamma_1, datAnavar$sel_gamma_2, datAnavar$sel_gamma_3);
  
}

sink();
dev.off();
