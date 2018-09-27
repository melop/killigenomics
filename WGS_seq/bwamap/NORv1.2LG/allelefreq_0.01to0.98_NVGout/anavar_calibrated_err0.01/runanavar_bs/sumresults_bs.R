#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_0.01to0.98_NVGout/anavar_calibrated_err0.01/runanavar_bs");
#sSrcDir <- "RACDRY.excludehybrid";
#sSrcDir <- "RACWET.excludehybrid";
#sSrcDir <- "ORTDRY";
#sSrcDir <- "ORTWET";

arrPops <- c('ORTDRY', 'ORTWET');
#arrPops <- c('RACDRY', 'RACWET');

arrXLim <- c(-20, 10);
arrYLim <- c(0, 0.30);

# arrXLim <- c(-1000, 0);
# arrYLim <- c(0.5, 1);


sOutStem <- paste( 'out', paste(arrPops, collapse = '-') , 'x',paste(arrXLim, collapse = '-'), 'y',  paste(arrYLim, collapse = '-'), sep="_");
sOutStem2 <- paste( 'out', paste(arrPops, collapse = '-') , sep="_");

lsCol <- list( c(rgb(1,3/16,0,1/3), rgb(1,9/16,3/16,1/2), rgb(12/16,0,3/12,1/3) ), c(rgb(0,3/16,1,1/3), rgb(6/16,0,1,1/3), rgb(3/16,12/16,0,1/3) ) );
arrRunTypes <- c('full');
sRunTypeOrig <- 'full'
options(digits=16)
datOutStats <- data.frame();
  
bPlotMade <- F;

datAnavarPrev <- NULL;
datOriginalPrev <- NULL;
pdf(file=paste(sOutStem,'.pdf', sep=''), width=5, height=5 )
for( i in 1:length(arrPops) ) {

  sSrcDir <- arrPops[[i]];
  arrCols <- lsCol[[i]];
  arrOrigFile <- Sys.glob(paste("../runanavar/", arrPops[[i]],"*.out_",sRunTypeOrig ,".txt", sep=""));
  datOriginal <- read.table(arrOrigFile[1], header=T, sep="\t");
  
  datAllBsReps <- data.frame();
  
  if (sRunTypeOrig == 'full') {
    sRunType <- '';
  }
  
  arrFiles <- Sys.glob(paste(sSrcDir,'/',sRunType, 'out.bsrep*.txt', sep=""));
  datAnavar <- NULL;
  for(sFile in arrFiles) {
    datAnavar <-  read.table(sFile, header=T, skip = 6, sep="\t");
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
    datAllBsReps <- rbind(datAllBsReps, datBest);
  }
  
  
    sRunType1 <- sRunType;
    datAnavar <- datAllBsReps;
    datAnavar[, c('sel_theta_prop_1','sel_theta_prop_2', 'sel_theta_prop_3')] <- datAnavar[, grep('sel_theta_', colnames(datAnavar)) ] / rowSums(datAnavar[, grep('sel_theta_', colnames(datAnavar)) ]);
    
    if (!bPlotMade) {
      plot(datAnavar[,'sel_gamma_1'], datAnavar[,'sel_theta_prop_1'], col=arrCols[1], xlim=arrXLim, ylim=arrYLim, pch=16, cex=1, xlab="Selection Coefficient (Gamma)", ylab="Proportion of Mutations");
    } else {
      points(datAnavar[,'sel_gamma_1'], datAnavar[,'sel_theta_prop_1'], col=arrCols[1], pch=16, cex=1);
      
    }
    points(datAnavar[,'sel_gamma_2'], datAnavar[,'sel_theta_prop_2'], col=arrCols[2] , pch=16, cex=1);
    points(datAnavar[,'sel_gamma_3'], datAnavar[,'sel_theta_prop_3'], col=arrCols[3], pch=16, cex=1);
    points(datOriginal[1, 'sel_gamma_1'], datOriginal[1, 'sel_theta_prop_1'], col='black', pch=16,cex=2);
    points(datOriginal[1, 'sel_gamma_2'], datOriginal[1, 'sel_theta_prop_2'], col='black', pch=16,cex=2);
    points(datOriginal[1, 'sel_gamma_3'], datOriginal[1, 'sel_theta_prop_3'], col='black', pch=16,cex=2);
    
    datThisStat <- data.frame();
    arrColumns <-c('sel_theta_prop_1' , 'sel_gamma_1','sel_theta_prop_2' , 'sel_gamma_2','sel_theta_prop_3' , 'sel_gamma_3' );
    for(sColumn in arrColumns) {
      datThisStat <- rbind(datThisStat, c(datOriginal[1, sColumn], sd(datAnavar[,sColumn]), quantile(datAnavar[,sColumn], c(0.025,0.975) )));
    }
    rownames(datThisStat) <- arrColumns;
    colnames(datThisStat) <- paste(sSrcDir, c('mean', 'sd', '0.025', '0.975'), sep='_' );
    
    if (!bPlotMade) {
      datOutStats <-  datThisStat;
      
    } else {
      #perform test
      datOutStats <- cbind(datOutStats, datThisStat);
      arrP12s <- NULL;
      arrP21s <- NULL;
      arrPEmp <- NULL;
      
      for(sParam in arrColumns) {
        # oTest <- wilcox.test(datAnavarPrev[, sParam], datAnavar[, sParam], paired = F);
        # arrPs <- c(arrPs, oTest$p.value);
        nPrevMean <- datOriginalPrev[1, sParam];
        nCurrMean <- datOriginal[1, sParam];
        arrPrevBS <- datAnavarPrev[, sParam];
        arrCurrBS <- datAnavar[, sParam];
        
        if (nCurrMean > nPrevMean) {
          nProb12 <- length(arrPrevBS[arrPrevBS>nCurrMean])/length(arrPrevBS);
          nProb21 <- length(arrCurrBS[arrCurrBS<nPrevMean])/length(arrCurrBS);
          nProbOverlap <- length( arrPrevBS[arrPrevBS>=min(arrCurrBS)] ) / length(arrCurrBS);
          
        } else {
          nProb12 <- length(arrPrevBS[arrPrevBS<nCurrMean])/length(arrPrevBS);
          nProb21 <- length(arrCurrBS[arrCurrBS>nPrevMean])/length(arrCurrBS);
          nProbOverlap <- length( arrPrevBS[arrPrevBS<=max(arrCurrBS)] ) / length(arrCurrBS);
          
        }
        
        arrP12s <- c(arrP12s , nProb12);
        arrP21s <- c(arrP21s , nProb21);
        arrPEmp <- c(arrPEmp, nProbOverlap)
        
      }
      # datOutStats <- cbind(datOutStats, data.frame(wilcoxP = arrPs ) );
      datOutStats <- cbind(datOutStats, data.frame( P12 = arrP12s, P21 = arrP21s, POverlap = arrPEmp ) );
      
    }
    bPlotMade <- T;
    datAnavarPrev <- datAnavar;
    datOriginalPrev <- datOriginal;
}
dev.off();
write.table(datOutStats, file=paste(sOutStem2, ".stats.txt",sep=''), col.names = T, row.names = T, sep="\t", quote = F);
