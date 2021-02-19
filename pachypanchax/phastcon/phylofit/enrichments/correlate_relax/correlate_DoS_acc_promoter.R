setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/correlate_relax");

arrDistance <- c(3, seq(10, 100, 10) );
#arrDistance <- seq(0, 50, 10) ;

sK <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/PLPv1.2LG/allelefreq_0.05to1/MKTEST_syn/PLP.PDC.pergene.alpha.txt"

datK <- read.table(pipe(paste("{ head -n1 ", sK," ; grep -v GeneName ", sK,"; }", sep='')) , header=T, sep="\t" , quote="", fill = T, stringsAsFactors = F);
#datK <- datK[datK$Success == 'Success',];
#arrRelax <- datK[datK$DoS < 0 , 'OrthoID'];
#arrNotRelax <- datK[datK$DoS >= 0 , 'OrthoID'];


fnPlot <- function(bBootstrap) {
  datRet <- NULL;
  for(nDist in arrDistance) {
    sTarget <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/GO/",nDist,"kb/target.txt" , sep="");
    sBg <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/GO/",nDist, "kb/bg.txt" , sep="");
    #sTarget <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/GO/",nDist,"-",nDist+10, "kb/target.txt" , sep="");
    #sBg <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/GO/",nDist,"-",nDist+10, "kb/bg.txt" , sep="");
    
    datTarget <- read.table(sTarget, header=F, sep="\t");
    datBg <- read.table(sBg, header=F, sep="\t");
    
    
    if (bBootstrap) {
      datTarget <- datTarget[ sample(1:nrow(datTarget), size = nrow(datTarget), replace = T), ];
      datBg <- datBg[ sample(1:nrow(datBg), size = nrow(datBg), replace = T), ];
      
    }
    
    datBg <- datBg[!(datBg$V1 %in% datTarget$V1), ];
    
    datBg <- merge(datBg, datK[ , c('OrthoID', 'DoS') ], by.x = 'V1', by.y = 'OrthoID');
    datTarget <- merge(datTarget, datK[ , c('OrthoID', 'DoS') ], by.x = 'V1', by.y = 'OrthoID');

    nAvgBGDoS <- median(datBg$DoS);
    nAvgTargetDoS<- median(datTarget$DoS);
    
    oTest <- wilcox.test(datBg$DoS, datTarget$DoS)
    datRow <- data.frame(AccElementFromGeneCutoff=nDist, fisherP=oTest$p.value, BgRelaxRatio=nAvgBGDoS, TargetRelaxRatio=nAvgTargetDoS );
    datRet <- rbind(datRet, datRow);
    
    # arrKTarget <- datK$K_q_free[datK$GeneName %in% datTarget$V1];
    # arrKBg <- datK$K_q_free[datK$GeneName %in% datBg$V1];
    # median(arrKBg);
    # median(arrKTarget);
    # wilcox.test(arrKTarget, arrKBg)
    
  }
  
  if (!bBootstrap) {
    plot(datRet$AccElementFromGeneCutoff , datRet$TargetRelaxRatio, col='red', type = 'l', lwd=3, ylim=c(-0.6, -0.1), xlab="Distance from start codon", ylab="Median DoS"  );
    lines(datRet$AccElementFromGeneCutoff , datRet$BgRelaxRatio, col='blue', lwd=3)
    #View(datRet);
    datRetSig <- datRet[datRet$fisherP<0.05,]
    points(datRetSig$AccElementFromGeneCutoff , datRetSig$TargetRelaxRatio - 0.01, pch=8, col='red');
  } else {
    lines(datRet$AccElementFromGeneCutoff , datRet$TargetRelaxRatio, col=rgb(1,0,0,1/20), lwd=1 );
    lines(datRet$AccElementFromGeneCutoff , datRet$BgRelaxRatio, col=rgb(0,0,1,1/20), lwd=1);
    
  }
}

fnPlot(F);

for(i in 1:500) {
  cat(i," ...\n");
  fnPlot(T)
}


