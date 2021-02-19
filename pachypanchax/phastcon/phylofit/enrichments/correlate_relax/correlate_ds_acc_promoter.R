setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/correlate_relax");

arrDistance <- c(3, seq(10, 100, 10) );
#arrDistance <- seq(0, 50, 10) ;

sK <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/PML_PLP_dNdS/dNdS_PML_PLP/ret_part0.of.1.txt"

datK <- read.table(sK , header=F, sep="\t" , quote="", fill = T, stringsAsFactors = F);
colnames(datK) <- c("OrthoID", "Success", "GName", "len", "t", 'S', 'N', 'dNdS', 'dN', 'dS', 'tmpdir','fullname');
#datK <- datK[datK$Success == 'Success',];
#arrRelax <- datK[datK$dS < 0 , 'OrthoID'];
#arrNotRelax <- datK[datK$dS >= 0 , 'OrthoID'];


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
    
    datBg <- merge(datBg, datK[ , c('OrthoID', 'dS') ], by.x = 'V1', by.y = 'OrthoID');
    datTarget <- merge(datTarget, datK[ , c('OrthoID', 'dS') ], by.x = 'V1', by.y = 'OrthoID');

    nAvgBGdS <- median(datBg$dS);
    nAvgTargetdS<- median(datTarget$dS);
    
    oTest <- wilcox.test(datBg$dS, datTarget$dS)
    datRow <- data.frame(AccElementFromGeneCutoff=nDist, fisherP=oTest$p.value, BgRelaxRatio=nAvgBGdS, TargetRelaxRatio=nAvgTargetdS );
    datRet <- rbind(datRet, datRow);
    
    # arrKTarget <- datK$K_q_free[datK$GeneName %in% datTarget$V1];
    # arrKBg <- datK$K_q_free[datK$GeneName %in% datBg$V1];
    # median(arrKBg);
    # median(arrKTarget);
    # wilcox.test(arrKTarget, arrKBg)
    
  }
  
  if (!bBootstrap) {
    print(datRet);
    plot(datRet$AccElementFromGeneCutoff , datRet$TargetRelaxRatio, col='red', type = 'l', lwd=3, ylim=c(0.1,0.13), xlab="Distance from start codon", ylab="Median dS"  );
    lines(datRet$AccElementFromGeneCutoff , datRet$BgRelaxRatio, col='blue', lwd=3)
    #View(datRet);
    datRetSig <- datRet[datRet$fisherP<0.05,]
    points(datRetSig$AccElementFromGeneCutoff , datRetSig$TargetRelaxRatio - 0.01, pch=8, col='red');
  } else {
    lines(datRet$AccElementFromGeneCutoff , datRet$TargetRelaxRatio, col=rgb(1,0,0,1/20), lwd=1 );
    lines(datRet$AccElementFromGeneCutoff , datRet$BgRelaxRatio, col=rgb(0,0,1,1/20), lwd=1);
    
  }
}

pdf(file="corr.dS.acc.pdf", width=10, height=6)
fnPlot(F);

for(i in 1:500) {
  cat(i," ...\n");
  fnPlot(T)
}


dev.off();