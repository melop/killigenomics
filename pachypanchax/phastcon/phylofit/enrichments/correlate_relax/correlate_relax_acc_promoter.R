setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment/correlate_relax");
nPRelaxCutoff <- 0.05;

arrDistance <- c(3, seq(10, 100, 10) );
#arrDistance <- seq(0, 50, 10) ;

sK <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/relaxq_46_spp/sum_Pachypanchax.txt"

datK <- read.table(pipe(paste("{ head -n1 ", sK," ; grep -v GeneName ", sK,"; }", sep='')) , header=T, sep="\t" , quote="", fill = T);
datK <- datK[datK$Success == 'Success',];
arrRelax <- datK[datK$P_relax < nPRelaxCutoff & datK$K_q_free < 1 , 'GeneName'];
arrNotRelax <- datK[datK$P_relax >= nPRelaxCutoff, 'GeneName'];


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
    

    
    nBgRelax <- length(intersect(datBg$V1 , arrRelax));
    nBgNotRelax <- length(intersect(datBg$V1 , arrNotRelax));
    nTargetRelax <- length(intersect(datTarget$V1 , arrRelax));
    nTargetNotRelax <- length(intersect(datTarget$V1 , arrNotRelax));
    
    if (!bBootstrap) {
      write.table(datTarget[datTarget$V1 %in% intersect(datTarget$V1 , arrRelax), ], file = paste("relax_and_loss_of_regelement_",nDist,"_kb.txt", sep='' ) , sep="\t", col.names = T, row.names = F, quote=F);
    }
    
    nBgRelaxRatio <- nBgRelax / (nBgRelax + nBgNotRelax)
    nTargetRelaxRatio <- nTargetRelax / (nTargetRelax+ nTargetNotRelax);
    matTest <- matrix(c(nBgRelax, nBgNotRelax, nTargetRelax, nTargetNotRelax) , nrow=2);
    #matTest
    oTest <- fisher.test(matTest)
    datRow <- data.frame(AccElementFromGeneCutoff=nDist, fisherP=oTest$p.value, BgRelaxRatio=nBgRelaxRatio, BgRelax=nBgRelax, BgNotRelax=nBgNotRelax, TargetRelaxRatio=nTargetRelaxRatio , TargetRelax=nTargetRelax, TargetNotRelax=nTargetNotRelax);
    datRet <- rbind(datRet, datRow);
    
    # arrKTarget <- datK$K_q_free[datK$GeneName %in% datTarget$V1];
    # arrKBg <- datK$K_q_free[datK$GeneName %in% datBg$V1];
    # median(arrKBg);
    # median(arrKTarget);
    # wilcox.test(arrKTarget, arrKBg)
    
  }
  
  if (!bBootstrap) {
    plot(datRet$AccElementFromGeneCutoff , datRet$TargetRelaxRatio, col='red', type = 'l', lwd=3, ylim=c(0.03, 0.07), xlab="Distance from start codon", ylab="Proportion of relaxed genes"  );
    lines(datRet$AccElementFromGeneCutoff , datRet$BgRelaxRatio, col='blue', lwd=3)
    datRetSig <- datRet[datRet$fisherP<0.05,]
    points(datRetSig$AccElementFromGeneCutoff , datRetSig$TargetRelaxRatio + 0.01, pch=8, col='red');
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


