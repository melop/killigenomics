setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_above0.2_NVGout/MKTEST_syn_exclCMD");

#use DoS > avg : DoS < avg ratio

sBlackList <- "../../relax_codeml/consurf/alignment_blacklist.txt"; #problems in alignment potentially. exclude.

sAlpha <- "CombinedPops.pergene.alpha.txt";
sAlphaCol <- "avgDoS";#"avgDoS" #"avgalpha";

sOutPDFName <-"correlate_mk_relax.childexc.excludeMNM.avg.plotDoSMean.kvsknonsig.DoSAll.pdf";
#sOutPDFName <-"correlate_mk_relax.childexc.excludeMNM.avg.plotDoSMean.kvsknonsig.DoSlessthan0.pdf";

datBlackList <- read.table(sBlackList, header=F, stringsAsFactors = F);


sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0/sum_Nothobranchius_childrenexcl.txt";
sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp_exclude_MNM/rerun_omega0/sum_Nothobranchius_childrenexcl.txt";



#sCodeML <- "~/killifish_genomes/hyphyrelax/codeml_46_spp/sum_Nothobranchius.txt";
sCodeML <- "~/killifish_genomes/hyphyrelax/codeml_46_spp/rerun_exclude_MNM/sum_Nothobranchius.txt";


datCodeml <- read.table(sCodeML, header=F, quote='', sep="\t" ,fill = TRUE);
arrCodemlSig <- datCodeml[datCodeml$V5 < 0.05, 1]


datMK <- read.table(sAlpha, header=T, quote='', sep="\t");
datK <- read.table(sK, header=F, quote='', sep="\t");

datK <- datK[ !(datK$V1 %in% datBlackList$V1), ];
datMK <- datMK[ !(datMK$OrthoID %in% datBlackList$V1), ];

#View(datMK[ (datMK$OrthoID %in% datBlackList$V1), ]);

datMerge <- merge(datMK, datK, by.x = "OrthoID", by.y = "V1");

arrPCutoffs <- seq(0.01, 0.95, 0.01);
arrFisherPs <- rep(0, length(arrPCutoffs) );
arrWilcoxPs <- rep(0, length(arrPCutoffs) );

arrFisherRatio1 <- arrFisherPs;
arrFisherRatio2 <- arrFisherPs;
arrMeanDoS1 <- arrFisherPs;
arrMeanDoS2 <- arrFisherPs;



arrFisherPs2 <- rep(0, length(arrPCutoffs) );
arrPosRatio1 <- arrPosRatio2 <- arrFisherPs;

for (nPCutoffCount in 1:length( arrPCutoffs) ) {
  nPCutoff <- arrPCutoffs[nPCutoffCount];
   datMergeFilter <- datMerge;
   datMergeFilter <- datMergeFilter[!is.na(datMergeFilter[, sAlphaCol]), ];
  # datMergeFilter <- datMergeFilter[datMergeFilter[, sAlphaCol]<=0, ]; # exclude pos selection
  nMeanAlpha <- mean(datMergeFilter[, sAlphaCol], na.rm = T);
  
  datIntensified <- datMergeFilter[datMergeFilter$V5 >= nPCutoff, ];
  datRelaxed <- datMergeFilter[datMergeFilter$V9 < 1 & datMergeFilter$V5 < nPCutoff , ];


  nIntenseRelax <- length(datIntensified[datIntensified[, sAlphaCol] < nMeanAlpha, 1]);
  nIntenseIntense <- length(datIntensified[datIntensified[, sAlphaCol] > nMeanAlpha, 1]);
  nRelaxIntense <- length(datRelaxed[datRelaxed[, sAlphaCol] > nMeanAlpha , 1]);
  nRelaxRelax <- length(datRelaxed[datRelaxed[, sAlphaCol] < nMeanAlpha, 1]);
  
  nMeanIntenseDoS <- mean(datIntensified[, sAlphaCol]);
  nMeanRelaxDoS<- mean(datRelaxed[, sAlphaCol]);
  
  oTestWil <- wilcox.test(datIntensified[, sAlphaCol], datRelaxed[, sAlphaCol])
  arrWilcoxPs[nPCutoffCount] <- oTestWil$p.value

  
  nRelaxIntensePos <- length(datRelaxed[datRelaxed[, sAlphaCol] > 0 &  datRelaxed$OrthoID %in% arrCodemlSig , 1]); 
  nRelaxIntenseNotPos <- nRelaxIntense - nRelaxIntensePos;
  nRelaxRelaxPos <- length(datRelaxed[datRelaxed[, sAlphaCol] < 0 &  datRelaxed$OrthoID %in% arrCodemlSig , 1]); 
  nRelaxRelaxNotPos <- nRelaxRelax - nRelaxRelaxPos;
  
  
  oMat <- matrix(c(nIntenseRelax, nIntenseIntense, nRelaxRelax, nRelaxIntense), nrow = 2 )
  oTest <- fisher.test(oMat);
  arrFisherPs[nPCutoffCount] <- oTest$p.value;
  # arrFisherRatio1[nPCutoffCount] <- nIntenseRelax / nIntenseIntense;
  # arrFisherRatio2[nPCutoffCount] <- nRelaxRelax / nRelaxIntense;
   arrFisherRatio1[nPCutoffCount] <- 1 / (nIntenseRelax / nIntenseIntense);
   arrFisherRatio2[nPCutoffCount] <- 1/ ( nRelaxRelax / nRelaxIntense);
  
  arrMeanDoS1[nPCutoffCount] <-nMeanIntenseDoS ;
  arrMeanDoS2[nPCutoffCount] <- nMeanRelaxDoS ;
  
  
  oMat <- matrix(c(nRelaxIntensePos, nRelaxIntenseNotPos, nRelaxRelaxPos, nRelaxRelaxNotPos), nrow = 2 )
  oTest <- fisher.test(oMat);
  arrFisherPs2[nPCutoffCount] <- oTest$p.value;
  
  arrPosRatio1[nPCutoffCount] <- nRelaxIntensePos / nRelaxIntenseNotPos;
  arrPosRatio2[nPCutoffCount] <- nRelaxRelaxPos / nRelaxRelaxNotPos;
  
}

pdf(sOutPDFName, width=3.5, height=4);
 plot(arrPCutoffs , arrMeanDoS2, pch=16, col='blue', ylim=c( min(c(arrMeanDoS1, arrMeanDoS2))-0.03,  max(c(arrMeanDoS1, arrMeanDoS2)) ), xlab="Relax test P value cutoff" , ylab="Mean DoS | Wilcox test");
 points(arrPCutoffs , arrMeanDoS1, pch=16, col='black');

par(new=T)
plot(arrPCutoffs , arrFisherPs, col=rgb(0,0,0,0), pch=16, ylim=c(0,1), axes=F , xlab="", ylab="");
points(arrPCutoffs , arrWilcoxPs, col="grey", pch=16);

axis(side = 4)
# legend('topright', # places a legend at the appropriate place c(“Health”,”Defense”), # puts text in the legend
#        
#        lty=c(1,1,1), # gives the legend appropriate symbols (lines)
#        
#        lwd=c(2.5,2.5,2.5),col=c('blue','red', 'grey'), legend = c( paste("K < 1:",sAlphaCol,"> 0 /",sAlphaCol,"< 0"), paste("K > 1:",sAlphaCol,"> 0 /",sAlphaCol,"< 0"), "Fisher's test") ) # gives the legend lines the correct color and width
# 
# #       lwd=c(2.5,2.5,2.5),col=c('blue','red', 'grey'), legend = c( paste("K < 1:",sAlphaCol,"< 0 /",sAlphaCol,"> 0"), paste("K > 1:",sAlphaCol,"< 0 /",sAlphaCol,"> 0"), "Fisher's test") ) # gives the legend lines the correct color and width


dev.off();

