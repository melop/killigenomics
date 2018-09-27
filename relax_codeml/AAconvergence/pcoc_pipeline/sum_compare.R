
nMinPosteriorProb <- 0.99;

sCodeMLRet <- "../codeml_46_spp";
sRelaxRet <- "../relax_46_spp/rerun_omega0";
sDir <- "./pcoc_out/" ; #real data
arrGenera <- c("Nothobranchius", "Callopanchax");
arrAbb <- c("Nothos", "Callo");

  # sDir <- "./pcoc_out_Aphy_Scriptaphy/"
  # arrGenera <- c("Aphyosemion", "Scriptaphyosemion");
  # arrAbb <- c("Aphy", "Scriptaphy");

sDir <- "./pcoc_codemlsim_out/" #simulated
arrGenera <- c("Nothobranchius", "Callopanchax");
arrAbb <- c("Nothos", "Callo");
sCodeMLRet <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline/sim_codeml_relax_files/codeml/";
sRelaxRet <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline/sim_codeml_relax_files/relax/";


sDir <- "../../AAconvergence_sim_noposselection/pcoc_pipeline/pcoc_codemlsim_out/" #simulated
arrGenera <- c("Nothobranchius", "Callopanchax");
arrAbb <- c("Nothos", "Callo");
sCodeMLRet <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline/sim_codeml_relax_files/codeml/";
sRelaxRet <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline/sim_codeml_relax_files/relax/";



pConn <- pipe( paste("zcat -f ", sDir, "/*.gz | grep -v OrthoId", sep=""), 'r' ) ;

datPcoc <- read.table(pConn, header =F);
colnames(datPcoc) <- c( 'OrthoID', 'AASite',  'Indel_prop',      'Indel_prop_ConvLeaves',   'PCOC'  ,  'PC'   ,   'OC'  ,    'Path');

datOrthologs <- read.table("/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/map2humanensemblID/killi_orthologs_table.ensemblinformed.humanEns.txt" , sep="\t", quote="", fill = T, header=T);
datPcoc <- merge(datPcoc , datOrthologs[, c(1,2,3, ncol(datOrthologs) )], by.x="OrthoID", by.y = 'Group_Id', all.x = T);

datCodeml <- read.table( paste(sCodeMLRet, "/sum_", arrGenera[1],".txt",sep="" ), quote = '', header = F, sep="\t", fill=T)
datCodeml2 <- read.table(paste(sCodeMLRet, "/sum_", arrGenera[2],".txt",sep="" ), quote = '', header = F, sep="\t", fill=T)

colnames(datCodeml)[5] <- sCodemlCol1 <- paste(arrAbb[1],".codeml.p" ,sep="");
colnames(datCodeml2)[5] <- sCodemlCol2 <- paste(arrAbb[2],".codeml.p" ,sep="");

datRelax <- read.table( paste(sRelaxRet, "/sum_", arrGenera[1],".txt", sep=''), quote = '', header = F, sep="\t", fill=T)
datRelax2 <- read.table( paste(sRelaxRet, "/sum_", arrGenera[2],".txt" , sep=''), quote = '', header = F, sep="\t", fill=T)
colnames(datRelax)[5] <- sRelaxPCol1 <- paste(arrAbb[1],".relax.p" ,sep="");
colnames(datRelax2)[5] <- sRelaxPCol2 <- paste(arrAbb[2],".relax.p" ,sep="");
colnames(datRelax)[9] <- sRelaxKCol1 <- paste(arrAbb[1],".relax.k" ,sep="");
colnames(datRelax2)[9] <- sRelaxKCol2 <-paste(arrAbb[2],".relax.k" ,sep="");


datPcoc <- merge(datPcoc, datCodeml[, c(1,5)], by.x="OrthoID", by.y="V1", all.x = T);
datPcoc <- merge(datPcoc, datCodeml2[, c(1,5)], by.x="OrthoID", by.y="V1", all.x = T);
datPcoc <- merge(datPcoc, datRelax[, c(1,5,9)], by.x="OrthoID", by.y="V1", all.x = T);
datPcoc <- merge(datPcoc, datRelax2[, c(1,5,9)], by.x="OrthoID", by.y="V1", all.x = T);



write.table(datPcoc[datPcoc$PCOC>=nMinPosteriorProb,], paste(sDir,"/pcoc.sig.txt",sep=""), col.names = T, row.names = F, quote=F, sep="\t" );

sink(file = paste(sDir,"/fisher.test.txt",sep="") );
cat( length(datPcoc[datPcoc$PCOC>=nMinPosteriorProb,'OrthoID']) , " sites convergent out of ");
cat( length( unique(datPcoc[datPcoc$PCOC>=nMinPosteriorProb,'OrthoID'])) , " genes (before overlap with codeml/relax results)\n\n ");

datPcoc <- datPcoc[complete.cases(datPcoc), ];

cat( length(datPcoc[datPcoc$PCOC>=nMinPosteriorProb,'OrthoID']) , " sites convergent out of ");
cat( length( unique(datPcoc[datPcoc$PCOC>=nMinPosteriorProb,'OrthoID'])) , " genes\n\n ");

arrPassGenes <- intersect( intersect(datCodeml$V1, datCodeml2$V1), intersect(datRelax$V1, datRelax2$V1) );
nTotalPassGenes <- length(arrPassGenes)
cat("Total complete cases: (genes)" , nTotalPassGenes, "\n");
arrCodemlSig1 <- intersect(datCodeml[datCodeml[sCodemlCol1] < 0.05 , 'V1' ], arrPassGenes)
arrCodemlNonSig1 <- intersect(datCodeml[datCodeml[sCodemlCol1] >= 0.05 , 'V1' ], arrPassGenes)
arrCodemlSig2 <- intersect(datCodeml2[datCodeml2[sCodemlCol2] < 0.05 , 'V1' ], arrPassGenes)
arrCodemlNonSig2 <- intersect(datCodeml2[datCodeml2[sCodemlCol2] >= 0.05 , 'V1' ], arrPassGenes)
arrCodemlBothSig <- intersect(arrCodemlSig1, arrCodemlSig2)
arrCodemlEitherSig <- union(arrCodemlSig1, arrCodemlSig2);
arrCodemlBothNonSig <- intersect(arrCodemlNonSig1, arrCodemlNonSig2)
arrCodemlEitherNonSig <- union(arrCodemlNonSig1, arrCodemlNonSig2)


arrRelaxSig1 <- intersect(datRelax[datRelax[sRelaxPCol1] < 0.05 & datRelax[sRelaxKCol1] < 1 , 'V1' ], arrPassGenes)
arrRelaxSig1 <- arrRelaxSig1[!(arrRelaxSig1 %in% arrCodemlSig1) ]; #remove anything sig in codeml
arrRelaxNonSig1 <- intersect(datRelax[datRelax[sRelaxPCol1] >= 0.05  , 'V1' ], arrPassGenes)
arrRelaxSig2 <- intersect(datRelax2[datRelax2[sRelaxPCol2] < 0.05 & datRelax2[sRelaxKCol2] < 1 , 'V1' ], arrPassGenes)
arrRelaxSig2 <- arrRelaxSig2[!(arrRelaxSig2 %in% arrCodemlSig2) ]; #remove anything sig in codeml
arrRelaxNonSig2 <- intersect(datRelax2[datRelax2[sRelaxPCol2] >= 0.05  , 'V1' ], arrPassGenes)
arrRelaxBothSig <- intersect(arrRelaxSig1, arrRelaxSig2)
arrRelaxEitherSig <- union(arrRelaxSig1, arrRelaxSig2)
arrRelaxBothNonSig <- intersect(arrRelaxNonSig1, arrRelaxNonSig2)
arrRelaxEitherNonSig <- union(arrRelaxNonSig1, arrRelaxNonSig2)



cat("###################### Pos sel in both genera ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 & datPcoc[sCodemlCol2]<0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 | datPcoc[sCodemlCol2]>=0.05) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 & datPcoc[sCodemlCol2]<0.05 ), 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 & datPcoc[sCodemlCol2]>=0.05), 'OrthoID'])) - nNotPosConvergeAA;
nPosNonConvergeAA <-  length(arrCodemlBothSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlEitherNonSig) - nNotPosConvergeAA;


matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )

cat("###################### relaxed sel in both genera ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )  & (datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[ (datPcoc[sRelaxPCol1]>=0.05   | datPcoc[sRelaxPCol2]>=0.05 )  & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )  & (datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 ) , 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol1]>=0.05 & datPcoc[sCodemlCol1]>=0.05 )  & (datPcoc[sRelaxPCol2]>=0.05 & datPcoc[sCodemlCol2]>=0.05 ), 'OrthoID'])) - nNotPosConvergeAA ;
nPosNonConvergeAA <- length(arrRelaxBothSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxEitherNonSig) - nNotPosConvergeAA ;


matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )

cat("###################### Pos sel in either genera ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 | datPcoc[sCodemlCol2]<0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 & datPcoc[sCodemlCol2]>=0.05) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 | datPcoc[sCodemlCol2]<0.05 ), 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 | datPcoc[sCodemlCol2]>=0.05), 'OrthoID'])) - nNotPosConvergeAA;
nPosNonConvergeAA <- length(arrCodemlEitherSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlBothNonSig) - nNotPosConvergeAA;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )


cat("###################### relaxed sel in either genera ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )  | (datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol1]>=0.05  & datPcoc[sRelaxPCol2]>=0.05  )  & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )  | (datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 ) , 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol1]>=0.05 & datPcoc[sCodemlCol1]>=0.05 )  & (datPcoc[sRelaxPCol2]>=0.05 & datPcoc[sCodemlCol2]>=0.05 ), 'OrthoID'])) - nNotPosConvergeAA ;
nPosNonConvergeAA <- length(arrRelaxEitherSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxBothNonSig) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )




cat("###################### Pos sel in genus 1 ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol1]<0.05 ), 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol1]>=0.05 ), 'OrthoID'])) - nNotPosConvergeAA;
nPosNonConvergeAA <- length(arrCodemlSig1) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlNonSig1) - nNotPosConvergeAA;


matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )

cat("###################### Pos sel in genus 2 ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol2]<0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol2]>=0.05 ) & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sCodemlCol2]<0.05 ), 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sCodemlCol2]>=0.05 ), 'OrthoID'])) - nNotPosConvergeAA;
nPosNonConvergeAA <- length(arrCodemlSig2) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlNonSig2) - nNotPosConvergeAA;


matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )


cat("###################### relaxed sel in genus 1 ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )   & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol1]>=0.05 )   & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol1]<0.05 & datPcoc[sRelaxKCol1]< 1 & datPcoc[sCodemlCol1]>=0.05 )  , 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol1]>=0.05 & datPcoc[sCodemlCol1]>=0.05 ) , 'OrthoID'])) - nNotPosConvergeAA ;
nPosNonConvergeAA <- length(arrRelaxSig1) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxNonSig1) - nNotPosConvergeAA ;


matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )

cat("###################### relaxed sel in genus 2 ######################\n");

nPosConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 )   & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
nNotPosConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol2]>=0.05  )   & datPcoc$PCOC >= nMinPosteriorProb, 'OrthoID']));
#nPosNonConvergeAA <- length(unique( datPcoc[(datPcoc[sRelaxPCol2]<0.05 & datPcoc[sRelaxKCol2]< 1 & datPcoc[sCodemlCol2]>=0.05 )  , 'OrthoID'])) - nPosConvergeAA;
#nNotPosNonConvergeAA <- length(unique(datPcoc[(datPcoc[sRelaxPCol2]>=0.05 & datPcoc[sCodemlCol2]>=0.05 ) , 'OrthoID'])) - nNotPosConvergeAA ;
nPosNonConvergeAA <- length(arrRelaxSig2) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxNonSig2) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )




sink();
