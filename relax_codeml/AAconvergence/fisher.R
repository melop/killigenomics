
#sCodeMLRet <- "relax_codeml/codeml_46_spp";
#sRelaxRet <- "relax_codeml/relax_46_spp/rerun_omega0";

sCodeMLRet <- "relax_codeml/codeml_46_spp/rerun_exclude_MNM";
sRelaxRet <- "relax_codeml/relax_46_spp_exclude_MNM/rerun_omega0";


#sDir <- "./annotate_aa_changes/" ; 

sDir <- "./annotate_aa_changes_excludeCMD/" ; #The folder of annotate_aa_changes output

arrGenera <- c("Nothobranchius", "Callopanchax");
arrAbb <- c("Nothos", "Callo");


datOrthologs <- read.table("/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/map2humanensemblID/killi_orthologs_table.ensemblinformed.humanEns.txt" , sep="\t", quote="", fill = T, header=T);

datConvergentAA <- read.table(paste(sDir,"/ret_part0.of.1.txt", sep="" ), header=F, sep="\t", quote='');
datConvergentAAMerge <- merge(datConvergentAA, datSiteLK, by = c('V1', 'V2') )
colnames(datConvergentAAMerge) <- c('OrthoID','Pos','ConvergentAA', 'ConvergentTaxa', 'OutgroupTaxa', 'SpeciesTreeLK', 'ForcedMonophylyLK','DeltaLK');
datConvergentAAMerge <- merge(datConvergentAAMerge , datOrthologs[, c(1,2,3, ncol(datOrthologs) )], by.x="OrthoID", by.y = 'Group_Id', all.x = T);


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


datConvergentAAMerge <- merge(datConvergentAAMerge, datCodeml[, c(1,5)], by.x="OrthoID", by.y="V1", all.x = T);
datConvergentAAMerge <- merge(datConvergentAAMerge, datCodeml2[, c(1,5)], by.x="OrthoID", by.y="V1", all.x = T);
datConvergentAAMerge <- merge(datConvergentAAMerge, datRelax[, c(1,5,9)], by.x="OrthoID", by.y="V1", all.x = T);
datConvergentAAMerge <- merge(datConvergentAAMerge, datRelax2[, c(1,5,9)], by.x="OrthoID", by.y="V1", all.x = T);



write.table(datConvergentAAMerge, paste(sDir,"/annotatedAA.withLK.withCodemlRelax.sig.txt",sep=""), col.names = T, row.names = F, quote=F, sep="\t" );

sink(file = paste(sDir,"/fisher.test.txt",sep="") );
cat( length(datConvergentAAMerge[,'OrthoID']) , " sites convergent out of ");
cat( length( unique(datConvergentAAMerge[,'OrthoID'])) , " genes (before overlap with codeml/relax results)\n\n ");

datConvergentAAMerge <- datConvergentAAMerge[complete.cases(datConvergentAAMerge), ];

cat( length(datConvergentAAMerge[,'OrthoID']) , " sites convergent out of ");
cat( length( unique(datConvergentAAMerge[,'OrthoID'])) , " genes\n\n ");

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

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]<0.05 & datConvergentAAMerge[sCodemlCol2]<0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]>=0.05 | datConvergentAAMerge[sCodemlCol2]>=0.05) , 'OrthoID']));
nPosNonConvergeAA <-  length(arrCodemlBothSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlEitherNonSig) - nNotPosConvergeAA;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )

cat("###################### relaxed sel in both genera ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]<0.05 & datConvergentAAMerge[sRelaxKCol1]< 1 & datConvergentAAMerge[sCodemlCol1]>=0.05 )  & (datConvergentAAMerge[sRelaxPCol2]<0.05 & datConvergentAAMerge[sRelaxKCol2]< 1 & datConvergentAAMerge[sCodemlCol2]>=0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]>=0.05 )  | (datConvergentAAMerge[sRelaxPCol2]>=0.05 )  , 'OrthoID']));
nPosNonConvergeAA <- length(arrRelaxBothSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxEitherNonSig) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )

cat("###################### Pos sel in either genera ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]<0.05 | datConvergentAAMerge[sCodemlCol2]<0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]>=0.05 & datConvergentAAMerge[sCodemlCol2]>=0.05) , 'OrthoID']));
nPosNonConvergeAA <- length(arrCodemlEitherSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlBothNonSig) - nNotPosConvergeAA;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )


cat("###################### relaxed sel in either genera ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]<0.05 & datConvergentAAMerge[sRelaxKCol1]< 1 & datConvergentAAMerge[sCodemlCol1]>=0.05 )  | (datConvergentAAMerge[sRelaxPCol2]<0.05 & datConvergentAAMerge[sRelaxKCol2]< 1 & datConvergentAAMerge[sCodemlCol2]>=0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]>=0.05)  & (datConvergentAAMerge[sRelaxPCol2]>=0.05  )  , 'OrthoID']));
nPosNonConvergeAA <- length(arrRelaxEitherSig) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxBothNonSig) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )




cat("###################### Pos sel in genus 1 ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]<0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol1]>=0.05 ) , 'OrthoID']));
nPosNonConvergeAA <- length(arrCodemlSig1) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlNonSig1) - nNotPosConvergeAA;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )

cat("###################### Pos sel in genus 2 ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol2]<0.05 ) , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sCodemlCol2]>=0.05 ) , 'OrthoID']));
nPosNonConvergeAA <- length(arrCodemlSig2) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrCodemlNonSig2) - nNotPosConvergeAA;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Pos' , 'NotPos');
matTest
fisher.test(matTest )


cat("###################### relaxed sel in genus 1 ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]<0.05 & datConvergentAAMerge[sRelaxKCol1]< 1 & datConvergentAAMerge[sCodemlCol1]>=0.05 )   , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol1]>=0.05 )   , 'OrthoID']));
nPosNonConvergeAA <- length(arrRelaxSig1) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxNonSig1) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )

cat("###################### relaxed sel in genus 2 ######################\n");

nPosConvergeAA <- length(unique( datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol2]<0.05 & datConvergentAAMerge[sRelaxKCol2]< 1 & datConvergentAAMerge[sCodemlCol2]>=0.05 )   , 'OrthoID']));
nNotPosConvergeAA <- length(unique(datConvergentAAMerge[(datConvergentAAMerge[sRelaxPCol2]>=0.05 )   , 'OrthoID']));
nPosNonConvergeAA <- length(arrRelaxSig2) - nPosConvergeAA;
nNotPosNonConvergeAA <- length(arrRelaxNonSig2) - nNotPosConvergeAA ;

matTest <- matrix(c(nPosConvergeAA, nNotPosConvergeAA, nPosNonConvergeAA, nNotPosNonConvergeAA) , nrow = 2 );
colnames(matTest) <- c('ConvergeAA' , 'NoConvergeAA');
rownames(matTest) <- c('Relax' , 'NotRelax');
matTest
fisher.test(matTest )




sink();
