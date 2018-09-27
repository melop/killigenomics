setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/MKTEST_syn");
sPop <- "RACWET.excludehybrid";
#sPop <- "ORTDRY.excludehybrid";
#sPop <- "ORTWET";
#sPop <- "RACDRY.excludehybrid";


sRNACoord <- "mRNA.coord.txt";
s4DSites <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites/ret_4fold_table.txt" ,sep=""); #these are synonymous sites, not 4D, for 4D , see ../MKTEST
sNonSynSites <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites/ret_nonsyn_table.txt" ,sep=""); #these are non-syn sites computed based on actual mutations.

cat("Loading gene coordinates...\n");
datRNACoord <- read.table(sRNACoord, header=T, quote='', sep="\t");
datRNACoord <- datRNACoord[with(datRNACoord, order(chr, start)), ]

cat("Loading synonymous sites...\n");
dat4DSites <- read.table(s4DSites , header=F);
cat("Loading non-synonymous sites...\n");
datNonSynSites <- read.table(sNonSynSites , header=F);

sSNPFolder <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/", sPop, sep=''); # <- JUST CHANGE THIS FOLDER FOR DIFFERENT POPS!
sOut <- paste(basename(sSNPFolder) , ".pergene.dNdSpNpS.txt" , sep="");
sPolySites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.txt", sep=""); #polymorphic sites
sDivergentSites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.invariable_der.txt", sep=""); #divergent sites
cat("Loading polymorphic sites...\n");
datPolySites <- read.table(sPolySites, header=T);
cat("Loading divergent sites...\n");
datDivergentSites <- read.table(sDivergentSites, header=T);

cat("Computing...\n");

bHeaderWritten <- F;
#go over each gene
for(nRow in 1:nrow(datRNACoord)) {
  #nRow <- 10;
  sGeneChr <- as.character( datRNACoord[nRow, 'chr']);
  nGeneStart <- datRNACoord[nRow, 'start'];
  nGeneEnd <- datRNACoord[nRow, 'end'];
  #select sites
  
  arrGene4DSites <- dat4DSites[dat4DSites$V1 == sGeneChr & dat4DSites$V3 >= nGeneStart & dat4DSites$V3 <= nGeneEnd , 3 , drop=T];
  arrGeneNonSynSites <-  datNonSynSites[datNonSynSites$V1 == sGeneChr & datNonSynSites$V3 >= nGeneStart & datNonSynSites$V3 <= nGeneEnd , 3 , drop=T];
  arrGenePolySites <- datPolySites[datPolySites$scaffold == sGeneChr & datPolySites$site >= nGeneStart & datPolySites$site <= nGeneEnd , 2 , drop=T];
  arrGeneDivergentSites <-  datDivergentSites[datDivergentSites$scaffold == sGeneChr & datDivergentSites$site >= nGeneStart & datDivergentSites$site <= nGeneEnd , 2 , drop=T];
  
  dN <- length(intersect(arrGeneNonSynSites , arrGeneDivergentSites))/length(arrGeneNonSynSites);
  dS <- length(intersect(arrGene4DSites , arrGeneDivergentSites))/length(arrGene4DSites);
  pN <- length(intersect(arrGeneNonSynSites , arrGenePolySites))/length(arrGeneNonSynSites);
  pS <- length(intersect(arrGene4DSites , arrGenePolySites))/length(arrGene4DSites);
  
  
  datOut <- datRNACoord[nRow, ,drop=F]
  
  datOut$dN <- dN;
  datOut$dS <- dS;
  datOut$pN <- pN;
  datOut$pS <- pS;
 
  write.table(datOut  , file=sOut, quote=F, sep="\t", row.names = F, col.names = (!bHeaderWritten) , append = bHeaderWritten);
  bHeaderWritten <- T;
}
