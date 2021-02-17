
sPop <- "PLP.PDC";
nMinFreq <- 0.2; #min derived allele freq to use.

sRNACoord <- "mRNA.coord.txt";
s4DSites <- paste("/beegfs_old/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_PLP_0.05to1/ret_pop_", sPop, "_syn_table.txt" ,sep=""); #these are synonymous sites, not 4D, for 4D , see ../MKTEST
sNonSynSites <- paste("/beegfs_old/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_PLP_0.05to1/ret_pop_", sPop, "_nonsyn_table.txt" ,sep=""); #these are non-syn sites computed based on actual mutations.

sink(file=paste(sPop, ".pergene.alpha.log", sep=""));

cat("Loading gene coordinates...\n");
datRNACoord <- read.table(sRNACoord, header=T, quote='', sep="\t");
datRNACoord <- datRNACoord[with(datRNACoord, order(chr, start)), ]

cat("Loading synonymous sites...\n");
dat4DSites <- read.table(s4DSites , header=F);
cat("Loading non-synonymous sites...\n");
datNonSynSites <- read.table(sNonSynSites , header=F);

sSNPFolder <- paste("../", sPop, sep=''); # <- JUST CHANGE THIS FOLDER FOR DIFFERENT POPS!
sOut <- paste(basename(sSNPFolder) , ".pergene.alpha.txt" , sep="");
sPolySites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.txt", sep=""); #polymorphic sites
sDivergentSites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.invariable_der.txt", sep=""); #divergent sites
cat("Loading polymorphic sites...\n");
datPolySites <- read.table(sPolySites, header=T);
cat("Loading divergent sites...\n");
datDivergentSites <- read.table(sDivergentSites, header=T);

datPolySites <- datPolySites[datPolySites$derived_freq >= nMinFreq, ];

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
  
  dN <- length(intersect(arrGeneNonSynSites , arrGeneDivergentSites));
  dS <- length(intersect(arrGene4DSites , arrGeneDivergentSites));
  pN <- length(intersect(arrGeneNonSynSites , arrGenePolySites));
  pS <- length(intersect(arrGene4DSites , arrGenePolySites));
  
  if (dN < 1 | dS < 1 | pN < 1 ) {
    cat( as.character(datRNACoord[nRow, c('OrthoID' )]), as.character(datRNACoord[nRow, c('GeneSymbol' )]), " have too few sites, skip\n");
    next;
  }
  
  
  oConTable <- matrix(c(dN , dS, pN , pS ), nrow = 2)
  DoS <- dN / (dN+dS) - pN/(pN+pS);

  datOut <- datRNACoord[nRow, ,drop=F]
  datOut$alpha <- 0;
  datOut$p <- 0;
  datOut$dN <- dN;
  datOut$dS <- dS;
  datOut$pN <- pN;
  datOut$pS <- pS;
  datOut$DoS <- DoS;

  if (pS == 0) {
    pS <- 1; #set pS to 1 if it's 0 , so that alpha can be computed. the alpha computed this way should be conservative for looking at relaxation.
  }
  alpha <- 1 - ( (dS * pN) / (dN * pS) );
  oTest <- fisher.test(oConTable);

  datOut$alpha <- alpha;
  datOut$p <- oTest$p.value;

  
  write.table(datOut  , file=sOut, quote=F, sep="\t", row.names = F, col.names = (!bHeaderWritten) , append = bHeaderWritten);
  bHeaderWritten <- T;
}
sink();
