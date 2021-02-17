source("asym.MK.R");

sPop <- "PLP.PDC";


sRNACoord <- "mRNA.coord.txt";
s4DSites <- paste("/beegfs_old/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_PLP_0.05to1/ret_pop_", sPop, "_syn_table.txt" ,sep=""); #these are synonymous sites, not 4D, for 4D , see ../MKTEST
sNonSynSites <- paste("/beegfs_old/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_PLP_0.05to1/ret_pop_", sPop, "_nonsyn_table.txt" ,sep=""); #these are non-syn sites computed based on actual mutations.

cat("Loading gene coordinates...\n");
datRNACoord <- read.table(sRNACoord, header=T, quote='', sep="\t");
datRNACoord <- datRNACoord[with(datRNACoord, order(chr, start)), ]

cat("Loading synonymous sites...\n");
dat4DSites <- read.table(s4DSites , header=F);
cat("Loading non-synonymous sites...\n");
datNonSynSites <- read.table(sNonSynSites , header=F);

sSNPFolder <- paste("../", sPop, sep=''); # <- JUST CHANGE THIS FOLDER FOR DIFFERENT POPS!
sOut <- paste(basename(sSNPFolder) , ".global.alpha.txt" , sep="");
sPolySites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.txt", sep=""); #polymorphic sites
sDivergentSites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.invariable_der.txt", sep=""); #divergent sites
cat("Loading polymorphic sites...\n");
datPolySites <- read.table(sPolySites, header=T);
cat("Loading divergent sites...\n");
datDivergentSites <- read.table(sDivergentSites, header=T);

cat("Computing...\n");
dN <- dS <- pN <- pS <- 0;

arrAlleleFreq <- seq(0.08, 0.90, 0.05); #2% bins, RIGHT BOUNDARY
arrPN <- rep(0, length(arrAlleleFreq) );
arrPS <- rep(0, length(arrAlleleFreq) );
arrScflds <- levels(dat4DSites$V1);
for(sScf in arrScflds) {
  sGeneChr <- sScf;
  if (!grepl("chr", sScf )) { #only use genes assigned onto LGs. to save time
    break;
  }
  cat(sScf,"\n");
  arrGene4DSites <- dat4DSites[dat4DSites$V1 == sGeneChr , 3 , drop=T];
  arrGeneNonSynSites <-  datNonSynSites[datNonSynSites$V1 == sGeneChr , 3 , drop=T];
  datGenePolySites <- datPolySites[datPolySites$scaffold == sGeneChr , c( 'site' , 'derived_freq') , drop=F];
  datGeneDivergentSites <-  datDivergentSites[datDivergentSites$scaffold == sGeneChr  , c( 'site' ) , drop=F];
  
  dN <- dN + length(intersect(arrGeneNonSynSites , datGeneDivergentSites$site));
  dS <- dS + length(intersect(arrGene4DSites , datGeneDivergentSites$site));
  pN <- pN + length(intersect(arrGeneNonSynSites , datGenePolySites$site));
  pS <- pS + length(intersect(arrGene4DSites , datGenePolySites$site));
  
  for(nFreqCount in 1:length(arrAlleleFreq)) {
    nFreqLowCut <- 0;
    if (nFreqCount > 1) {
      nFreqLowCut <- arrAlleleFreq[nFreqCount-1];
    }
    nFreqUpperCut <- arrAlleleFreq[nFreqCount];
    arrGenePolySites <- datGenePolySites[datGenePolySites$derived_freq > nFreqLowCut & datGenePolySites$derived_freq <= nFreqUpperCut , 'site', drop=T ];
    arrPS[nFreqCount] <-  arrPS[nFreqCount] + length(intersect(arrGene4DSites , arrGenePolySites));
    arrPN[nFreqCount] <-  arrPN[nFreqCount] + length(intersect(arrGeneNonSynSites , arrGenePolySites));
  }
  
  
  
}

sink(file = sOut)
alpha <- 1 - ( (dS * pN) / (dN * pS) );
oConTable <- matrix(c(dN , dS, pN , pS ), nrow = 2);
oTest <- fisher.test(oConTable);

cat("dN\tdS\tpN\tpS\n");
cat(dN, "\t", dS, "\t", pN, "\t", pS, "\n");
cat("alpha global = ", alpha);
oTest

arrAsymMK <- asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.15, xhigh=0.95, output='table')
arrAsymMK

n10percIdx <- which(arrAlleleFreq == 0.10);
alpha_10perc <- 1 - ( (dS * arrPN[n10percIdx] ) / (dN * arrPS[n10percIdx] ) );
alpha_asymptotic <- arrAsymMK['alpha_asymptotic'];

alphaDiff <- alpha_asymptotic - alpha_10perc;

cat("alpha_15%\talpha_original\talpha_asymptotic\talphaDiff\n");
cat(alpha_10perc, "\t",  as.numeric(arrAsymMK['alpha_original']) , "\t", as.numeric(alpha_asymptotic), "\t", as.numeric(alphaDiff), "\n");


sink();
pdf(file=paste(sOut, ".pdf",sep="") , width=6, height=6)
asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.2, xhigh=0.9, output='default')
abline(h=alpha_10perc , col='blue');
abline(h= as.numeric(arrAsymMK['alpha_original'])  , col='green');
abline(h= 0  , col='black');

dev.off();
  
