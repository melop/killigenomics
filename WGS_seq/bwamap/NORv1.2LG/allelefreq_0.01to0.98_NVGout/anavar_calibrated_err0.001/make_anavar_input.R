#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_0.01to0.98_NVGout/anavar_calibrated_err0.001/");

nChromosomes <- 50; #set chromosome to 50. this is "n" in the anavar configuration file
sPop <- "RACDRY.excludehybrid";
sCalibration <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1_hicov_err0.001/allelefreq_0.01to0.98_NVGout/accuracy/RACDRY.acc.txt";

#sPop <- "RACWET.excludehybrid";
#sCalibration <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1_hicov_err0.001/allelefreq_0.01to0.98_NVGout/accuracy/RACWET.acc.txt";

#sPop <- "ORTDRY.excludehybrid";
#sCalibration <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1_hicov_err0.001/allelefreq_0.01to0.98_NVGout/accuracy/ORTDRY.acc.txt";

#sPop <- "ORTWET";
#sCalibration <- "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1_hicov_err0.001/allelefreq_0.01to0.98_NVGout/accuracy/ORTWET.acc.txt";


datCalibration <- read.table(sCalibration , header = T);
datCalibration$allele_freq_int <- as.integer(datCalibration$allele_freq * 1e6);

s4DSites <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_cov_filter/ret_pop_", sPop, "_4fold_table.txt" ,sep=""); #these are 4D
sNonSynSites <- paste("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/classify_sites_cov_filter/ret_pop_", sPop, "_nonsyn_table.txt" ,sep=""); #these are non-syn sites computed based on actual mutations.

cat("Loading synonymous sites...\n");
dat4DSites <- read.table(s4DSites , header=F);
cat("Loading non-synonymous sites...\n");
datNonSynSites <- read.table(sNonSynSites , header=F);

sSNPFolder <- paste("../", sPop, sep=''); # <- JUST CHANGE THIS FOLDER FOR DIFFERENT POPS!
sOut <- paste(basename(sSNPFolder) , ".SFS.for.anavar.txt" , sep="");
sPolySites <- paste(sSNPFolder , "/ancestral_polarized_for_allelefreq.txt", sep=""); #polymorphic sites

cat("Loading polymorphic sites...\n");
datPolySites <- read.table(sPolySites, header=T);

cat("Computing...\n");


datPolySites$PosID <- paste(datPolySites$scaffold, datPolySites$site, sep = "_");
dat4DSites$PosID <- paste(dat4DSites$V1, dat4DSites$V3, sep="_");
datNonSynSites$PosID <- paste(datNonSynSites$V1, datNonSynSites$V3, sep="_");

arr4DPoly <- as.numeric(datPolySites[datPolySites$PosID %in% dat4DSites$PosID, 'derived_freq' ]);
arrNonSynPoly <- as.numeric(datPolySites[datPolySites$PosID %in% datNonSynSites$PosID, 'derived_freq' ]);

nMinFreq <- min(min(arr4DPoly) , min(arrNonSynPoly));
nMaxFreq <- max(max(arr4DPoly) , max(arrNonSynPoly));
arrBreaks <- seq(0, 1-1/nChromosomes , 1/nChromosomes );
if (nMinFreq - arrBreaks[2] > 0) {
  arrBreaks <- arrBreaks + (nMinFreq - arrBreaks[2]+0.0005 );
}
arrBreaks[length(arrBreaks)] <- 1;

#recalculate the calibration factors for each break;

arrCalibrationFactors <- c();
arrCacheSum <- NULL;
for(i in 2:length(arrBreaks)) {
  nLower <- arrBreaks[i-1];
  datSelRows <- data.frame();
  nPlusI <- 0;
  while(nrow(datSelRows)==0) {
    nUpper <- arrBreaks[i+ nPlusI];
    nPlusI <-  nPlusI + 1;
    datSelRows <- datCalibration[datCalibration$allele_freq > nLower & datCalibration$allele_freq <= nUpper, ];
  }
  if (nrow(datSelRows) == 0) {
    arrSums <- arrCacheSum;
  } else {
    arrCacheSum <- arrSums <- colSums(datSelRows[, c('total_called' , 'total_real')])
  }
  nCalibrationFactor <- arrSums[2]/arrSums[1];
  arrCalibrationFactors <- c(arrCalibrationFactors, nCalibrationFactor);
}

pdf(file=paste(sOut, ".pdf",sep=""), width=7,height=5);
hist0d <- hist( arrNonSynPoly, breaks=arrBreaks, col=rgb(0,0,1,1/3), main=sPop, xlab="Derived Allele Freq");
hist4d <- hist( arr4DPoly, breaks=arrBreaks, col=rgb(1,0,0,1/3), add=TRUE);

hist0d$counts <-  hist0d$counts * arrCalibrationFactors
hist4d$counts <-  hist4d$counts * arrCalibrationFactors

hist0d$density <- hist0d$counts / sum(hist0d$counts);
hist4d$density <- hist4d$counts / sum(hist4d$counts);

plot(hist0d, col=rgb(0,0,1,1/3), main=sPop, xlab="Derived Allele Freq");
plot(hist4d, col=rgb(1,0,0,1/3), add=TRUE);

dev.off();
length(hist0d$counts)
length(hist4d$counts)
sink(file=sOut);#output
cat("model: neutralSNP_vs_selectedSNP\n");
cat("n: ", nChromosomes,"\n");
cat("folded: false\n"); 
cat("r_range: 0.02, 50\n");
cat("neu_m: ", length(dat4DSites$PosID), "\n");
cat("neu_sfs: ", paste(hist4d$counts, collapse = ", "), "\n");
cat("neu_theta_range: 1e-5, 0.1\n");
cat("neu_e_range: 0, 0.5\n");
cat("sel_m: " , length(datNonSynSites$PosID), "\n");
cat("sel_sfs: ", paste(hist0d$counts, collapse = ", "), "\n"); 
sink();
