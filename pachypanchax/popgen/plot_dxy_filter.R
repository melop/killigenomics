setwd("/beegfs/group_dv/home/RCui/killifish_genomes/plp_dxy");

pIn <- pipe("zcat -f PLP_merged.Fst.Dxy.pi.filter.csv.gz");
dat <- read.csv(pIn, header=T)
dat <- dat[grepl('chr', dat$scaffold),]; #only use assigned scaffolds

arrMedian <- apply(dat[, 6:ncol(dat)], 2,  function(x) {median(x,na.rm = T)} )

arrMean <- apply(dat[, 6:ncol(dat)], 2,  function(x) {mean(x,na.rm = T)} )

arrPopOrder <- c('A', 'B', 'C', 'F', 'D',  'E', 'G');
arrPops <- paste("PLP", arrPopOrder, sep="");

matOutput <- matrix(nrow=length(arrPops), ncol = length(arrPops) )
colnames(matOutput) <- arrPops;
rownames(matOutput) <- arrPops;

for (nPop1 in 1:length(arrPops) ) {
  sPop1 <- arrPops[nPop1];
  for(nPop2 in nPop1:length(arrPops)) {
    sPop2 <- arrPops[nPop2];
    if (sPop1 == sPop2) {
      matOutput[sPop1, sPop2] <- arrMean[paste("pi_", sPop1, sep="")];
    } else {
      matOutput[sPop1, sPop2] <- arrMean[paste("Fst_", sPop1, "_", sPop2, sep="")];
      if (is.na(matOutput[sPop1, sPop2])) {
        matOutput[sPop1, sPop2] <- arrMean[paste("Fst_", sPop2, "_", sPop1, sep="")];
        
      }
      matOutput[sPop2, sPop1] <- arrMean[paste("dxy_", sPop1, "_", sPop2, sep="")];
      if (is.na(matOutput[sPop2, sPop1])) {
        matOutput[sPop2, sPop1] <- arrMean[paste("dxy_", sPop2, "_", sPop1, sep="")];
        
      }
    }
  }
}

write.table(matOutput, file="fst_dxy_matrix.filter.txt", sep = "\t", col.names = T, row.names = T, quote=F);
