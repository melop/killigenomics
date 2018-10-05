#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/codonbias/Nothobranchius");

sCodonTable <- "codontable.txt";
#sCodonUsage <- "Codon_usage_highlyexpressed_genes_X10NVS.txt"; #"Codon_usage_highlyexpressed_genes_X10NVS.txt"; #"Codon_usage_highlyexpressed_genes_NotfurRef.txt"; #"functional_wobble.txt";# "Codon_usage_highlyexpressed_genes_NotfurRef.txt" #"functional.txt"; #""; plp_functional.txt# this is codon usage of tRNA, can be replaced by codon usage table from highly expressed genes

arrCodonUsage <- Sys.glob("Codon_usage_*.txt");

arrFolds <- c(2,3,4,5,6); #which fold of redundancy to include

datCodonTable <- read.table(sCodonTable, header=T, stringsAsFactors = F);

lsCounts <- list();

for( nLn in 1:nrow(datCodonTable)) {
  
  sAA <- datCodonTable[nLn,2];
  sCodons <- datCodonTable[nLn,3];
  arrCodons <- unlist( strsplit(sCodons,",") );
  
  if (sAA == "START" || sAA == "STOP") {
    next;
  }
  if (! (length(arrCodons) %in% arrFolds)) {
    next;
  }
  
  lsCounts[[sAA]] <- as.list( rep(0, length(arrCodons)));
  names( lsCounts[[sAA]] ) <- arrCodons;
}

# now read in the codon usage table:
for(sCodonUsage in arrCodonUsage) {
  if (grepl(".RSCU.txt", sCodonUsage)) {
	next;
  }
datCodonUsage <- read.table(sCodonUsage, header=F);

#NOW Fill in the counts
for( sAA in names(lsCounts) ) {

  lsCodons <- lsCounts[[sAA]];
  for (sCodon in names(lsCodons) ) {
    if (! (sCodon %in% datCodonUsage$V1) ) {
      lsCounts[[sAA]][[sCodon]] <- 0.5;
    } else {
      lsCounts[[sAA]][[sCodon]] <- datCodonUsage[datCodonUsage$V1==sCodon, 2];
    }
  }
}

datRSCU <- NULL;
for( sAA in names(lsCounts) ) {
  
  lsCodons <- lsCounts[[sAA]];
  arrCodons <- unlist(lsCodons);
  arrNormalized <- arrCodons / ( sum(arrCodons)/length(arrCodons));
  arrNormalized[arrNormalized==0] <- 0.0001;
  arrW <- arrNormalized / max(arrNormalized);
  datRSCU <- rbind(datRSCU, cbind(sAA , names(lsCodons), arrNormalized, arrW ));
}

colnames(datRSCU) <- c("AminoAcid" , "Codon", "RSCU", "w");

write.table(datRSCU, file=paste(sCodonUsage, ".RSCU.txt", sep=""), sep="\t", quote=F, row.names = F);

}
