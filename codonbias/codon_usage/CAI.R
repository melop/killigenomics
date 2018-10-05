#this script computes CAI for each gene
#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/codonbias/Nothobranchius");
sCodonUsageCounts <- "codon_use_counts.txt";
#sCodonUsageCounts <- "codon_use_counts_frameshift_1.txt";
#sCodonUsageCounts <- "codon_use_counts_frameshift_2.txt";
#sCodonUsageCounts <- "codon_use_counts_revcomp.txt";
#sCodonUsageCounts <- "codon_use_counts_frameshift_1_revcomp.txt";
#sCodonUsageCounts <- "codon_use_counts_frameshift_2_revcomp.txt";

#definition of a single ref table is no longer allowed. The species specific ref table will be used.
#sRefTable <- "./Codon_usage_highlyexpressed_genes_NOR.txt.RSCU.txt"; #use highly exp genes as reference.
#sRefTable <- "./Codon_usage_frameshift_1_highlyexpressed_genes_NotfurRef.txt.RSCU.txt"; #use highly exp genes as reference.
#sRefTable <- "./Codon_usage_frameshift_2_highlyexpressed_genes_NotfurRef.txt.RSCU.txt"; #use highly exp genes as reference.
#sRefTable <- "./Codon_usage_revcomp_highlyexpressed_genes_NotfurRef.txt.RSCU.txt"; #use highly exp genes as reference.
#sRefTable <- "./Codon_usage_frameshift_1_revcomp_highlyexpressed_genes_NotfurRef.txt.RSCU.txt"; #use highly exp genes as reference.
#sRefTable <- "./Codon_usage_frameshift_2_revcomp_highlyexpressed_genes_NotfurRef.txt.RSCU.txt"; #use highly exp genes as reference.

          
          


nMinCodonLength <- 100; # at least need to have 100 codons to include in analysis
sOut <- "CAI_min_100AA.txt";
#sOut <- "CAI_min_100AA_frameshift_1_highexpgene_ref.txt";
#sOut <- "CAI_min_100AA_frameshift_2_highexpgene_ref.txt";
#sOut <- "CAI_min_100AA_revcomp_highexpgene_ref.txt";
#sOut <- "CAI_min_100AA_frameshift_1_revcomp_highexpgene_ref.txt";
#sOut <- "CAI_min_100AA_frameshift_2_revcomp_highexpgene_ref.txt";

datCodonUsageCounts <- read.table(sCodonUsageCounts , header=T, stringsAsFactors = F);


arrGenes <- unique(datCodonUsageCounts$GeneName);
arrSpecies <- colnames(datCodonUsageCounts);
arrSpecies <- arrSpecies[3:length(arrSpecies)];

lsSpeciesSpecificRefs <- list();
for(sSp in arrSpecies) {
	lsSpeciesSpecificRefs[[sSp]] <- read.table(paste('./Codon_usage_highlyexpressed_genes_', sSp, '.txt.RSCU.txt' ,sep='') , header=T, stringsAsFactors = F); #load specifs specific ref table.
	lsSpeciesSpecificRefs[[sSp]]$lnW <- log(lsSpeciesSpecificRefs[[sSp]]$w);
}



cat("GeneName", paste( rep(arrSpecies , each=2), c("CAI", "AALen"),sep="_") , "\n", file=sOut, sep="\t");

fnComputeCAI <- function(datC, sSp) {
  colnames(datC) <- c("Codon", "X");
  datMatchedCodon <- merge( lsSpeciesSpecificRefs[[sSp]] , datC, by="Codon");
  nL <- sum(datMatchedCodon$X);
  if (nL < nMinCodonLength) {
    return(c(NA , nL ) );
  }
  nCAI <- exp( sum(datMatchedCodon$X * datMatchedCodon$lnW) / nL);
  return(c(nCAI, nL ) );
}
#go over each gene now

for(i in 1:length(arrGenes) ) { #length(arrGenes)
  lsRet <- list();
  sGene <- arrGenes[i];
  lsRet$Gene <- sGene;
  datGeneCodonUse <- datCodonUsageCounts[datCodonUsageCounts$GeneName==sGene,];
  #now go through each species
  for(sSp in arrSpecies) {
   arrCAIRet  <- fnComputeCAI(datGeneCodonUse[, c("Codon", sSp)], sSp);
   lsRet[[paste(sSp , "CAI" , sep="_")]] <- arrCAIRet[1];
   lsRet[[paste(sSp , "AALen" , sep="_")]] <- arrCAIRet[2];
  }
  
  write.table( rbind(NULL  , lsRet) , file=sOut , append = T, sep="\t", quote=F, row.names = F, col.names = F);
  
}

