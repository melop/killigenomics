#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0");

arrFiles <- Sys.glob("sum_*[^err].txt");

for(sFile in arrFiles) {
  datTable <- read.table(sFile, header=F, sep="\t", quote = "", fill=T);
  datTable <- datTable[complete.cases(datTable),];
  datTable$FDR <- p.adjust(datTable$V5, method="fdr");
  datTableTrim <- datTable[, c(1,3,23,4, 9, 5,25)];
  colnames( datTableTrim ) <- c("Ortholog_ID" , "GeneSymbol", "GeneName", "AALen", "K", "P", "FDR");
  datRelax <- datTableTrim[datTableTrim$P <= 0.05 & datTableTrim$K < 1, ];
  datIntense <- datTableTrim[datTableTrim$P <= 0.05 & datTableTrim$K > 1, ];
  sRelaxOut <- paste("relax_",sFile, sep="");
  sIntenseOut <- paste("intense_",sFile, sep="");
  
  write.table(datRelax[with(datRelax, order(P)), ], file=sRelaxOut, col.names = T, row.names = F, quote = F, sep="\t");
  write.table(datIntense[with(datIntense, order(P)), ], file=sIntenseOut, col.names = T, row.names = F, quote = F, sep="\t");
  
}