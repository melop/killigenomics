#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/");

arrGenera <- c('Nothobranchius', 'Callopanchax', 'Aphyosemion', 'Scriptaphyosemion', 'PronothoNothos', 'Fundulopanchax', 'Archiaphyosemion', 'Epiplatys');
arrFiles <- paste("sum_", arrGenera, ".txt", sep="")
for(sFile in arrFiles) {
  datTable <- read.table(sFile, header=F, sep="\t", quote = "", fill=T);
  datTable <- datTable[complete.cases(datTable),];
  datTable$FDR <- p.adjust(datTable$V5, method="fdr");
  datTableTrim <- datTable[, c(1,3,9,4, 5,10)];
  colnames( datTableTrim ) <- c("Ortholog_ID" , "GeneSymbol", "GeneName", "AALen", "P", "FDR");
  datRelax <- datTableTrim[datTableTrim$P <= 0.05 , ];
  sRelaxOut <- paste("codeml_",sFile, sep="");

  write.table(datRelax[with(datRelax, order(P)), ], file=sRelaxOut, col.names = T, row.names = F, quote = F, sep="\t");

}