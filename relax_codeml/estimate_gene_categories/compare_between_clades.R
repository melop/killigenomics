setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/estimate_gene_categories/");

sGenus1 <- "Nothobranchius"
sGenus2 <- "Aphyosemion"

#sGenus1 <- "Callopanchax"
#sGenus2 <- "Scriptaphyosemion"



datResampled1 <- read.table(paste("estimated_gene_categories_resampled_", sGenus1, ".txt", sep=""), header=T);
datResampled2 <- read.table(paste("estimated_gene_categories_resampled_", sGenus2, ".txt", sep="") , header=T);

fnOverlap <- function(arr1, arr2) {
  arrLarger <- arr1;
  arrSmaller <- arr2;
  if (mean(arr1) < mean(arr2) ) {
    arrLarger <- arr2;
    arrSmaller <- arr1;
  }
  nOverlap1 <- length(arrLarger[arrLarger<= max(arrSmaller)]);
  nOverlap2 <- length(arrSmaller[arrSmaller>= min(arrLarger)]);
  #cat(nOverlap1, nOverlap2)
  return( (nOverlap1+nOverlap2) / (length(arr1) + length(arr2) ) );
}

arrP <- c();
for(sCol in colnames(datResampled1)) {
  arrP <- c(arrP , fnOverlap(datResampled1[, sCol], datResampled2[, sCol]) );
}

names(arrP) <- colnames(datResampled1);
write.table(t(arrP) , file=paste("compare_", sGenus1, "_", sGenus2, ".txt", sep=""), col.names = T, row.names = F, quote = F, sep="\t")
