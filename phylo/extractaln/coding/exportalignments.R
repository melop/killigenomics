#Change the following inputs

#sFile <- "ret_codon12_table.txt";
#sFile <- "ret_4fold_table.txt";
sFile <- "ret_allcodons_table.txt";

sFastaOut <- paste(sFile,".fasta",sep="");

dat <- read.table(sFile, header=T, as.is=T);
length(unique(dat$GeneName));

arrTaxa <- colnames(dat);
arrTaxa <- arrTaxa[4:length(arrTaxa)];
cat("" , file=sFastaOut, append = F);

for(sTaxon in arrTaxa) {
  
  cat(">" , sTaxon, "\n", sep="", gsub('(.{1,100})', "\\1\n", paste(as.character(dat[,sTaxon]), collapse="") ), file=sFastaOut, append = T);
  
}
