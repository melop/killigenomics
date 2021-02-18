#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/size_est/mappedtoBUSCO");

sOutFile <-"cov.all.txt";
arrSkip <- c();
datRet <- NULL;

if (file.exists(sOutFile)) {
  datRet <- read.table(sOutFile, header=T, stringsAsFactors = F);
  arrSkip <- datRet$species;
}


arrCountFolders <- Sys.glob("basecount/*");

for( sCountFolder in arrCountFolders ) {
  sSp <- basename(sCountFolder);
  cat("do ", sSp, " ");
  if (sSp %in% arrSkip) {
    cat("Skip", sSp, "...\n");
    next;
  }
  
  if (length(Sys.glob(paste( "mapped/ref*/",sSp,".cov.ok", sep="") )) !=1 ) {
    cat("bedtools coverage output error ", sSp, "\n");
    next;
    }
  arrBedCovFile <- Sys.glob(paste( "mapped/ref*/",sSp,".cov.txt", sep="") );
  if (length(arrBedCovFile) !=1) next;
  sBedCovFile <- arrBedCovFile[1];
  cat("reading ", sBedCovFile, "\n");
  datCov <- read.table(sBedCovFile  , header=F);
  if (nrow(datCov)==0) {
    cat("bedtools coverage output empty ", sSp, "\n");
    next;
  }
  nAvgCov <- mean(datCov$V5)
  
  datBaseCount <- read.table( paste(sCountFolder, "/counts.txt" ,sep=""), header=T );
  nTotalBase <- datBaseCount$Both_bp[1];
  datRet <- rbind(datRet, list(species=sSp, cov=nAvgCov, readbase=nTotalBase, genomesize=nTotalBase/nAvgCov) );
}

write.table(datRet, file=sOutFile, quote=F, row.names = F, col.names = T, sep="\t" );
