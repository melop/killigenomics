#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_0.01to0.98_NVGout/anavar_calibrated_err0.0033/runanavar_bs");

arrInFiles <- Sys.glob("../*.SFS.for.anavar.txt")
nReps <- 100;

fnResample <- function(arr) {
  nSum <- sum(arr);
  arrProb <- arr / nSum;
  return(as.numeric(rmultinom(1, size = nSum, prob = arrProb )));
}

for(sFile in arrInFiles) {
    sBaseName <- basename(sFile)
    arrF <- unlist(strsplit(sBaseName,'\\.'));
    sPop <- arrF[1];
    sOutDir <- paste('bs/', sPop, sep="");
    dir.create(sOutDir ,F, T );
    arrLines <- readLines(sFile);
    for(nRep in 1:nReps) {
      hOut <- file( paste(sOutDir, '/rep', nRep, '.txt', sep='') , "w");
      arrOut <- c();
      for(nLn in 1:length(arrLines) ) {
        sLn <- arrLines[nLn]
        arrF <- unlist(strsplit(sLn,':'));
        if (arrF[1] == 'neu_sfs' || arrF[1] == 'sel_sfs') { #resample this line
          arrCounts <- as.numeric(unlist(strsplit(arrF[2],'\\,')));
          arrOut <- c(arrOut, paste(arrF[1], ': ', paste( fnResample(arrCounts) , collapse = ', ') , sep="" ));
        } else {
          arrOut <- c(arrOut, sLn);
        }
      }
      writeLines(text = arrOut , con = hOut)
      close(hOut);
    }
}

