arrPops <- c( "ORTWET.pergene.alpha.txt", "RACWET.excludehybrid.pergene.alpha.txt", "RACDRY.excludehybrid.pergene.alpha.txt", "ORTDRY.excludehybrid.pergene.alpha.txt");
sOut <- "CombinedPops.pergene.alpha.txt";
nPopCount <- length(arrPops);
#Get the list of genes that agree in the sign of alpha (>0 or <0) in all populations, and output the lowest p and a multiplied p in all pops:
datCombined <- NULL;
for(sPop in arrPops) {
  datMK <- read.table(sPop, header=T, quote='', sep="\t");
  if (is.null(datCombined) ) {
    datCombined <- datMK[ , c(1,2,3,4,5,6,7,8,9,14)];
    next;
  }
  datCombined <- merge(x=datCombined, y=datMK[, c(1,2,3,4,5,6,7,8,9,14)], by=1:7, all=T);
}

arrAlphaCols <- seq(8,nPopCount*3+7, 3 );
arrPCols <- arrAlphaCols + 1;
arrDoSCols <- arrAlphaCols + 2;

for(nRow in 1:nrow(datCombined)) {
  bAllAgree <- T;
  bPrevNeg <- NULL;
  arrAlphas <- c();
  arrPs <- c();
  arrDoSs <- c();
  for(nPop in 1:nPopCount) {

    if (is.na(datCombined[nRow , arrPCols[nPop] ] ) ) {
      next;
    }

    #check if the p value is =1
    #if (datCombined[nRow , arrPCols[nPop] ] == 1) {#
	#next;
    #}

    if (is.null(bPrevNeg)) {
      bPrevNeg <- (datCombined[nRow , arrDoSCols[nPop]] <0);
      arrAlphas <- datCombined[nRow , arrAlphaCols[nPop]];
      arrPs <- datCombined[nRow , arrPCols[nPop] ];
      arrDoSs <- datCombined[nRow , arrDoSCols[nPop]];
      
      next;
    }
    
    if (bPrevNeg != (datCombined[nRow , arrDoSCols[nPop]] <0) ) {
      bAllAgree <- F;
      break;
    }
    
    arrAlphas <- c( arrAlphas, datCombined[nRow , arrAlphaCols[nPop]]);
    arrPs <- c(arrPs, datCombined[nRow , arrPCols[nPop] ]);
    arrDoSs <- c( arrDoSs, datCombined[nRow , arrDoSCols[nPop]]);
    
    
  }
  
  if (bAllAgree & length(arrAlphas)>=1 ) {
    datCombined[nRow , 'avgalpha'] <- mean(arrAlphas);
    datCombined[nRow , 'lowestp'] <- min(arrPs);
    datCombined[nRow , 'multipliedp'] <- prod(arrPs);
    datCombined[nRow , 'avgDoS'] <- mean(arrDoSs);
    
  }
}

datCombined <- datCombined[!is.na(datCombined$avgalpha), ];
write.table(datCombined, file=sOut, quote=F, row.names = F, col.names=T, sep="\t");
