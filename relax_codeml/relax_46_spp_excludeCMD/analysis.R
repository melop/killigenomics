setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/relax_46_spp/");

nPerGeneMaxDiff <- 50;
nMaxMedianDiff <- 50;

datSum <- read.table("sum_Scriptaphyosemion.txt", header=F, fill=T, stringsAsFactors = F, sep = "\t", quote='');
datSum <- datSum[datSum$V2 == 'Success', ];

datSum$V4 <- as.numeric(datSum$V4)
datSum$V5 <- as.numeric(datSum$V5)
datSum$V9 <- as.numeric(datSum$V9)

datSum$FDR <- p.adjust(datSum$V5, method="fdr")

datIntense <- datSum[datSum$V9>1, ] 
datRelax<- datSum[datSum$V9<1, ] 
hist(as.numeric(datRelax$V5), breaks=50 , col=rgb(0,0,1,1/4), xlab="Relax p value", main="", ylim=c(0,800))
hist(as.numeric(datIntense$V5), breaks=50, col=rgb(1,0,0,1/4) , add=T)

plot(datRelax$V4, datRelax$V5)
summary(lm(datRelax$V5 ~ datRelax$V4))

datRelaxSig <- datRelax[datRelax$FDR < 0.1, ]
View(datRelaxSig);

datIntenseSig <- datIntense[datIntense$FDR < 0.1, ]
View(datIntenseSig);


arrTestSet <- as.numeric(datRelaxSig$V4);
names(arrTestSet) <- datRelaxSig$V1;
arrTestSet <- arrTestSet[!is.na(arrTestSet)];

arrBGSet <- as.numeric(datSum$V4);
names(arrBGSet) <- datSum$V1;
arrBGSet <- arrBGSet[!is.na(arrBGSet)];


hist(log10(arrBGSet) , breaks=30, col=rgb(1,0,0,1/4));
hist(log10(arrTestSet) , breaks=30, col=rgb(0,0,1,1/4),  add=T);

#hist( datSum$V4 , breaks=30, col=rgb(1,0,0,1/4));
#hist(datRelaxSig$V4 , breaks=30, col=rgb(0,0,1,1/4),  add=T);



fnMatchGOBackgroundLength <- function(arrTestSet, arrBGSet) {
  arrLogTestSet <- sort(log10(arrTestSet) );
  arrDrawnBg <- sort(arrLogTestSet); #first, put all test sets into the background.
  arrRemainSet <- log10(sort(arrBGSet[!(names(arrBGSet) %in% names(arrLogTestSet) ) ]));
  nCurrExceedingCount <- 0;
  
  
  set.seed(100);
  while( abs(median(10^arrLogTestSet) - median(10^arrDrawnBg)) <= nMaxMedianDiff ) {
    nRandLen <- sample(arrLogTestSet, size = 1);#rnorm(n = 1, mean=nTestMean, sd=nTestSd);
    arrDiff <- sort( abs(arrRemainSet - nRandLen) );
    sGene <- names(arrDiff[1]);
    if ( 10^arrDiff[1] > nPerGeneMaxDiff ) next;
    arrDrawnBg <- c(arrDrawnBg, arrRemainSet[sGene]);
    arrRemainSet <- arrRemainSet[names(arrRemainSet) != sGene];
    if (length(arrRemainSet) == 0 ) break;
  }
  
  hist(arrDrawnBg , breaks=30, col=rgb(1,0,0,1/4));
  hist(arrLogTestSet , breaks=30, col=rgb(0,0,1,1/4),  add=T);
  
  return(arrDrawnBg)
}

arrDrawnBg <- fnMatchGOBackgroundLength(arrTestSet, arrBGSet);


