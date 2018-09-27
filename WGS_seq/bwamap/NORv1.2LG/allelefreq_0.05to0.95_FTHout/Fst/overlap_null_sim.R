# nWindows1 <- 96580;
# nSig1 <- 800;
# 
# nWindows2 <- 96104;
# nSig2 <- 776;

nWindows1 <- 96580/5.625;
nSig1 <- 800/5.625;

nWindows2 <- 96104/5.625;
nSig2 <- 776/5.625;
arrOverlapSim <- c();

arr1 <- 1:nWindows1;
arr2 <- 1:nWindows2;

for(i in 1:10000) {
  arrSample1 <- sample(x = arr1, size = nSig1);
  arrSample2 <- sample(x = arr2, size = nSig2);
  arrSample2 <- c(arrSample2 , arrSample2+1,arrSample2-1);
  #arrSample1 <- c(arrSample1 );
  #arrSample2 <- c(arrSample2 , arrSample2+1,  arrSample2+2,arrSample2+3,arrSample2+4,arrSample2+5, arrSample2-1,  arrSample2-2, arrSample2-3,  arrSample2-4, arrSample2-5);
  arrOverlapSim <- c(arrOverlapSim, length(intersect(arrSample1, arrSample2)));
  
}