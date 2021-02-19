setwd("/beegfs/group_dv/home/RCui/killifish_genomes/plp_treemix");
sF <- "cov_ref_PML"
sF <- "cov_ref_A.PraslinPasquereRiver"
sF <- "cov_ref_B.PraslinPlaineHollandaise"
sF <- "cov_ref_C.PraslinWest"
sF <- "cov_ref_D.LaDigue"
sF <- "cov_ref_D.LaDigue"
sF <- "cov_ref_E.MaheNorth"
sF <- "cov_ref_F.Curieuse"
sF <- "cov_ref_G.MaheSouth"
sF <- "cov_PLP.PDC.2.txt"



dat <- read.table(sF, header=T)

plot(dat$Depth, dat$Count, xlim=c(1,200));
nPrevCount <- nPeakCount <- max(dat$Count[dat$Depth>0]);

nPeakDepth <- dat$Depth[[which(dat$Count ==nPeakCount)]]


#search downhill
nLowerBound <- 0;
for(nDP in nPeakDepth:1) {
  nCurrCount <- dat[dat$Depth==nDP, 'Count' ];
  if (nCurrCount <= nPrevCount) {
    nPrevCount <- nCurrCount;
  } else {
    nLowerBound <- nDP+1;
    break;
  }
}


#search downhill
# nPrevCount <- nPeakCount <- max(dat$Count);
# nUpperBound <- 0;
# for(nDP in nPeakDepth:max(dat$Depth)) {
#   nCurrCount <- dat[dat$Depth==nDP, 'Count' ];
#   if (nCurrCount <= nPrevCount) {
#     nPrevCount <- nCurrCount;
#   } else {
#     nUpperBound <- nDP+1;
#     break;
#   }
# }

nUpperBound <- nPeakDepth + abs(nPeakDepth-nLowerBound);
if (nUpperBound < 50) {
  nUpperBound <- 50
}

datGood <- dat[dat$Depth>=nLowerBound & dat$Depth<=nUpperBound, ];
nSd <-  sqrt( sum( (datGood$Depth - nPeakDepth)^2 * datGood$Count ) / (sum(datGood$Count)-1) )

datGood$Freq <- datGood$Count / sum(datGood$Count)
plot(datGood$Depth, datGood$Freq, xlim=c(1,200));
lines(datGood$Depth, dnorm(datGood$Depth,nPeakDepth,nSd));

arrSim <- rnorm(10000, mean = nPeakDepth, sd = nSd);

arrQuantiles <- quantile(arrSim, c(0.025, 0.975) )

sFout <- paste(sF, ".cov.cutoffs.txt", sep="");
write.table(data.frame(pop=sF, peak=nPeakDepth, low=arrQuantiles[1], high=arrQuantiles[2]), file = sFout, col.names = F, row.names = F, sep="\t" , quote=F);
