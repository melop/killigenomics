library(qqman)
setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/Fst");

#sPop <- "RAC";
sPop <- "ORT";
nPcutoff <-  0.01;


nWindowSize <- 10000;
nStepSize <- 5000;
nMinSites <- 10;

pIn <- pipe(paste("zcat -f ", sPop, ".slideFST.win",nWindowSize, ".step", nStepSize, ".minsite", nMinSites , ".txt.gz" , sep=""));
datFST <- read.table(pIn, header=T, sep="\t", fill = T)
datFST <- datFST[grepl("chr", datFST$chr), ];#ignore scaffolds.
datFST$ChrNum <- as.numeric(gsub('[^0-9]*([0-9]+)','\\1', as.character( datFST$chr) ));

#pdf(file=paste( sPop, ".slideFST.win",nWindowSize, ".step", nStepSize, ".minsite", nMinSites , ".pdf", sep="") , width=10,height=4);
peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=200)
  ins <- dh[["density"]]
  nbins <- length(ins)
  ss <- which(rank(ins) %in% seq(from=nbins-30,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

m <- mean(peakfinder(datFST$fst)); #arr50To84[1]; #mean(datFST$fst[datFST$fst >=arr50To84[1] &  datFST$fst <=arr50To84[2]]);
arrTmpFst <- datFST$fst[datFST$fst>=m];
hist(datFST$fst, freq = F, prob=TRUE, density=0, breaks=100);
arrQuantiles <- quantile(arrTmpFst , c(0.3413*2), na.rm=T)
std <- arrQuantiles[1] - m;

curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

datFST$p <- pnorm(datFST$fst, mean=m, sd=std);
datFST$p[datFST$fst <= m] <- 1;#datFST$p[datFST$p <= 0.5] * 2; # ignore lowly divergent regions;
datFST$p[datFST$fst > m] <- (1-datFST$p[datFST$fst > m] ) * 2;

manhattan(datFST, chr="ChrNum", bp = 'midpoint', p = "p", logp = T, genomewideline = 4)

datOutBed <- datFST[datFST$p < nPcutoff,];

datOutBed$Info <- paste("fsthigherthanavg:", datOutBed$fst > m , "fst:",datOutBed$fst, ";","p:",datOutBed$p, sep="");

datOutBed$left <- datOutBed$left - 1;
write.table(datOutBed[, c('chr', 'left', 'right', 'Info')], file=paste( sPop, ".slideFST.win",nWindowSize, ".step", nStepSize, ".minsite", nMinSites , ".sig",nPcutoff,".bed", sep=""), row.names = F, col.names = F, quote=F, sep="\t");
#dev.off()
nrow(datFST);
nrow(datOutBed);

