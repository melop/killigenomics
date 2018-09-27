setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/Fst/pergene/");
sPop <- "RAC";
nMinSites <- 5;
nPcutoff <-  1e-3;

datGeneNames <- read.table("~/killifish_genomes/annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt" , header=T, sep="\t" , quote="");
datGeneNames <- datGeneNames[, c(1,2,3)];

sOut <- paste("zcat -f ", sPop, ".pergene.", ".minsite", nMinSites , ".txt.gz" , sep="");
pOut <- pipe(sOut)

datPerGeneFST <- read.table(pOut , header=T);

peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=200)
  ins <- dh[["density"]]
  nbins <- length(ins)
  ss <- which(rank(ins) %in% seq(from=nbins-30,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

m <- mean(peakfinder(datPerGeneFST$fst)); #arr50To84[1]; #mean(datPerGeneFST$fst[datPerGeneFST$fst >=arr50To84[1] &  datPerGeneFST$fst <=arr50To84[2]]);
#m <- mean(datPerGeneFST$fst);
#std <- sd(datPerGeneFST$fst);
arrTmpfst <- datPerGeneFST$fst[datPerGeneFST$fst>=m];

pdf(paste( sPop, ".pergene.", ".minsite", nMinSites , ".pdf" , sep=""), width=5,height=5);
hist(datPerGeneFST$fst, freq = F, prob=TRUE, density=0, breaks=100);
arrQuantiles <- quantile(arrTmpfst , c(0.3413*2), na.rm=T)
std <- arrQuantiles[1] - m;

curve(dnorm(x, mean=m, sd=std),       col="darkblue", lwd=2, add=TRUE, yaxt="n")


datPerGeneFST$p <- pnorm(datPerGeneFST$fst, mean=m, sd=std);
datPerGeneFST$p[datPerGeneFST$fst <= m] <- datPerGeneFST$p[datPerGeneFST$fst <= m] * 2; #datPerGeneFST$p[datPerGeneFST$p <= 0.5] * 2; # ignore lowly divergent regions;
datPerGeneFST$p[datPerGeneFST$fst > m] <- (1-datPerGeneFST$p[datPerGeneFST$fst > m] ) * 2;

#manhattan(datPerGeneFST, chr="ChrNum", bp = 'midpoint', p = "p", logp = T, genomewideline = 3, ylim=c(0, 8))

datOutBed <- datPerGeneFST[datPerGeneFST$p < nPcutoff,];

#datOutBed$Info <- paste("fst:",datOutBed$fst, ";","p:",datOutBed$p, sep="");

datOutBed <- merge(datOutBed, datGeneNames, by.x='OrthoId', by.y='Group_Id');

datOutBed$left <- datOutBed$left - 1;
write.table(datOutBed[, c('chr', 'left', 'right', 'fst', 'p', 'OrthoId', 'Gene_Symbol', "Gene_Name")], file=paste( sPop, ".minsite", nMinSites , ".sig",nPcutoff,".bed", sep=""), row.names = F, col.names = T, quote=F, sep="\t");
nrow(datPerGeneFST);
nrow(datOutBed);
dev.off();
      