#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref/enrichment");
datTarget <- read.table("plp.nfzref.upstreamcds.gff", header=F, sep="\t");
datBg<- read.table("annualconserved.nfzref.upstreamcds.gff", header=F, sep="\t");

pdf(file="dist.pdf", width=5, height=5);
hBg <- hist( log10(abs(datBg$V19)), col = rgb(0,0,1,1/3), add=F, freq=F, breaks=50, )
hTg <- hist( log10(abs(datTarget$V19)), col = rgb(1,0,0,1/3) , add=T, freq=F, breaks=50)
dev.off();
cat("median background: ", median(datBg$V19));
cat("median plp: ", median(datTarget$V19));

wilcox.test(datBg$V19 , datTarget$V19);

 datM <- merge(data.frame(Mid=round(as.numeric(hBg$mids) , 4), BgDensity=hBg$density), data.frame(Mid= round(as.numeric(hTg$mids) , 4), TgDensity=hTg$density), by="Mid", all = F)
datM$diff <- datM$BgDensity - datM$TgDensity

for( i in length(datM$diff):1 ) {
  if (datM$diff[i]>0) {
    
    cat("Chosen distance cutoff is: 10^", datM$Mid[i], "  = " , 10^datM$Mid[i], "\n");
    break;
  } 
}