setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq/Fst/pergene/");
sPop <- "RAC";
sPop2<- "ORT"
nMinSites <- 5;
nPcutoff <-  0.1; ##find genes where in both species, the FSTs are abnormal


peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=200)
  ins <- dh[["density"]]
  nbins <- length(ins)
  ss <- which(rank(ins) %in% seq(from=nbins-30,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

sOut <- paste("zcat -f ", sPop, ".pergene.", ".minsite", nMinSites , ".txt.gz" , sep="");
pOut <- pipe(sOut)

datPerGeneFST <- read.table(pOut , header=T);

m <- mean(peakfinder(datPerGeneFST$fst)); #arr50To84[1]; #mean(datPerGeneFST$fst[datPerGeneFST$fst >=arr50To84[1] &  datPerGeneFST$fst <=arr50To84[2]]);
arrTmpfst <- datPerGeneFST$fst[datPerGeneFST$fst>=m];
arrQuantiles <- quantile(arrTmpfst , c(0.3413*2), na.rm=T)
std <- arrQuantiles[1] - m;
datPerGeneFST$p <- pnorm(datPerGeneFST$fst, mean=m, sd=std);
datPerGeneFST$p[datPerGeneFST$fst <= m] <- datPerGeneFST$p[datPerGeneFST$fst <= m] * 2; #datPerGeneFST$p[datPerGeneFST$p <= 0.5] * 2; # ignore lowly divergent regions;
datPerGeneFST$p[datPerGeneFST$fst > m] <- (1-datPerGeneFST$p[datPerGeneFST$fst > m] ) * 2;
datPerGeneFST<- datPerGeneFST[datPerGeneFST$p < nPcutoff , ];
datPerGeneFST$fstHigherThanAvg <- (datPerGeneFST$fst > m);


sOut2 <- paste("zcat -f ", sPop2, ".pergene.", ".minsite", nMinSites , ".txt.gz" , sep="");
pOut2 <- pipe(sOut2)

datPerGeneFST2 <- read.table(pOut2 , header=T);

m <- mean(peakfinder(datPerGeneFST2$fst)); #arr50To84[1]; #mean(datPerGeneFST$fst[datPerGeneFST$fst >=arr50To84[1] &  datPerGeneFST$fst <=arr50To84[2]]);
arrTmpfst <- datPerGeneFST2$fst[datPerGeneFST2$fst>=m];
arrQuantiles <- quantile(arrTmpfst , c(0.3413*2), na.rm=T)
std <- arrQuantiles[1] - m;
datPerGeneFST2$p <- pnorm(datPerGeneFST2$fst, mean=m, sd=std);
datPerGeneFST2$p[datPerGeneFST2$fst <= m] <- datPerGeneFST2$p[datPerGeneFST2$fst <= m] * 2; #datPerGeneFST$p[datPerGeneFST$p <= 0.5] * 2; # ignore lowly divergent regions;
datPerGeneFST2$p[datPerGeneFST2$fst > m] <- (1-datPerGeneFST2$p[datPerGeneFST2$fst > m] ) * 2;

datPerGeneFST2<- datPerGeneFST2[datPerGeneFST2$p < nPcutoff , ];
datPerGeneFST2$fstHigherThanAvg <- (datPerGeneFST2$fst > m);


datMerge <- merge(datPerGeneFST , datPerGeneFST2[,c(1,6,7,8)], by="OrthoId");

summary(lm(datMerge$fst.x ~ datMerge$fst.y));

plot(datMerge$fst.y , datMerge$fst.x);

sRelax <- "/beegfs/group_dv/home/RCui/killifish_genomes/figures/Fig4/NothoRelax/sum.txt";
datRelax <- read.table(sRelax , header=F, sep="\t", quote="");
datRelax <- datRelax[complete.cases(datRelax),];
colnames(datRelax)[c(3,23,5,9)] <- c("genesymbol", "genename", "relax_p", "relax_k");
datMerge <- merge(datMerge, datRelax[,c('V1', 'genesymbol', 'genename', 'relax_p', 'relax_k')], by.x="OrthoId", by.y ="V1")


sCodeml <- "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/sum_Nothobranchius.txt"
datCodeml <- read.table(sCodeml , header=F, sep="\t", quote="", fill=T);
datCodeml <- datCodeml[complete.cases(datCodeml),];
colnames(datCodeml)[c(5)] <- c( "codeml_p");
datMerge <- merge(datMerge, datCodeml[,c('V1', 'codeml_p')], by.x="OrthoId", by.y ="V1")

datCoolGenes <- datMerge[datMerge$fstHigherThanAvg.x == datMerge$fstHigherThanAvg.y,];
datCoolGenes$relax_k[datCoolGenes$relax_k ==0] <- min(datCoolGenes$relax_k[datCoolGenes$relax_k > 0])
View(datCoolGenes[datCoolGenes$relax_p <= 0.05 & datCoolGenes$codeml_p <=0.05,]);
#boxplot( log10(datCoolGenes$relax_k) ~ datCoolGenes$fstHigherThanAvg.x);

#t.test( log10(datCoolGenes[datCoolGenes$fstHigherThanAvg.x, 'relax_k']) , log10( datCoolGenes[!datCoolGenes$fstHigherThanAvg.x, 'relax_k']) , paired = F)

