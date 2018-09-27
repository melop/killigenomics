#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/finalrun/");
dat <- read.table("demultiplexed/reports.txt", header=T);
datInc <- read.table("includesamples.txt", header=T);
dat2 <- merge(dat,datInc, by= "SampleID")
dat3 <- dat2[dat2$included_in_pool == "yes", ];

dat3$estcov <- dat3$TotalBases / 1548592001;

sum(as.numeric( dat3$TotalBases) );

boxplot(dat3$estcov ~ dat3$Species)

datAgg <- aggregate(dat3[, c("estcov") ] , by = list(dat3$Species), mean);
datAgg$NormX <- datAgg$x / min(datAgg$x)
datAgg

plot(dat3$pool_vol,dat3$TotalBases )



require(ggplot2)    
p <- ggplot(data = dat3, aes(x=estcov)) 
p <- p + geom_histogram(aes(weights=estcov, fill=Species))
p <- p + scale_fill_brewer(palette="Set3")
p <- p + facet_wrap( ~ Species, ncol=1)
p

quit();

dat3$MappedCov <- 0;
for(i in 1:nrow(dat3) ) {
  dat3[i, "MappedCov"] <- fnGetRealCov(dat3[i,]);
}

plot(dat3$estcov, dat3$MappedCov)

#GENOME SIZE ESTIMATES:
mean(dat3[dat3$Species=="ORTWET", "estcov"] / dat3[dat3$Species=="ORTWET", "MappedCov"] ) * 1.1
mean(dat3[dat3$Species=="ORTDRY", "estcov"] / dat3[dat3$Species=="ORTDRY", "MappedCov"] ) * 1.1
mean(dat3[dat3$Species=="RACWET", "estcov"] / dat3[dat3$Species=="RACWET", "MappedCov"] ) * 1.1
mean(dat3[dat3$Species=="RACDRY", "estcov"] / dat3[dat3$Species=="RACDRY", "MappedCov"] ) * 1.1

#[1] 1.424889
#[1] 1.428899
#[1] 1.318414
#[1] 1.320455



require(ggplot2)    
p <- ggplot(data = dat3, aes(x=MappedCov)) 
p <- p + geom_histogram(aes(weights=MappedCov, fill=Species))
p <- p + scale_fill_brewer(palette="Set3")
p <- p + facet_wrap( ~ Species, ncol=1)
p

datAgg <- aggregate(dat3[, c("MappedCov") ] , by = list(dat3$Species), mean);
datAgg$NormX <- datAgg$x / min(datAgg$x)
datAgg

datGC <- read.table("testmapping/ref.fa.gc_out.txt", header=T);

nMeanN <- mean( datGC[datGC$Total_Count > 50000, "N_Count"] / datGC[datGC$Total_Count > 50000, "Total_Count"])

dat3$MappedCovCorrN <- dat3$MappedCov / (1-nMeanN); #adjust for the "N"s, these do not receive mapping.

datAgg <- aggregate(dat3[, c("MappedCovCorrN") ] , by = list(dat3$Species), mean);
datAgg$NormX <- datAgg$x / min(datAgg$x)
datAgg

p <- ggplot(data = dat3, aes(x=MappedCovCorrN)) 
p <- p + geom_histogram(aes(weights=MappedCov, fill=Species))
p <- p + scale_fill_brewer(palette="Set3")
p <- p + facet_wrap( ~ Species, ncol=1)
p

fnGetRealCov <- function(datRow) {
  sFileName <- paste("testmapping/bams/", datRow$Species, "_", datRow$Pop, "_", datRow$SampleID, "_", datRow$i5, "_" , datRow$i7, ".cov.txt", sep="")
  cat("doing ", sFileName,"\n");
  if (!file.exists(sFileName)) {
    return(0);
  }
  
  datCov <- read.table(sFileName, header=F);
  colnames(datCov) <- c("scf" , "cov_bin", "sitecount" , "len", "freq");
  datCov$depth <- datCov$cov_bin * datCov$sitecount;
  datCovAgg <- aggregate(datCov[, c("depth")], list(datCov$scf), sum );
  colnames(datCovAgg) <- c("scf", "sumdepth");
  
  datCovMerge <- unique( merge(datCovAgg, datCov[, c('scf', 'len') ], by='scf'));
  datCovMerge$cov <-  datCovMerge$sumdepth / datCovMerge$len;
  datCovMerge <- datCovMerge[datCovMerge$cov < 10 & datCovMerge$len < 0.5e9, ]
  #plot(datCovMerge$len , datCovMerge$cov , cex=0.5)
  #boxplot(datCovMerge$cov)
  return(mean(datCovMerge[datCovMerge$len > 50000, "cov"] ))
}
