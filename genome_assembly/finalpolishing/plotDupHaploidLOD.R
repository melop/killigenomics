dat <- read.table("HaploidLODdistr.txt" , header=F, na.strings=c("NAN" , "INF", "-INF") );
hist(dat$V1,breaks=50, xlab="Haploid LOD score of duplicated segments");
