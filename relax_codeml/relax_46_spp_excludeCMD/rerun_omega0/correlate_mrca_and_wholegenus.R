setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0");

s1 <- "sum_Nothobranchius.txt";
s2 <- "sum_Nothobranchius_mrca.txt";

dat1 <- read.table(s1, header=F, sep="\t" , quote="");
dat2 <- read.table(s2, header=F, sep="\t" , quote="");

datMerge <- merge(dat1[, c(1,5,9)], dat2[, c(1,5,9)], by=1)

summary(lm(datMerge$V9.x ~ datMerge$V9.y));
summary(lm(datMerge$V5.x ~ datMerge$V5.y));

datMerge$logK1 <- log10(datMerge$V9.x); 
datMerge$logK2 <- log10(datMerge$V9.y); 
