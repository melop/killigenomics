setwd("/beegfs/group_dv/home/RCui/killifish_genomes/codonbias/46spp/codon_usage");
sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0/sum_Nothobranchius.txt";
#sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Callopanchax.txt";
#sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Aphyosemion.txt";
#sK <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Scriptaphyosemion.txt";

sCAI <- "CAI_min_100AA.txt";
sSp <- 'NOR';
sOutGroup <- 'AAU';

datK <- read.table(sK, header=F, sep="\t", quote='');
datCAI <- read.table(sCAI, header=T);

datK <- datK[datK$V5 <= 1 & datK$V9 >= 1/30 & datK$V9 <= 30, c('V1', 'V9') ];
datCAI <- datCAI[, c('GeneName', paste(sSp, '_CAI', sep='') , paste(sOutGroup, '_CAI', sep='')) ];

datMerge <- merge(datK, datCAI, by.x='V1' , by.y='GeneName');
colnames(datMerge) <- c('OrthoID', "K", "CAI", "Outgroup");
# 
# datMerge$CAIContrast <- datMerge$CAI - datMerge$Outgroup
# summary(lm( datMerge$CAI ~ datMerge$K ))
# summary(lm( datMerge$CAIContrast ~ datMerge$K ))
# 
# 
# plot( log10(datMerge$K) , log10(datMerge$CAI), pch=16, cex=0.5, col=rgb(0,0,1,1/8));
# plot(datMerge$CAI, datMerge$K , pch=16, cex=0.5, col=rgb(0,0,1,1/8));

oModel <- summary(lm(  log10(datMerge$K) ~ datMerge$CAI ))

oModel

plot(  datMerge$CAI  , log10(datMerge$K), pch=16, col=rgb(0,0,1,1/16), ylab="Log10(K)", xlab="Codon Adaptation Index (CAI)");
abline(a = oModel$coefficients[1,1], b= oModel$coefficients[2,1] , col='red', lwd=4);

