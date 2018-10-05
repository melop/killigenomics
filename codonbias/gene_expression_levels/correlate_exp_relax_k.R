setwd("~/killifish_genomes/codonbias/gene_expression_levels/")
sSp <- "EPI";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Epiplatys.txt";
sSp <- "ARC";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Archiaphyosemion.txt";
sSp <- "SCR";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Scriptaphyosemion.txt";
sSp <- "FUN";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Fundulopanchax.txt";
sSp <- "AKN";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_AKN.txt";
sSp <- "CTO";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Callopanchax.txt";
sSp <- "AAU";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Aphyosemion.txt";
sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Aphyosemion.txt";

#sSp <- "NOR";
#sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/sum_Nothobranchius.txt";
#sKFile <- "~/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0/sum_Nothobranchius.txt";


datExp <- read.table("all.expression.level.txt", header= T, sep="\t", quote = '');
datK <- read.table(sKFile, header=F, sep="\t", quote='');

sink(file = paste("correlate_relaxK_rpkm_", sSp, ".txt", sep=''));
pdf(file = paste("correlate_relaxK_rpkm_", sSp, ".pdf", sep=''), width=5, height=5);
datK <- datK[datK$V5 <=1 & datK$V9 > (1/30) &  datK$V9 < 30 , ];
datExp <- datExp[datExp$SdOverAvg < 0.5, ];
datMerge <- merge(datK[, c(1, 9)] ,  datExp, by.x='V1', by.y='OrthoID');
colnames(datMerge)[1] <- 'OrthoID';
colnames(datMerge)[2] <- 'relaxk';

oModel <- summary(lm( datMerge$relaxk ~ datMerge$AvgRPKM  ))
oModel

plot(  datMerge$AvgRPKM, datMerge$relaxk  , pch=16, col=rgb(0,0,1,1/16), ylab=" K ", xlab=" Average RPKM ");
abline(a = oModel$coefficients[1,1], b= oModel$coefficients[2,1] , col='red');


datMerge$AvgRPKM[datMerge$AvgRPKM ==0] <- 0.01; # do not allow 0 , set it to a small number to allow log transformation

datMerge$LogAvgRPKM <- log10(datMerge$AvgRPKM);
datMerge$LogRelaxK <- log10(datMerge$relaxk);

oModel <- summary(lm(datMerge$LogRelaxK ~ datMerge$LogAvgRPKM  ))
oModel

plot(  datMerge$LogAvgRPKM, datMerge$LogRelaxK  , pch=16, col=rgb(0,0,1,1/16), ylab="Log10(K)", xlab="log10(Average RPKM)");
abline(a = oModel$coefficients[1,1], b= oModel$coefficients[2,1] , col='red');

sink();
dev.off();
