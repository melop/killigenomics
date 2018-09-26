setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/rerun_exclude_MNM");

sRet1 <- "sum_Nothobranchius.txt";
sRet2 <- "../sum_Nothobranchius.txt";

dat1 <- read.table(sRet1, header=F, sep="\t", quote="", fill = T);
dat2 <- read.table(sRet2, header=F, sep="\t", quote="", fill = T);
datMerge <- merge(dat1, dat2, by=1)

datBothSig <- datMerge[datMerge$V5.x < 0.05 & datMerge$V5.y < 0.05, ];
datOnly1Sig <- datMerge[datMerge$V5.x < 0.05 & datMerge$V5.y >= 0.05, ];
datOnly2Sig <- datMerge[datMerge$V5.x >= 0.05 & datMerge$V5.y < 0.05, ];

arrMat <- matrix( c(nrow(datBothSig), nrow(datOnly1Sig), nrow(datBothSig), nrow(datOnly2Sig)), byrow = T, nrow = 2 )
colnames(arrMat) <- c("Overlap", "NonOverlap");
rownames(arrMat) <- c("Study1", "Study2");
arrMat
fisher.test(arrMat)
