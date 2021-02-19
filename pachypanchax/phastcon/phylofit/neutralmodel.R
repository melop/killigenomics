setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref");

library(rphast)

#args = commandArgs(trailingOnly=TRUE)

nChr <- 1; #use chromosome 1 to infer the neutral model

sOutStem <- paste("neutralmodel.chr",nChr, sep="");


sRefSp <- 'NFZ'
tree <- "((austrofundulus,kryptolebias),(PLP,(CTO,(AAU,(NOR,NFZ)))))";

arrAnnuals <- c( 'NOR', 'CTO', 'austrofundulus' );
arrNonAnnuals <- c('AAU', 'PLP', 'kryptolebias');

treeAnnuals <- "(austrofundulus,(CTO,(NOR,NFZ)))";

feats <- read.feat(paste("chr",nChr,"_cds.gff", sep=""))
alnAnnual <- read.msa(filename = paste("../../cactus/split_NFZ_maf/NFZ_ref_5spp_chr",nChr,".maf", sep=""), format = "MAF", ordered = T, seqnames = c(sRefSp, arrAnnuals) )

names(alnAnnual)
feats$seqname <- sRefSp;
table(feats$feature)

aln4d <- get4d.msa(alnAnnual, feats)
neutralMod <- phyloFit(aln4d, tree=tree, subst.mod="REV")
write.tm(neutralMod , paste(sOutStem, ".annualspp.mod", sep=""));
alnAll <- read.msa(filename = paste("../../cactus/split_NFZ_maf/NFZ_ref_5spp_chr",nChr,".maf", sep=""), format = "MAF", ordered = T )
aln4dAll <- get4d.msa(alnAll, feats)
neutralModAll <- phyloFit(aln4dAll, tree=tree, subst.mod="REV")
write.tm(neutralModAll , paste(sOutStem, ".allspp.mod", sep=""));

