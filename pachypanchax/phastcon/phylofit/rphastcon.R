setwd("/beegfs/group_dv/home/RCui/killifish_genomes/phastcon/phylofit/NFZ_ref");

library(rphast)
args = commandArgs(trailingOnly=TRUE)
nChr <- args[1]

sChr <- paste("chr", nChr,sep="");

sRefSp <- 'NFZ'
tree <- "((austrofundulus,kryptolebias),(PLP,(CTO,(AAU,(NOR,NFZ)))))";

arrAnnuals <- c( 'NOR', 'CTO', 'austrofundulus' );
arrNonAnnuals <- c('AAU', 'PLP', 'kryptolebias');

treeAnnuals <- "(austrofundulus,(CTO,(NOR,NFZ)))";



fnPlotSim <- function(nonParaPhyloP,obsPhyloP) {
  layout(matrix(c(1,2), nrow=2), heights=c(0.7, 0.3))
  par(mar=c(4.5, 4, 4, 2), mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
  qqplot(nonParaPhyloP$lnlratio,obsPhyloP$lnlratio,
         xlim=c(0,15),ylim=c(0,15), xlab="Simulated likelihood ratio",
         ylab="Observed likelihood ratio")
  abline(0, 1, lty=2)
  par(mar=c(4,4,1,2))
  plot(density(obsPhyloP$lnlratio,adjust=3), lty=1,xlim=c(0,15),
       xlab="Likelihood Ratio",
       ylab="Density",main="", col="red")
  lines(density(nonParaPhyloP$lnlratio,adjust=3), lty=1,
        col="black",xlim=c(0,15))
}

empirical.pval <- function(x, dist) {
   sum(x <= dist)/length(dist)
}



feats <- read.feat(paste(sChr, "_cds.gff", sep=""))
alnAnnual <- read.msa(filename = paste("../../cactus/split_NFZ_maf/NFZ_ref_5spp_", sChr, ".maf" , sep=""), format = "MAF", ordered = T, seqnames = c(sRefSp, arrAnnuals) )

names(alnAnnual)
feats$seqname <- sRefSp;
table(feats$feature)

aln4d <- get4d.msa(alnAnnual, feats)
neutralMod <- phyloFit(aln4d, tree=tree, subst.mod="REV")
pcEM_Annual <- phastCons(alnAnnual, neutralMod, viterbi=TRUE, estimate.transitions=TRUE, estimate.expected.length = TRUE, estimate.rho = TRUE)

alnAll <- read.msa(filename = paste("../../cactus/split_NFZ_maf/NFZ_ref_5spp_", sChr, ".maf" , sep=""), format = "MAF", ordered = T )
aln4dAll <- get4d.msa(alnAll, feats)
neutralModAll <- phyloFit(aln4dAll, tree=tree, subst.mod="REV")

hasNonAnnuals <- informative.regions.msa(alnAll, 3, arrNonAnnuals);
hasAtLeast6 <- informative.regions.msa(alnAll, 6)

informativeElements <- coverage.feat(pcEM_Annual$most.conserved, hasNonAnnuals, hasAtLeast6, get.feats=TRUE)
coverage.feat(informativeElements)/coverage.feat(pcEM_Annual$most.conserved)

splitLength <- 50

splitElements <- split.feat(informativeElements, f=splitLength, drop=TRUE)

obsPhyloP_PLP <- phyloP(neutralModAll, msa=alnAll, mode="ACC", features=splitElements, subtree="PLP")
obsPhyloP_AAU <- phyloP(neutralModAll, msa=alnAll, mode="ACC", features=splitElements, subtree="AAU")

#do non parametric bootstrap:
elementAlign <- extract.feature.msa(copy.msa(alnAll), informativeElements, pointer.only=TRUE)
nrep <- 100000
simMsa <- sample.msa(elementAlign, nrep*splitLength, replace=TRUE)
# produce features allowing long alignment to be interpreted as
 # concatenation of shorter alignments
 startIdx <- seq(from=1, by=splitLength, length.out=nrep)
features <- feat(seqname=names.msa(simMsa)[1], src="sim", feat=".", start=startIdx, end=startIdx+splitLength-1)
nonParaPhyloP_PLP <- phyloP(neutralModAll, msa=simMsa, mode="ACC", features=features, subtree="PLP")
nonParaPhyloP_AAU <- phyloP(neutralModAll, msa=simMsa, mode="ACC", features=features, subtree="AAU")

fnPlotSim(nonParaPhyloP = nonParaPhyloP_PLP, obsPhyloP = obsPhyloP_PLP )
fnPlotSim(nonParaPhyloP = nonParaPhyloP_AAU, obsPhyloP = obsPhyloP_AAU )

nonParaPval_PLP <- sapply(obsPhyloP_PLP$lnlratio, empirical.pval, nonParaPhyloP_PLP$lnlratio)
nonParaFDR_PLP <- p.adjust(nonParaPval_PLP, method="BH")

nonParaPval_AAU <- sapply(obsPhyloP_AAU$lnlratio, empirical.pval, nonParaPhyloP_AAU$lnlratio)
nonParaFDR_AAU <- p.adjust(nonParaPval_AAU, method="BH")

nonParaSigFeats_PLP <- splitElements[nonParaFDR_PLP < 0.05,]
nrow(nonParaSigFeats_PLP)
nonParaSigFeats_PLP$feature <- "PLPacc"
nonParaSigFeats_PLP$score <- obsPhyloP$lnlratio[nonParaFDR_PLP < 0.05]
nonParaSigFeats_PLP$seqname <- sChr; 
write.feat(nonParaSigFeats_PLP, paste("PLP_NFZref_",sChr,".gff",sep="") )


nonParaSigFeats_AAU <- splitElements[nonParaFDR_AAU< 0.05,]
nrow(nonParaSigFeats_AAU)
nonParaSigFeats_AAU$feature <- "AAUacc"
nonParaSigFeats_AAU$score <- obsPhyloP$lnlratio[nonParaFDR_AAU < 0.05]
nonParaSigFeats_AAU$seqname <- sChr; 
write.feat(nonParaSigFeats_AAU, paste("AAU_NFZref_",sChr,".gff",sep="") )


