setwd("/beegfs/group_dv/home/RCui/killifish_genomes/assemble_mt/mitos_annotation/alignments/aln_Callopanchax_nad2-0/raxml_syn");

library(phytools);
library(ggplot2)

oTree <- read.tree(file = "RAxML_result.syn");
oOutGroup <- fastMRCA(oTree, 'PLP', 'APL');
oTree <- reroot(tree = oTree, node.number = oOutGroup)

lsCompare <- list();
lsCompare[['Western']] <- list();
lsCompare[['Western']][['Annual']] <- c('CSB', 'CTO', 'CMR');
lsCompare[['Western']][['NonAnnual']] <- c('SSC', 'SGG', 'SBT', 'SCV');

lsCompare[['Eastern']] <- list();
lsCompare[['Eastern']][['Annual']] <- c('NKF', 'NVG', 'NVS', 'NFS', 'NOR','NFZ','NRC', 'NOC');
lsCompare[['Eastern']][['NonAnnual']] <- c('ACM', 'ACL', 'ACG', 'ACY', 'AGM', 'AAU', 'AKM');

datRet <- data.frame(Sp=character(), Clade=character(), Annualism=character(), NodeHeigh=c(), stringsAsFactors=FALSE);

for (sClade in names(lsCompare) ) {
  for(sAnnual in names(lsCompare[[sClade]]) ) {
    for(sSp in lsCompare[[sClade]][[sAnnual]]) {
      nHeight <- fastHeight(oTree, sSp, sSp);
      datRet <- rbind( datRet, data.frame(Sp=sSp, Clade=sClade,Annualism=sAnnual, NodeHeigh=nHeight, stringsAsFactors=F) );
    }
  }
  
}

pdf(file="annual_vs_nonannual_syn_branchlen.pdf", width=5, height=4)
par(mfrow=c(1,2))

for (sClade in names(lsCompare) ) {
  boxplot(NodeHeigh ~ Annualism ,data = datRet[datRet$Clade==sClade,], main=sClade, ylab="Synonymous site branch length", col=c('red','blue'));
}
dev.off();