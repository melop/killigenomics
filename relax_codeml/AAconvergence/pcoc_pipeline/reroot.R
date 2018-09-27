library(ape)
library(phytools)
args = commandArgs(trailingOnly=TRUE)

arrOutgroupSpp <- c('PLP' , 'APL');

oTree <- read.tree(file = args[1] );

if (is.null(oTree)) {
  cat("Error reading tree \n");
}


oTree <- unroot( oTree );
arrOutgroupSppFound <- intersect(arrOutgroupSpp, oTree$tip.label);
if (length(arrOutgroupSppFound) ==0) {
  cat("ERROR: All outgroup species are missing, cannot root");
}

if (length(arrOutgroupSppFound ) == 1) {
  cat( " only one outgroup \n");
  
  oTree <- reroot(tree = oTree, node.number = which(oTree$tip.label == arrOutgroupSppFound[1]));
} else {
  oTree <- reroot(tree = oTree, node.number = findMRCA(oTree, tips=arrOutgroupSppFound,type = 'node' ));
}

sOutTree <- paste(args[1], '.reroot.tre', sep="") ;
if (length(args) >=2) {
	sOutTree <- args[2];
}

write.tree(oTree , sOutTree)
#oTree <- read.tree(file = paste(args[1], '.reroot.tre', sep="") );

#plot(oTree)
#nodelabels()
#tiplabels()
