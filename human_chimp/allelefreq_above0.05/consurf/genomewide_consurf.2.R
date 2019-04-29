args = commandArgs(trailingOnly=TRUE)
sPop <- args[1];
#sPop <- "ptt";
#sPop <- "pts";


sOutNonSyn <- paste(sPop, ".nonsyn.AF.Consurf.txt", sep="");
sOutSyn <- paste(sPop, ".syn.AF.Consurf.txt", sep="");
sConsurfScores <- pipe(paste("zcat -f", " /beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/consurf/human_GRCH37/consurf_scores_flat.txt.gz") )

s4DSites <- paste("../classify_sites/ret_pop_", sPop, "_syn_table.txt" ,sep=""); #these are synonymous sites, not 4D, for 4D , see ../MKTEST
sNonSynSites <- paste("../classify_sites/ret_pop_", sPop, "_nonsyn_table.txt" ,sep=""); #these are non-syn sites computed based on actual mutations.

cat("Loading synonymous sites...\n");
dat4DSites <- read.table( pipe(paste("zcat -f ", s4DSites)) , header=F);
dat4DSites$coord <- paste(dat4DSites$V1, dat4DSites$V3, sep=":"); # 1-based coord.
dat4DSites <- dat4DSites[ , c(8,9,10) ];
colnames(dat4DSites) <- c('AncCodon', 'DerCodon', 'coord');

cat("Loading non-synonymous sites...\n");
datNonSynSites <- read.table(pipe(paste("zcat -f ", sNonSynSites))  , header=F);
datNonSynSites$coord <- paste(datNonSynSites$V1, datNonSynSites$V3, sep=":"); # 1-based coord.
datNonSynSites <- datNonSynSites[ , c(8,9,10) ];
colnames(datNonSynSites) <- c('AncCodon', 'DerCodon', 'coord');

sSNPFolder <- "../"; # <- JUST CHANGE THIS FOLDER FOR DIFFERENT POPS!

sPolySites <- pipe(paste("n=0; for i in ", sSNPFolder , "/polarized/", sPop,".chr*.txt.gz; do n=$(( n+ 1 )); if (( $n == 1)); then zcat $i; else zcat $i | tail -n +2; fi; done", sep="")); #polymorphic sites
sDivergentSites <- pipe(paste("n=0; for i in ", sSNPFolder , "/polarized/der_", sPop,".chr*.txt.gz; do n=$(( n+ 1 )); if (( $n == 1)); then zcat $i; else zcat $i | tail -n +2; fi; done", sep="")); #polymorphic sites



cat("Loading polymorphic sites...\n");
datPolySites <- read.table(sPolySites, header=T, colClasses = c("factor",	"numeric", "factor", "factor", 	rep("numeric", 8)	) );
cat("Loading divergent sites...\n");
datDivergentSites <- read.table(sDivergentSites, header=T, colClasses = c("factor",	"numeric", "factor", "factor", 	rep("numeric", 3)	) );


cat("Reading Consurf scores...\n");
datConsurf <- read.table(sConsurfScores, header=F, sep="\t", colClasses = c("factor", "numeric", "numeric", "factor") );
datConsurf$coord <- paste(datConsurf$V1, datConsurf$V2, sep=":");


datAllFreq <- datPolySites[ , c("scaffold", "site" , "derived_freq") ];
datDivergentFreq <- datDivergentSites[, c("scaffold", "site" ) ];
datDivergentFreq$derived_freq <- 1;

datAllFreq <- rbind(datAllFreq , datDivergentFreq);
datAllFreq$coord <- paste(datAllFreq$scaffold, datAllFreq$site, sep=":");

datAllFreqFilter <- datAllFreq[datAllFreq$coord %in% datConsurf$coord, ];
datAllFreqFilter <- merge(datAllFreqFilter[, c('coord', 'derived_freq')], datConsurf[, c('coord', 'V3', 'V4')], by='coord')
names(datAllFreqFilter)[3] <- "ConsurfScore";
names(datAllFreqFilter)[4] <- "Gene";


cat(nrow(datAllFreqFilter),"\n");

datNonSynAllFreq <- unique(merge(datAllFreqFilter, datNonSynSites, by="coord", all.y=F));
datSynAllFreq <- unique(merge(datAllFreqFilter, dat4DSites, by="coord", all.y=F));

cat(nrow(datNonSynAllFreq), "\n");
cat(nrow(datSynAllFreq), "\n");


write.table(datNonSynAllFreq, file=sOutNonSyn, row.names = F, col.names = T, sep="\t", quote=F);
write.table(datSynAllFreq, file=sOutSyn, row.names = F, col.names = T, sep="\t", quote=F);



