setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline");

sPCOC <- "pcoc_out/pcoc.sig.txt";
sAnnotation <- "../annotate_aa_changes/with.deltaLK.txt";

sPCOC <- "pcoc_codemlsim_out/pcoc.sig.txt";
sAnnotation <- "../../AAconvergence_sim/annotate_aa_changes/with.deltaLK.txt";


datPCOC <- read.table(sPCOC, header=T, sep="\t", quote="");
datAAAnnotation <- read.table(sAnnotation, header=T, sep="\t", quote="");

datAAAnnotation$Pos = datAAAnnotation$Pos + 1; #adjust 0 based to 1 base
datMerge <- merge(datAAAnnotation , datPCOC[,c(1:8) ], by.x=c('OrthoID', 'Pos'), by.y=c('OrthoID', 'AASite') );
write.table(datMerge, paste(sPCOC, ".AAAnnote.overlap.txt",sep=""), col.names = T, row.names = F, quote=F, sep="\t" );
