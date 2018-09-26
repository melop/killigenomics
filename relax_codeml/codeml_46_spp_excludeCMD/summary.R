#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/rerun_exclude_MNM");
#summarize current and previous runs.
args = commandArgs(trailingOnly=TRUE)
sGenus <- "Callopanchax";

if (length(args) == 1) {
  sGenus <- args[1];
}

cat("summarizing ", sGenus, "\n");
sOut <- paste("sum_" , sGenus, ".txt" , sep="")
sErrOut <- paste("sum_" , sGenus, "_err.txt" , sep="")

sOld <- paste("../Codeml_" , sGenus, "/ret*", sep="");
sNew <- paste("./Codeml_" , sGenus, "/ret*", sep="");

system(paste("cat ", sOld ," | grep -v -P ",'"Success\\t"', " | sort -t $'\\t' -k 5n,5n  > ",sErrOut , sep=""))
system(paste("wc -l ",sErrOut , sep=""))


oPOld <- pipe(paste("cat ", sOld ," | grep -P ",'"Success\\t"', " | sort -t $'\\t' -k 5n,5n  ", sep="") );
oPNew <- pipe(paste("cat ", sNew ," | grep -P ",'"Success\\t"', " | sort -t $'\\t' -k 5n,5n  ", sep="") );

datOld <- read.table(oPOld , header = F, sep="\t", quote = "", fill = T, stringsAsFactors = F)
datNew <- read.table(oPNew , header = F, sep="\t", quote = "", fill = T, stringsAsFactors = F)

datUnChanged <- datOld[!(datOld$V1 %in% datNew$V1), ]

cat("Unchanged tests: ", nrow(datUnChanged), "\n" );
datNew <- rbind(datNew, datUnChanged)
datNew <- datNew[order(datNew$V5) , ]
datNew$FDR <- p.adjust(datNew$V5, method="fdr")

write.table(datNew , file=sOut, col.names = F, row.names = F, quote=F, sep="\t");
system(paste("wc -l ",sOut , sep=""))
