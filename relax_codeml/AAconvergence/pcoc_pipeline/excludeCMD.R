setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/AAconvergence/pcoc_pipeline");

sPCOC <- "pcoc_out/pcoc.sig.txt";
sPCOCexcludeCMD <- "pcoc_out/pcoc.sig.excludeCMD.txt";

sCMDList <- "../CMDlist/CMDlist.txt";

datPCOCSig <- read.table(sPCOC , header=T, sep="\t", quote="");
datCMDList <- read.table(sCMDList, header=F, sep="\t");
datCMDList$V2 <- datCMDList$V2 + 1; # adjust to 1-based index, because that's what the pcoc output is

datCMDList$Pos <- paste(datCMDList$V1, datCMDList$V2 , sep=":");

datPCOCSigFilter <- datPCOCSig[ ! (paste(datPCOCSig$OrthoID, datPCOCSig$AASite, sep=":") %in% datCMDList$Pos),  ];
nrow(datPCOCSig);
nrow(datPCOCSigFilter);
write.table(datPCOCSigFilter , file=sPCOCexcludeCMD, sep="\t", row.names = F, col.names = T, quote=F);
