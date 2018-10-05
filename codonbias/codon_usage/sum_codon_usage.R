#this script computes CAI for each gene
sCodonUsageCounts <- 'codon_use_counts.txt';#"codon_use_counts.txt";
#sCodonUsageCounts <- "codon_use_counts_frameshift_1_revcomp.txt"
#sCodonUsageCounts <- "codon_use_counts_frameshift_1.txt"
#sCodonUsageCounts <- "codon_use_counts_frameshift_2_revcomp.txt"
#sCodonUsageCounts <- "codon_use_counts_frameshift_2.txt"
#sCodonUsageCounts <- "codon_use_counts_revcomp.txt"



sGeneExpTable <- "../../gene_expression_levels/high.expression.level.txt"; #use this table to select genes that are highly expressed

sOut <- "Codon_usage_highlyexpressed_genes";
#sOut <- "Codon_usage_frameshift_1_revcomp_highlyexpressed_genes";
#sOut <- "Codon_usage_frameshift_1_highlyexpressed_genes";
#sOut <- "Codon_usage_frameshift_2_revcomp_highlyexpressed_genes";
#sOut <- "Codon_usage_frameshift_2_highlyexpressed_genes";
#sOut <- "Codon_usage_revcomp_highlyexpressed_genes";

datCodonUsageCounts <- read.table(sCodonUsageCounts , header=T, stringsAsFactors = F, quote='', sep="\t");
datExp <- read.table(sGeneExpTable , header=T, stringsAsFactors =F, quote='', sep="\t");


arrSelectedGenes <- datExp$OrthoID;

datCodonUsageSelected <- datCodonUsageCounts[datCodonUsageCounts$GeneName %in% arrSelectedGenes, ];

arrCodons <- unique(datCodonUsageCounts$Codon);
arrSpecies <- colnames(datCodonUsageCounts);
arrSpecies <- arrSpecies[3:length(arrSpecies)];

arrAppend <- as.list(rep(F, length(arrSpecies)));
names(arrAppend) <- arrSpecies;
for(sCodon in arrCodons) {
  datSubset <- datCodonUsageSelected[datCodonUsageSelected$Codon == sCodon, ];
  arrSums <- colSums(datSubset[, 3:ncol(datSubset)]);
  for (sSp in arrSpecies) {
   cat( sCodon, "\t", arrSums[[sSp]], "\n", sep="", file=paste(sOut, "_", sSp, ".txt", sep=""), append=arrAppend[[sSp]]  );
    if (!arrAppend[[sSp]] ) {
      arrAppend[[sSp]]  <- T;
    }
  }
  
}
