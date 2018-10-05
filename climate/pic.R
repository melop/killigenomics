setwd("./bioclim");
#.libPaths("/beegfs/common/software/2017/modules/software/rlang/3.5.1/lib64/R/library")

library(caper)
library(ape)
library(geiger)
library(nlme)
library(phytools)


sTree <- "in.tre.txt";
sData <- "out.withMTGenomeSize.txt"

#sTestVar1 <- "Genomesize_BUSCO" #"Genomesize_BUSCO"; #"Genomesize_BUSCO"; #"Log_GenomeSize_nowarning" ; #"Log_GenomeSize_nowarning"; "Log_GenomeSize"; #GenomeSize_nowarning"; #arrTestVars1[nTestVar]; GenomeSize
sTestVar1 <- "mt_genome_size";

#sTestVar1 <- "median_relax_k";
#sTestVar2 <- "median_relax_k"; #precip_month_11"; #arrTestVars2[nTestVar];
#sTestVar2 <- "Annual"; #precip_month_11"; #arrTestVars2[nTestVar];
#sTestVar2 <- "bioclim_bio_14"; #precip_month_11"; #arrTestVars2[nTestVar];
#sTestVar2 <- "bioclim_bio_17"; #precip_month_11"; #arrTestVars2[nTestVar];


#  sTestVar1 <- "bioclim_bio_1"; 
#  sTestVar1 <- "bioclim_bio_2";
#  sTestVar1 <- "bioclim_bio_3";
#  sTestVar1 <- "bioclim_bio_4";
#  sTestVar1 <- "bioclim_bio_5";
#  sTestVar1 <- "bioclim_bio_6";
#  sTestVar1 <- "bioclim_bio_7";
#  sTestVar1 <- "bioclim_bio_8";
#  sTestVar1 <- "bioclim_bio_9";
#  sTestVar1 <- "bioclim_bio_10";
#  sTestVar1 <- "bioclim_bio_11";
#  sTestVar1 <- "bioclim_bio_12";
#  sTestVar1 <- "bioclim_bio_13";
#  sTestVar1 <- "bioclim_bio_14";
#  sTestVar1 <- "bioclim_bio_15";
#  sTestVar1 <- "bioclim_bio_16";
#  sTestVar1 <- "bioclim_bio_17";
#  sTestVar1 <- "bioclim_bio_18";
# sTestVar1 <- "bioclim_bio_19";
# 
# sTestVar2 <- "Annual";

# 
  sTestVar2 <- "bioclim_bio_1"; 
  sTestVar2 <- "bioclim_bio_2";
  sTestVar2 <- "bioclim_bio_3";
  sTestVar2 <- "bioclim_bio_4";
  sTestVar2 <- "bioclim_bio_5";
  sTestVar2 <- "bioclim_bio_6";
  sTestVar2 <- "bioclim_bio_7";
  sTestVar2 <- "bioclim_bio_8";
  sTestVar2 <- "bioclim_bio_9";
  sTestVar2 <- "bioclim_bio_10";
  sTestVar2 <- "bioclim_bio_11";
  sTestVar2 <- "bioclim_bio_12";
  sTestVar2 <- "bioclim_bio_13";
  sTestVar2 <- "bioclim_bio_14";
  sTestVar2 <- "bioclim_bio_15";
  sTestVar2 <- "bioclim_bio_16";
  sTestVar2 <- "bioclim_bio_17";
  sTestVar2 <- "bioclim_bio_18";
  sTestVar2 <- "bioclim_bio_19";

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter


#arrTestVars1 <- c("mk_alpha", "mk_alpha", "CAI_coef_Codon" , "CAI_coef_AA", "CAI_coef_AA");
#arrTestVars2 <- c("CAI_coef_Codon", "pi", "pi", "pi", "CAI_coef_Codon");

strcount <- function(x, pattern, split){
  
  unlist(lapply(
    strsplit(x, split),
    function(z) na.omit(length(grep(pattern, z)))
  ))
  
}



#for(nTestVar in 1:length(arrTestVars1)) {

sOut <- paste(sData , basename(sTree), "var1", sTestVar1, "var2", sTestVar2, ".txt", sep="_");

sink(file = sOut);

cat("Doing ", sTestVar1 , "vs.", sTestVar2, "\n");


oTree <- read.tree(file = sTree);
dat <- read.table(sData , header=T, sep="\t");
dat$Log_Genomesize_BUSCO <- log10(dat$Genomesize_BUSCO);
#dat$Annual <- factor(x=as.character(dat$Annual), levels=c('NonAnnual', 'SemiAnnual', 'Annual'), ordered=T )
boxplot(dat$Genomesize_BUSCO ~ dat$Annual, xlab="", ylab="Genome Size", col=c('blue', 'green', 'red'))

datClade1 <- dat[ dat$Genus %in% c('Aphyosemion', 'Nothobranchius', 'Fundulopanchax', 'Fundulosoma', "Pronothobranchius"), ];
boxplot(datClade1$Genomesize_BUSCO ~ datClade1$Annual, xlab="", ylab="Genome Size", col=c('blue', 'green', 'red'))


datClade2 <- dat[ dat$Genus %in% c('Scriptaphyosemion', 'Callopanchax'), ];
boxplot(datClade2$Genomesize_BUSCO ~ datClade2$Annual, xlab="", ylab="Genome Size", col=c('blue', 'green', 'red'))



arrNames <- colnames(dat);
arrNames[1] <- "taxa";
colnames(dat) <- arrNames;

dat <- dat[, c("taxa", sTestVar1, sTestVar2)];
#dat[paste(sTestVar2, "_ordered",sep=""), ] <- factor(x = dat[sTestVar2, ], levels = c("NonAnnual", "SemiAnnual", "Annual"), ordered = T  )


cdat <- NULL;
cdat <- comparative.data(data=dat, phy=oTree, names.col="taxa" , vcv=TRUE, vcv.dim=3)



######## method 1 
cat("raw lm without phylogenetic correction\n");
summary(lm(as.formula( paste( sTestVar1, " ~ ", sTestVar2 )) , data=cdat$data) );

cat("PGLS method, brownian\n");
pglsModel2 <- gls( as.formula( paste( sTestVar1, " ~ ", sTestVar2 )) , correlation = corBrownian(phy = cdat$phy),
                   data = cdat$data, method = "ML")
#print(anova(pglsModel2));
print(summary(pglsModel2));



######## method 1 end


cat("PGLS method with model selection\n");


arrParams <- c("kappa", "lambda" , "delta");
#perform model selection with AIC
arrFitModels <- list();
arrAIC <- c();
arrBIC <- c();
for (nFreeParams in 0:3) { #test combinations of parameters
  matAllComb <- combn(arrParams , nFreeParams);
  arrML <- list();
  
  
  for (nCol in 1:ncol(matAllComb) ) {
    arrMLTrue <- matAllComb[,nCol];
    arrML$kappa <- 1;
    arrML$lambda <- 1;
    arrML$delta <- 1;
    if (length(arrMLTrue) > 0) {
      arrML[arrMLTrue] <- "ML";
    }
    
    sID <- paste("ML_", arrMLTrue, collapse = "_", sep="");
    arrFitModels[[sID]] <- pgls( as.formula( paste(sTestVar1, " ~ ", sTestVar2)), cdat, kappa = arrML[["kappa"]], lambda=arrML[["lambda"]], delta=arrML[["delta"]])
    # summary(arrFitModels[[sID]]);
    arrAIC[sID] <- AIC(arrFitModels[[sID]]);
    arrBIC[sID] <- BIC(arrFitModels[[sID]]);
    
  }
  
}


print(sort(arrAIC));
nMinAIC <- min(arrAIC);
arrAICProb <- exp((nMinAIC - arrAIC)/2);
print(sort(arrAICProb));
nSuboptimalAIC <- nMinAIC+4 ; # 2 AIC units should be considered no sig. difference
arrAICCandidates <- sort(arrAIC[arrAIC <=nSuboptimalAIC]);
arrModelParams <-strcount(names(arrAICCandidates), "", "ML") - 1;
arrAICCandidatesFewParams <- which( arrModelParams ==  min(arrModelParams) );
sSelectedModel <- names( arrAICCandidates[arrAICCandidatesFewParams[1]]);
cat("Selected Model (AIC) is:", sSelectedModel, "\n");
print(summary(arrFitModels[[sSelectedModel ]]));
plot(arrFitModels[[sSelectedModel ]]);


print(sort(arrBIC));
nMinBIC <- min(arrBIC);
arrBICProb <- exp((nMinBIC - arrBIC)/2);
print(sort(arrBICProb));
nSuboptimalBIC <- nMinBIC+4 ; # 2 AIC units should be considered no sig. difference
arrBICCandidates <- sort(arrBIC[arrBIC <=nSuboptimalBIC]);
arrModelParams <-strcount(names(arrBICCandidates), "", "ML") - 1;
arrBICCandidatesFewParams <- which( arrModelParams ==  min(arrModelParams) );
sSelectedModel <- names( arrBICCandidates[arrBICCandidatesFewParams[1]]);
cat("Selected Model (BIC) is:", sSelectedModel, "\n");

print(summary(arrFitModels[[sSelectedModel ]]));


#plot( pgls.profile(mod, "kappa") );
#plot( pgls.profile(mod, "lambda") );
#plot( pgls.profile(mod, "delta") );

#plot(mod)
#cat("Independent contrast");

#cdat$data$AnnualOrdered <- factor(x = cdat$data$Annual, levels=c("NonAnnual", "SemiAnnual", "Annual"), ordered=T);

#crunchMod <- crunch(as.formula( paste(sTestVar1, " ~ ", paste(sTestVar2, "Ordered",sep="") )), data=cdat, factor.action = 'allow');
#print(summary(crunchMod));

#print(caic.diagnostics(crunchMod));
#plot(crunchMod$contrast.data$contr$explanatory , crunchMod$contrast.data$contr$response, xlab=sTestVar2 , ylab=sTestVar1 )
sink();
#}
