#This script extracts the bioclim variables from the .bil files given GPS coordinates
#The in.txt file should contain two columns, one named "longitude",the other "latitude"
#This script appends columns for each bioclim variables.

setwd("./bioclim/");
library("raster");

dat <- read.table("in.txt", header=T, sep="\t");
coordinates(dat) <- c("longitude", "latitude");

#points(dat$longitude, dat$latitude)

for(nBioClimFactor in 1:19) {
  sBioClimFactor <- nBioClimFactor;

  cat("Doing factor " , sBioClimFactor , "\n");
  arrBils <- Sys.glob(paste("bioclim_2.5/bio", sBioClimFactor,".bil", sep='') );
  
  arrRasters <- list();
  
  for(i in 1:length(arrBils) ) {
    arrRasters[[i]] <- raster(arrBils[i]);
  }
  

  oRet <- extract(arrRasters[[1]], dat, method='simple', buffer=1000, fun=mean, df=F) ## If df=TRUE, the results are returned as a 'dataframe'
  dat[[ paste("bioclim_bio_" , sBioClimFactor, sep='') ]] <- as.numeric( oRet);

}

write.table(dat, file="out.withBioClim.txt" , row.names = F, col.names = T, quote=F , sep="\t");

dat <- read.table("out.withBioClim.txt" , header=T, sep="\t");
lsLM <- list();
datMonth <- data.frame(Month=numeric(), R2=numeric(), P=numeric(), COEF=numeric());
for (nBioClimFactor in 1:19) {
  sBioClimFactor <- nBioClimFactor;

  oLm <- summary(lm( dat$Genomesize_BUSCO ~ dat[[ paste("bioclim_bio_" , sBioClimFactor, sep='') ]] ));
  lsLM[[nBioClimFactor]] <- oLm;
  datMonth <- rbind(datMonth, list(Month=nBioClimFactor , R2=oLm$r.squared, P=oLm$coefficients[2,4], COEF=oLm$coefficients[2,1] ) );
}

plot(datMonth$Month, datMonth$R2, xlab="BioClim Bio Factor", ylab= "R2 GenomeSize (bp) ~ BioFactor", type="l");
plot(datMonth$Month, datMonth$P, xlab="BioClim Bio Factor", ylab = "P-value GenomeSize (bp) ~ BioFactor", type="l");
abline(h=0.05, col='red')

#################### phylogenetic corrected version ###################
library(caper)
library(ape)
library(geiger)
library(nlme)
library(phytools)

dat <- read.table("out.withBioClim.txt" , header=T, sep="\t");


sTree <- "../../in.tre.txt";
oTree <- read.tree(file = sTree);

dat$taxa <- dat$Species;

lsLM <- list();
datMonth <- data.frame(Month=numeric(), R2=numeric(), P=numeric(), COEF=numeric());
for (nBioClimFactor in 1:19) {
  sBioClimFactor <- nBioClimFactor;

   sBioClimFactorCol <- paste("bioclim_bio_" , sBioClimFactor, sep='');
  datSub <- dat[, c("taxa", "Genomesize_BUSCO", sBioClimFactorCol )];
  cdat <- NULL;
  cdat <- comparative.data(data=datSub, phy=oTree, names.col="taxa" , vcv=TRUE, vcv.dim=3)
  
   
  pglsModel <- gls( as.formula( paste(  "Genomesize_BUSCO ~ ", paste("bioclim_bio_" , sBioClimFactor, sep='') )) , correlation = corBrownian(phy = cdat$phy),
    data = cdat$data, method = "ML")
  oLm <- anova(pglsModel);

  lsLM[[nBioClimFactor]] <- oLm;
  datMonth <- rbind(datMonth, list(Month=nBioClimFactor , P=oLm$`p-value`[2], COEF=pglsModel$coefficients[2]) );
}

plot(datMonth$Month, datMonth$P, xlab="BioClim Bio Factor", ylab = "P-value GenomeSize (bp) ~ BioFactor", type="l");
abline(h=0.05, col='red')




