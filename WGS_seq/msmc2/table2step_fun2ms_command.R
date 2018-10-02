setwd("../WGS_seq/bwamap/NORv1.2LG/msmc2/");
### txt2step_fun2ms_command ###
#takes txt table file from stdin, reads Popsize and generation, calculates step function
#and converts it into ms command

sSp <- "ORT";
nDivergenceTime <- 205.35e3;

sSp <- "RAC";
nDivergenceTime <- 258.2e3;
nsamp <- 100
nrep <- 1
nSamplesPop1 <- 51
mu <- 2.63E-9
locuslength <- 1;

args <- commandArgs(trailingOnly=TRUE)
sSp <- args[1]
nDivergenceTime <- as.numeric( args[2])
nsamp <- as.integer(args[3])
nrep <- as.integer(args[4])
nSamplesPop1 <- as.integer(args[5])
mu <- as.numeric(args[6])
locuslength <- as.numeric(args[7])

nSamplesPop2 <- as.numeric(nsamp) - as.numeric(nSamplesPop1);

sInDir <- "msmc2ret2/"
arrPops <- paste(sSp, c('DRY', "WET"), sep="" );
arrInFiles <- paste(sInDir, "/", arrPops, ".final.txt", sep="");

datMSMC1 <- read.table(arrInFiles[1], header=T, sep="\t" );
datMSMC1$gen <- datMSMC1$left_time_boundary/mu;
datMSMC1$popsize <- (1/datMSMC1$lambda)/(2*mu);

datMSMC2 <- read.table(arrInFiles[2], header=T, sep="\t" );
datMSMC2$gen <- datMSMC2$left_time_boundary/mu;
datMSMC2$popsize <- (1/datMSMC2$lambda)/(2*mu);

Popsize_at_Gen_t <- data.frame( Pop1Gen=datMSMC1$gen, Pop1Size=datMSMC1$popsize, Pop2Gen=datMSMC2$gen, Pop2Size=datMSMC2$popsize )

colnames(Popsize_at_Gen_t)[1] <- "Pop1_Gen"
colnames(Popsize_at_Gen_t)[2] <- "Pop1_Popsize"
colnames(Popsize_at_Gen_t)[3] <- "Pop2_Gen"
colnames(Popsize_at_Gen_t)[4] <- "Pop2_Popsize"

print(Popsize_at_Gen_t);




##convert table informations into ms parameters##

# fuction for calculating the populationsize for every generation of the step function relative to the N0
fnNirelativetoN0 <- function(PopsizeN0,PopsizeNi){

  return(PopsizeNi/PopsizeN0)
}

# function for calculating the time points of the coalescent model in ms by deviding every generation by N0
fnTimeFromPresent <- function(Gen,PopsizeN0) {
  return(Gen/(4*PopsizeN0))
}

# calculation of the ms parameters
Theta_ms <- (4*datMSMC1$popsize[1]*as.numeric(mu)*as.numeric(locuslength))

#pdf(paste(arrPops[1], arrPops[2], "smc2ms.pdf",sep="_"), width=6, height=6);
results_TimeFromPresent_Pop1 <- fnTimeFromPresent(Gen = datMSMC1$gen ,PopsizeN0 = datMSMC1$popsize[1] )
results_NitoN0_Pop1 <- fnNirelativetoN0(PopsizeN0 = datMSMC1$popsize[1] ,PopsizeNi = datMSMC1$popsize);
plot(log10(datMSMC1$gen) , log10(datMSMC1$popsize) );

results_TimeFromPresent_Pop2 <- fnTimeFromPresent(Gen = datMSMC2$gen ,PopsizeN0 = datMSMC1$popsize[1] )
results_NitoN0_Pop2 <- fnNirelativetoN0(PopsizeN0 = datMSMC1$popsize[1] ,PopsizeNi = datMSMC2$popsize);
plot(log10(datMSMC2$gen) , log10(datMSMC2$popsize) );

tCoalescent_Event= as.numeric(nDivergenceTime)/(4*datMSMC1$popsize[1])

##  writing the ms command line arguments which are outputted on screen by the cat function ##
Demography_values_Pop1 <- paste(paste("-en ",results_TimeFromPresent_Pop1,sep= " ")," 1 ",results_NitoN0_Pop1,  sep = " ")

# get rif of all populationsizes of Pop2 after coalescent event with Pop 1
results_TimeFromPresent_Pop2 <- results_TimeFromPresent_Pop2[results_TimeFromPresent_Pop2<=tCoalescent_Event]

Demography_values_Pop2 <- paste(paste("-en ",results_TimeFromPresent_Pop2,sep= " ")," 2 ",results_NitoN0_Pop2[1:length(results_TimeFromPresent_Pop2)],  sep = " ")


cat("ms ",nsamp," ",nrep," -t ",Theta_ms," -I 2 ",nSamplesPop1," ",nSamplesPop2," ",Demography_values_Pop1," ",Demography_values_Pop2," -ej  ",tCoalescent_Event,"2 1")






