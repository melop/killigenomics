### getting the Omega K values for each Omega categorie for significantly relaxed, positive selected and intensified genes ###
### only the K distribution needed is calculated ####

### read in the RELAX output ###
Input <- read.table("/beegfs/group_dv/home/LIasi/relax_output/sum_Nothobranchius.txt", header=F, fill = T, stringsAsFactors = F, sep = "\t", quote="")

### filter for complete cases ###
datIn <- Input[complete.cases(Input),]
### filter for significant results ###
datIn <- subset(datIn, V5 <0.05) 
### filter for relaxed genes with K smaller One ###
datIn_sign_Relaxed <- subset(datIn, V9 <1) 

#### calculate all the K's for Omega class 0 for significantly relaxed genes by the log(Omega-of-class-1-FG of base Omega-of-class-1-BG) ###
K_smaller_One <- log(datIn_sign_Relaxed$V17,datIn_sign_Relaxed$V11)
K_smaller_One <- K_smaller_One[complete.cases(K_smaller_One)]
K_smaller_One_Relaxed <- data.frame(K_smaller_One)

#### calculate all the K's for Omega class 1 for significantly relaxed genes ###
K_equal_One <-log(datIn_sign_Relaxed$V19,datIn_sign_Relaxed$V13)
K_equal_One <- K_equal_One[complete.cases(K_equal_One)]
K_equal_One_Relaxed <- data.frame(K_equal_One)


#### This section calculates all the K's for Omega class 2 for significantly positive selected genes ###
### therefore read in Codeml summary table ###
Input1 <- read.table("/beegfs/group_dv/home/LIasi/codeml_output/sum_Nothobranchius.txt", header=F, fill = T, stringsAsFactors = F, sep = "\t", quote="")
Input1 <- Input1[complete.cases(Input1),]
### filter for significant results ###
Input1_sign <- subset(Input1, V5 <0.05)

### read in RELAX output again (not really neccessary)###
Input2 <- read.table("/beegfs/group_dv/home/LIasi/relax_output/sum_Nothobranchius.txt", header=F, fill = T, stringsAsFactors = F, sep = "\t", quote="")
datIn2 <- Input2[complete.cases(Input2),]

### subset RELAX output by the Codeml signifiv=cant results ###
Codeml_sign2 <- subset(datIn2, V1 %in% Input1_sign$V1)

#### calculate all the K's for Omega class 2 for significantly positive selected genes ###
K_bigger_One <-log(Codeml_sign2$V21,Codeml_sign2$V15)
K_bigger_One <- K_bigger_One[complete.cases(K_bigger_One)]
K_bigger_One <- K_bigger_One[K_bigger_One>0]
K_bigger_One_positiv_selected <- data.frame(K_bigger_One)


### write the results for the K's in tables from which the Get_Omega_K.R programm can read from ###
write.table(file="K_smaller_One_Relaxed.txt",K_smaller_One_Relaxed)
write.table(file="K_equal_One_Relaxed.txt",K_equal_One_Relaxed)
write.table(file = "K_bigger_One_positiv_selected.txt",K_bigger_One_positiv_selected)

#### This section calculates all the K's for Omega class 0 and 1 for significantly intemsified  genes ###
### filter RELAX output for K bigger 1 ###
datIn_sign_Intensified <- subset(datIn, V9 >1) 

#### calculate all the K's for Omega class 0 for significantly intensified genes ###
K_smaller_One <- log(datIn_sign_Intensified$V17,datIn_sign_Intensified$V11)
K_smaller_One <- K_smaller_One[complete.cases(K_smaller_One)]
#get rid of inf values
K_smaller_One_Intensified <- K_smaller_One[is.finite(K_smaller_One)]
K_smaller_One_Intensified <- K_smaller_One_Intensified[K_smaller_One_Intensified>1]

K_smaller_One_Intensified <- data.frame(K_smaller_One_Intensified)


#### calculate all the K's for Omega class 1 for significantly intensified genes ###
K_equal_One <-log(datIn_sign_Intensified$V19,datIn_sign_Intensified$V13)
K_equal_One <- K_equal_One[complete.cases(K_equal_One)]
#get rid of inf values of K
K_equal_One_Intensified <- K_equal_One[is.finite(K_equal_One)]
K_equal_One_Intensified <- K_equal_One_Intensified[K_equal_One_Intensified>1]

K_equal_One_Intensified <- data.frame(K_equal_One_Intensified)

### write the results for the K's in tables from which the Get_Omega_K.R programm can read from ###
write.table(file="K_smaller_One_Intensified.txt",K_smaller_One_Intensified)
write.table(file="K_equal_One_Intensified.txt",K_equal_One_Intensified)


#### In this section the results can be plotted with an abline marking the upper (darkgreen) and lower (darkred) boundrys of the used K values ###
pdf("K smaller 1 for Relaxation.pdf")
hist(K_smaller_One, breaks = 100,xlab = "K ",xlim = c(0,5), main = "Distribution of K for Class smaller 1 Relaxed")
abline(v=1,col="darkred",lwd= 3)
dev.off()
pdf("K equal 1 for Relaxation.pdf")
hist(K_equal_One, breaks = 100,xlab = "K ",xlim = c(0,5), main = "Distribution of K for Class equal 1 Relaxed")
abline(v=1,col="darkred",lwd= 3)
dev.off()
pdf("K bigger 1 for Positive Selection.pdf")
hist(K_bigger_One, breaks = 1000,xlab = "K ", xlim = c(0,5),main = "Distribution of K for Class bigger 1 Positiv Selection")
abline(v=1,col="darkgreen",lwd= 3)
dev.off()

### same plot but now saved in 1 pdf document ###
pdf("K distribution used for the Simulations.pdf")
split.screen(c(3,1) )
screen(1)
hist(K_smaller_One, breaks = 100,xlab = "K ",xlim = c(0,5), main = "Distribution of K for Class smaller 1 Relaxed")
abline(v=1,col="darkred",lwd= 3)
screen(2)
hist(K_equal_One, breaks = 100,xlab = "K ",xlim = c(0,5), main = "Distribution of K for Class equal 1 Relaxed")
abline(v=1,col="darkred",lwd= 3)
screen(3)
hist(K_bigger_One, breaks = 1000,xlab = "K ", xlim = c(0,5),main = "Distribution of K for Class bigger 1 Positiv Selection")
abline(v=1,col="darkgreen",lwd= 3)
close.screen(all=T)
dev.off()
