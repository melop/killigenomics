### Tool to draw random numbers from Omega Class distributions ###

    K_smaller_One<- read.table(file="/beegfs/group_dv/home/LIasi/Simulation_Tools/K_smaller_One_Relaxed.txt")
    K_smaller_One <- K_smaller_One[K_smaller_One$K_smaller_One<=1,]

    K_equal_One<-read.table(file="/beegfs/group_dv/home/LIasi/Simulation_Tools/K_equal_One_Relaxed.txt")
    K_equal_One <- K_equal_One[K_equal_One$K_equal_One<=1,]

    K_bigger_One <- read.table(file="/beegfs/group_dv/home/LIasi/Simulation_Tools/K_bigger_One_positiv_selected.txt")
    K_bigger_One <- K_bigger_One[K_bigger_One$K_bigger_One>1,]

    K_smaller_One_Intensified <- read.table(file="/beegfs/group_dv/home/LIasi/Simulation_Tools/K_smaller_One_Intensified.txt")
    K_smaller_One_Intensified <- K_smaller_One_Intensified[K_smaller_One_Intensified$K_smaller_One_Intensified>1,]

    K_equal_One_Intensified <- read.table(file="/beegfs/group_dv/home/LIasi/Simulation_Tools/K_equal_One_Intensified.txt")
    K_equal_One_Intensified <- K_equal_One_Intensified[K_equal_One_Intensified$K_equal_One_Intensified>1,]

fnProcessFromPHPInput <- function (Input1) {
  Input <- read.table(text=Input1, header=F, fill = T, strip.white = F, stringsAsFactors = F, na.strings = "", sep = "\t")
  if (is.null(Input)) {
    cat("ERROE: No Input")
  } 
  

  if (Input$V1==1) {

   output <-  sample(K_smaller_One,1)
    
  }
  
  if (Input$V1==2) {
    

    output <- sample(K_equal_One,1)
    
  }
  
  if (Input$V1==3) {

    output <-  sample(K_bigger_One,1)
  }
 
  if (Input$V1==4) {
    
    output <-  sample(K_smaller_One_Intensified,1)
  }
  
  if (Input$V1==5) {
    
    output <-  sample(K_equal_One_Intensified,1)
  }
  
  if (Input$V1==6) {
    
   
    output <-  1
  }
 
  cat(output,"\n")
} 

#  else { output <-10^(rnorm(1, mean = m, sd = std))} 
#write((rnorm(1, mean = m, sd = std) ^ 10),file="")


fStdin <- file("stdin")

open(fStdin)


while(length(sLn <- readLines(fStdin,n=1, warn = T)) > 0) {
  fnProcessFromPHPInput(sLn)
}

