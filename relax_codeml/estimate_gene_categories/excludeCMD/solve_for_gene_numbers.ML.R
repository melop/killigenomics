setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/estimate_gene_categories/excludeMNM");
library(plotrix)
# simulated relax_k<1	   relax_k>1   	codeml	  intersect(relax_k<1, codeml)
# pos         a	            f           	k	         p
# relax	      b	            g	            l	         q
# intensified	c	            h	            m          r
# neutral	    d	            i	            n	         s
# relax_pos	  e	            j	            o	         t

# Positive only = P, Relax_only=R, Intensified_purification = I, NoChange = N, Relaxed_and_positive = B

# linear equations:

# 1. P + R + I + N + B = 13708
# 2. aP + bR + cI + dN + eB = 2729   relax test, p<0.05, k<1
# 3. fP + gR + hI + iN + jB = 1008   relax test, p<0.05, k>1
# 4. kP + lR + mI + nN + oB = 1464   codeml test, p < 0.05
# 5. pP + qR + rI + sN + tB = 450    intersect(relax test, p<0.05, k<1 , codeml test, p < 0.05) 

#solve linear equation to calculate Positive, Relax, Intensified, And unchanged genes

#A <- matrix(c(1,0.5107, 0.2352+0.1471, 0.3193,0.3193*0.5107, 1,0.0509, 0.2352, 0.3193, 0,   1, 0.5107, 0.1471, 0.0268, 0, 1, 0.1471, 0.5107, 0.0268, 0, 1, 0.0258, 0.0335, 0.0118,0), 5, 5)
sGenus <- "Nothobranchius"; #real data
#sGenus <- "Aphyosemion"; #real data
#sGenus <- "Callopanchax"; #real data
#sGenus <- "Scriptaphyosemion"; #real data


sGenus2 <- "Nothobranchius"; #sim matrix
#sGenus2 <- "Aphyosemion"; #real data
#sGenus2 <- "Callopanchax"; #real data
#sGenus2 <- "Scriptaphyosemion"; #real data


datReal <- read.table(paste("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/estimate_gene_categories/excludeMNM/real_ret_",sGenus, ".txt", sep=""), header=T);
nExcludedGenes <- datReal[1,'relaxexl'] ; #these are k>1 $11==0, omega2> omega1, simulation shows that it is likely that these are relaxed genes misclassified as intensified genes (pos sel)
                        #filter by cat ../../relax_46_sim/sum_Nothobranchius.txt |  awk '{if ( $11==0 && $15<$21 && $5 <0.05 && $9>1) print $0; }'
nExcludedGenes2<- datReal[1,'intenseexl'] ;#these are k>1 $11==0, omega2< omega1, simulation shows that it is likely that these are intensified genes misclassified as relaxed genes (pos sel)
#filter by cat ../../relax_46_sim/sum_Nothobranchius.txt |  awk '{if ( $11==0 && $15>$21 && $5 <0.05 && $9<1) print $0; }'


arrRawEst <- c(b[4], b[2], b[3], b[1]-sum(b[2:5]), b[5] ) / b[1];

A <- as.matrix(read.table(paste("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/estimate_gene_categories/excludeMNM/sum_ret_",sGenus2, ".txt", sep=""), header=T)[,2:6]);
A <- (cbind(c(1,1,1,1,1,1), A));
#total genes                         relax              intense                       pos                                             overlap          
b <- c(datReal[1,'relaxtotal'] -nExcludedGenes-nExcludedGenes2, datReal[1,'relaxsig'] ,datReal[1,'intensesig'] , datReal[1,'codemlsig'],  datReal[1,'overlapposrelax'], datReal[1,'overlapposintense'] ) #Nothobranchius

# A <- as.matrix(read.table(paste("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/sim_nothos_pos/sum_ret_",sGenus, ".txt", sep=""), header=T)[1:4,2:4]);
# A <- (cbind(c(1,1,1,1), A));
# 
# #total genes                         relax              intense                       pos                                             overlap          
# b <- c(datReal[1,'relaxtotal'] -nExcludedGenes-nExcludedGenes2, datReal[1,'relaxsig'] ,datReal[1,'intensesig'] , datReal[1,'codemlsig'] ) #Nothobranchius
# 

fnObjectiveFunc <- function(params, bDr) {
  #params are the real values, in the order of positive genes, relaxed, intensified pur, overlapped. neutral genes = total - sum(others)
  #              pos        relax      intense       neutral          overlap_pos_rel  overlap_pos_int
  #paramEst <- c(params[1], params[2], params[3], bDr[1] - sum(params) , params[4] , params[5]);
  paramEst <- params;
  #paramEst <- c(params[1], params[2], params[3], params[4]);
  arrEst <- paramEst %*% A;
  #compute eucliedian distance
  #return(sqrt(sum(( arrEst - bDr )^2 )));
  #return( sum( ( abs(arrEst - bDr )) ) );
  
  return(sum(( bDr - arrEst  )^2 / arrEst )); #chi-square stat
  #return(sum(( bDr - arrEst  ) / arrEst ));
}

fnSolve <- function(bDr) {
  init_val <- bDr[1]/6;
  lowerbound <- 0;
  upperbound <- bDr[1] ;
  #arrMLEstimates <- optim(par=c( init_val, init_val, init_val, init_val, init_val) , fn=fnObjectiveFunc , NULL, method = "L-BFGS-B" , lower=c( lowerbound , lowerbound, lowerbound,lowerbound,lowerbound) , upper=c( upperbound, upperbound, upperbound , upperbound, upperbound) , control=list(fnscale=1, ndeps=c( 1e-8 , 1e-8, 1e-8, 1e-8, 1e-8), parscale=c(0.1, 0.1,0.1,0.1,0.1), maxit=1000000, factr=1e8), bDr=bDr)
  arrMLEstimates <- optim(par=c( init_val, init_val, init_val, init_val, init_val, init_val) , fn=fnObjectiveFunc , NULL, method = "L-BFGS-B" , lower=c( lowerbound , lowerbound, lowerbound,lowerbound,lowerbound,lowerbound) , upper=c( upperbound, upperbound, upperbound , upperbound, upperbound, upperbound) , control=list(fnscale=1, ndeps=c( 1e-8 , 1e-8, 1e-8, 1e-8, 1e-8, 1e-8), parscale=c(0.1, 0.1,0.1,0.1,0.1,0.1), maxit=1000000, factr=1e8), bDr=bDr)
  return(arrMLEstimates$par)
  #return(c(arrMLEstimates$par[1:3], bDr[1]-sum(arrMLEstimates$par), arrMLEstimates$par[4], arrMLEstimates$par[5] ) );

  # init_val <- bDr[1]/4;
  # lowerbound <- 0;
  # upperbound <- bDr[1] ;
  # arrMLEstimates <- optim(par=c( init_val, init_val, init_val, init_val) , fn=fnObjectiveFunc , NULL, method = "L-BFGS-B" , lower=c( lowerbound , lowerbound, lowerbound, lowerbound) , upper=c( upperbound, upperbound, upperbound, upperbound ) , control=list(fnscale=1, ndeps=c( 1e-8 , 1e-8, 1e-8, 1e-8), parscale=c(0.1, 0.1,0.1,0.1), maxit=1000000, factr=1e8), bDr=bDr)
  # #return(c(arrMLEstimates$par[1:3], bDr[1]-sum(arrMLEstimates$par) ) );
  # 
  # return(arrMLEstimates);
}


nRep <- 1000; #100 bootstraps
arrSampledEst <- NULL;

arrProbs <- c( (b[1]-sum(b[2:6])) / b[1] ,  b[2]/b[1], b[3]/b[1], b[4]/b[1], b[5]/b[1] , b[6]/b[1]);
for(i in 1:nRep) {
  arrDrawn <- rmultinom(n = 1, size = b[1], prob = arrProbs);
  bDrawn <- c( sum(arrDrawn), arrDrawn[2:6] );
  oEst <- fnSolve(bDrawn)
  oEst <- oEst/sum(oEst);
  
  arrSampledEst <- rbind(arrSampledEst, oEst);
}
colnames(arrSampledEst) <- c("Pos", "Relaxed", "Intensified_purify", "NoChange", "Pos_and_Relaxed", "Pos_and_Intense");

#barplot(arrEstimates, ylim=c(0,1));

arrCols <- c(rgb(1,0,0,1/2), rgb(0,0,1,1/2), rgb(0,1,0,1/2), rgb(1/3,1/3,1/3,1/2), rgb(1,0,1,1/2), rgb(0,1,1,1/2) );
bFirst <- T;
# for(i in 1:ncol(arrSampledEst)) {
#     hist(arrSampledEst[,i], breaks=50, col=arrCols[i], add=(!bFirst), xlim=c(-0.05,1),lty="blank", freq = T, xlab="Estimated Proportion of Genes", main=sGenus);
#     points(x=arrEstimates[i], y=20*nRep/1000, pch=16, cex=1)
#  # points(x=arrRawEst[i], y=15*nRep/1000, pch=18 , col=arrCols[i], cex=2)
#       bFirst <- F;
# }
# 
# #solve linear equation to calculate Positive, Relax, Intensified, And unchanged genes
# 
# A <- matrix(c(1,0.0509, 0.2352, 0.3193, 1, 0.5107, 0.1471, 0.0268, 1, 0.1471, 0.5, 0.0268, 1, 0.0258, 0.0335, 0.0118), 4, 4)
# b <- c(13708,2729,1008, 1464)
# 
# solve(A,b)
arrEstimates <- fnSolve(b);
arrEstimates
arrEstimates %*% A;

arrEstimates <- arrEstimates/ sum(arrEstimates);

names(arrEstimates) <- colnames(arrSampledEst) ;

colMeans(arrSampledEst) * b[1]
colMeans(arrSampledEst)

#boxplot(arrSampledEst);
pdf(file=paste("estimated_gene_categories_", sGenus, ".pdf", sep=""), width=5, height=5);
violin_plot(arrSampledEst, x_axis_labels = colnames(arrSampledEst) , col=arrCols , show_mean = T, main=sGenus, ylim=c(0,1))
dev.off();

arrEstimates <- c(arrEstimates , sum(arrEstimates[c(1,5,6)]), sum(arrEstimates[c(2,5)]), sum(arrEstimates[c(3,6)])  );
names(arrEstimates)[7:9] <- c('Pos_all',	'Relax_all',	'Intense_all');
write.table(t(arrEstimates) , file=paste("estimated_gene_categories_", sGenus, ".txt", sep=""), col.names = T, row.names = F, quote = F, sep="\t")

arrSampledEst <- cbind(arrSampledEst , rowSums(arrSampledEst[, c(1,5,6)]) , rowSums(arrSampledEst[, c(2,5)]), rowSums(arrSampledEst[, c(3,6)])  );
colnames(arrSampledEst)[7:9] <- c('Pos_all',	'Relax_all',	'Intense_all');

write.table(arrSampledEst , file=paste("estimated_gene_categories_resampled_", sGenus, ".txt", sep=""), col.names = T, row.names = F, quote = F, sep="\t")
