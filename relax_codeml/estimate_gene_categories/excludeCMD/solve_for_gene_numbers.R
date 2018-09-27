setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/sim_nothos_pos");

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
sGenus <- "Nothobranchius";
datReal <- read.table(paste("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/sim_nothos_pos/real_ret_",sGenus, ".txt", sep=""), header=T);
nExcludedGenes <- datReal[1,'relaxexl'] ; #these are k>1 $11==0, omega2> omega1, simulation shows that it is likely that these are relaxed genes misclassified as intensified genes (pos sel)
                        #filter by cat ../../relax_46_sim/sum_Nothobranchius.txt |  awk '{if ( $11==0 && $15<$21 && $5 <0.05 && $9>1) print $0; }'
nExcludedGenes2<- datReal[1,'intenseexl'] ;#these are k>1 $11==0, omega2< omega1, simulation shows that it is likely that these are intensified genes misclassified as relaxed genes (pos sel)
#filter by cat ../../relax_46_sim/sum_Nothobranchius.txt |  awk '{if ( $11==0 && $15>$21 && $5 <0.05 && $9<1) print $0; }'

#total genes                         relax              intense                       pos                                             overlap          
b <- c(datReal[1,'relaxtotal'] -nExcludedGenes-nExcludedGenes2, datReal[1,'relaxsig'] ,datReal[1,'intensesig'] , datReal[1,'codemlsig'],  datReal[1,'overlapposrelax'], datReal[1,'overlapposintense'] ) #Nothobranchius
# 
#  sGenus <- "Aphyosemion";
#  
# nExcludedGenes <- 493 ; #
# nExcludedGenes2<- 179;
#    b <- c(13523-nExcludedGenes-nExcludedGenes2,1577-nExcludedGenes2,1044-nExcludedGenes, 1290, 192) #Aphyosemion

# sGenus <- "Callopanchax";
#  nExcludedGenes <- 316 ; #Callopanchax
#  b <- c(13545-nExcludedGenes,1097,785-nExcludedGenes, 307, 41) #Callopanchax ; count overlap with comm -12 callopanchax_pos.txt ../relax_46_spp/callopanchax_relaxed.txt | wc -l 

# sGenus <- "Scriptaphyosemion";
# 
#  nExcludedGenes <- 247 ; #Scriptaphyosemion
#  b <- c(13484-nExcludedGenes,546,1015-nExcludedGenes, 565, 24) #Scriptaphyosemion ; count overlap with comm -12 callopanchax_pos.txt ../relax_46_spp/callopanchax_relaxed.txt | wc -l 
 
 
arrRawEst <- c(b[4], b[2], b[3], b[1]-sum(b[2:5]), b[5] ) / b[1];

A <- as.matrix(read.table(paste("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/sim_nothos_pos/sum_ret_",sGenus, ".txt", sep=""), header=T, row.names = 1))[,1:5];
A <- t(cbind(c(1,1,1,1,1,1), A));

#b <- c(13708,2729,1008, 1464, 450)


arrEstimates <- solve(A,b) / b[1];
names(arrEstimates) <- c("Pos", "Relaxed", "Intensified_purify", "NoChange", "Pos_and_Relaxed");

nRep <- 10000; #100 bootstraps
arrSampledEst <- NULL;

arrProbs <- c( (b[1]-sum(b[2:5])) / b[1] ,  b[2]/b[1], b[3]/b[1], b[4]/b[1], b[5]/b[1]);
for(i in 1:nRep) {
  arrDrawn <- rmultinom(n = 1, size = b[1], prob = arrProbs);
  bDrawn <- c( sum(arrDrawn), arrDrawn[2:5] );
  oEst <- solve(A,bDrawn)/bDrawn[1];
  
  arrSampledEst <- rbind(arrSampledEst, oEst);
}

colnames(arrSampledEst) <- names(arrEstimates) ;
#barplot(arrEstimates, ylim=c(0,1));

arrCols <- c(rgb(1,0,0,1/2), rgb(0,0,1,1/2), rgb(0,1,0,1/2), rgb(1/3,1/3,1/3,1/2), rgb(1,0,1,1/2) );
bFirst <- T;
for(i in 1:ncol(arrSampledEst)) {
    hist(arrSampledEst[,i], breaks=100, col=arrCols[i], add=(!bFirst), xlim=c(-0.05,1),lty="blank", freq = T, xlab="Estimated Proportion of Genes", main=sGenus);
    points(x=arrEstimates[i], y=20*nRep/1000, pch=16, cex=1)
  points(x=arrRawEst[i], y=15*nRep/1000, pch=18 , col=arrCols[i], cex=2)
      bFirst <- F;
}
# 
# #solve linear equation to calculate Positive, Relax, Intensified, And unchanged genes
# 
# A <- matrix(c(1,0.0509, 0.2352, 0.3193, 1, 0.5107, 0.1471, 0.0268, 1, 0.1471, 0.5, 0.0268, 1, 0.0258, 0.0335, 0.0118), 4, 4)
# b <- c(13708,2729,1008, 1464)
# 
# solve(A,b)
solve(A,b)
arrEstimates
