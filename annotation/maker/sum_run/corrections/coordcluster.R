library(cluster)
#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/PLP_final_1/run_pilon_exonerateOn/corrections");
fnOutliers <- function(arr) {
  nMed <- median(arr);
  nSd <- sd(arr);
  nLower <- nMed - nSd;
  nUpper <- nMed + nSd;
  return(arr[arr<nLower | arr>nUpper ]);
}
args = commandArgs(trailingOnly=TRUE);

#dat <- read.table( "testcoord.txt", header=T );
dat <- read.table( args[1], header=T );
datUnique <- dat[!duplicated(dat[, c(1,2)]) , ];
k.max <- min(15 , nrow(datUnique));

data <- scale(dat[, -3]);
data[is.nan(data)] <- 0;
sil <- rep(0, k.max)
nIgnoreOverlap <- 10;


# Compute the average silhouette width for 
# k = 2 to k = 15
if (k.max >=2) {
  for(i in 2:k.max){
    tryCatch( {
    km.res <- kmeans(data, centers = i, nstart = 25)
    ss <- silhouette(km.res$cluster, dist(data))
    sil[i] <- mean(ss[, 3])
    },
    error = function(x) {}, warnings = function(x) {}
    )
  }

  # Plot the  average silhouette width
  #plot(1:k.max, sil, type = "b", pch = 19, 
  #     frame = FALSE, xlab = "Number of clusters k")
  #abline(v = which.max(sil), lty = 2)
  nBestCluster <- which.max(sil);
  km.res <- kmeans(data, centers = nBestCluster, nstart = 25);
  
  datOut <- cbind(dat, Cluster=km.res$cluster);
  #filter out outlier clusters
} else {
  datOut <- cbind(dat, Cluster=1);
  
}

datMedian <- aggregate(datOut[,c(1,2,4)], list(clu=datOut$Cluster), median);
datCount <- aggregate(datOut[,'Cluster'], list(clu=datOut$Cluster), length);

#now compare the clusters
arrExclude = c();

for(nClu1 in datMedian$clu) {
  if (nClu1 %in% arrExclude) next;
  
  datClu1 <- datMedian[datMedian$clu==nClu1, ];
  nStartClu1 <- datClu1$Start;
  nEndClu1 <- datClu1$End;
  nCountClu1 <- datCount[datCount$clu==nClu1, 'x'];
  
  for(nClu2 in datMedian$clu) {
    if (nClu1 >= nClu2) next;
    if (nClu2 %in% arrExclude || nClu1 %in% arrExclude) next;
    datClu2 <- datMedian[datMedian$clu==nClu2, ];
    nStartClu2 <- datClu2$Start;
    nEndClu2 <- datClu2$End;
    nCountClu2 <- datCount[datCount$clu==nClu2, 'x'];
    #check if overlapping:
    if ( (nStartClu1 + nIgnoreOverlap ) < nEndClu2 && (nStartClu2 + nIgnoreOverlap) < nEndClu1) { #found overlap
      if (nCountClu1 > nCountClu2) { #exclude the one with fewer counts
          arrExclude <- c(arrExclude , nClu2);
      } else if (nCountClu1 < nCountClu2) {
          arrExclude <- c(arrExclude , nClu1);
      } else { #when it is a tigh, take the lower ending position.
        if (nEndClu1 < nEndClu2) {
          arrExclude <- c(arrExclude , nClu2);
        } else {
          arrExclude <- c(arrExclude , nClu1);
        }
      }
    }
    
  }
}

datOut <- cbind(datOut , Keep=1);
for(sEx in arrExclude) {
  datOut[datOut$Cluster == sEx , 'Keep'] <- 0;
}

#now go over each cluster, and exclude the outliers
arrKeptClusters <- setdiff(datOut$Cluster, arrExclude);
for( nClu in arrKeptClusters) {
  datClu <- datOut[datOut$Cluster == nClu,];
  if (nrow(datClu) < 3) { #if too few alignments, skip
    next;
  }
  
  arrStarts <- datClu$Start;
  arrEnds <- datClu$End;
  arrOutlierStart <- fnOutliers(arrStarts);
  arrOutlierEnd <- fnOutliers(arrEnds);
  
  datOut[datOut$Cluster == nClu & (datOut$Start %in% arrOutlierStart) , 'Keep' ] <-0;
  datOut[datOut$Cluster == nClu & (datOut$End %in% arrOutlierEnd) , 'Keep' ] <-0;
  
}

datOutFilter <- datOut[datOut$Keep==1 ,];
datMedian <- aggregate(datOutFilter[ ,c(1,2,4)], list(clu=datOutFilter$Cluster), median);


write.table(datOut, file= args[2] , col.names = T, quote = F, row.names = F, sep="\t");
write.table(datMedian, file= args[3] , col.names = T, quote = F, row.names = F, sep="\t");

