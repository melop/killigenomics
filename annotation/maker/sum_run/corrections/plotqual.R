#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/NOR_final_1/run_pilon_exonerateOn/corrections");

args = commandArgs(trailingOnly=TRUE);

sCountTab <-args[1];  #"quality.refined.removeisoform.txt";  "frag2aln.quality.txt"; #

dat <- read.table(sCountTab , header=T);

arrBreaks <- seq(0, 1, 0.01);
arrLegendLabs <- c();
arrLegendCols <- c();

fnPlotAcc <- function(v, nStart=0, nEnd=1, nStep=0.01, ylow=0, yhigh=length(v), sCol="blue", sName="Completeness", ifadd=F) {
  dat1Cut <- cut(v, arrBreaks, right=TRUE) 
  dat1Cumfreq <- cumsum(table(dat1Cut));
  dat1RelCumfreq <- dat1Cumfreq / length(v);
  if (ifadd) {
    lines(arrBreaks[1:(length(arrBreaks)-1) ], dat1Cumfreq, xlim=c(nStart,nEnd), lwd=4, type = "l",col = sCol, xlab=sName, ylab="count", ylim=c(ylow, yhigh));
    
  } else {
    plot(arrBreaks[1:(length(arrBreaks)-1) ], dat1Cumfreq, xlim=c(nStart,nEnd),  lwd=4, type = "l",col = sCol, xlab=sName, ylab="count", ylim=c(ylow, yhigh));
    arrLegendLabs <<- c();
    arrLegendCols <<- c();
  }

  arrLegendLabs <<- c(arrLegendLabs , sName);
  arrLegendCols <<- c(arrLegendCols , sCol);
}

fnPlotLegend <- function() {
	legend('topleft'  , legend=arrLegendLabs , col=arrLegendCols, lwd=4);
}


pdf(file=paste(sCountTab,".pdf",sep=""))

fnPlotAcc(dat$completeness, sName="% complete")
fnPlotLegend();


fnPlotAcc(dat$X5primemissperc, sName="5 prime miss %", ifadd=F , sCol="orange")
fnPlotAcc(dat$X3primemissperc, sName="3 prime miss %" , ifadd=T, sCol="red")
fnPlotLegend();

fnPlotAcc(dat$X5primeuncertain / dat$score , sName="5 prime uncertain %" , ifadd=F, sCol="green")
fnPlotAcc(dat$X3primeuncertain / dat$score , sName="3 prime uncertain %" , ifadd=T, sCol="darkgreen")
fnPlotLegend();

fnPlotAcc(dat$smallinsertionlen / dat$score , sName="< 3AA insertion %" , ifadd=F, sCol="purple")
fnPlotAcc(dat$biginsertionlen / dat$score , sName="3AA - 10AA insertion %" , ifadd=T, sCol="violet")
fnPlotAcc(dat$hugeinsertionlen / dat$score , sName=">10AA insertion %" , ifadd=T, sCol="black")
fnPlotLegend();

fnPlotAcc(dat$smalldeletionlen / dat$score , sName="< 3AA deletion %" , ifadd=F, sCol="purple")
fnPlotAcc(dat$bigdeletionlen / dat$score , sName="3AA - 10AA deletion %" , ifadd=T, sCol="violet")
fnPlotAcc(dat$hugedeletionlen / dat$score , sName=">10AA deletion %" , ifadd=T, sCol="black")
fnPlotLegend();


fnPlotAcc(dat$hicov / dat$score , sName="abnormal coverage %" , ifadd=F, sCol="pink")
fnPlotLegend();


#plot(dat$X5primemissperc, dat$X3primemissperc);


dev.off();

sink(file=paste(sCountTab,".report.txt",sep=""));
arrCutoffs <- c(0.5, 0.7, 0.8, 0.9, 0.95);
n90 <- 0;

cat("Completeness statistics: \n");
cat("no filter: " , length(dat$completeness), "\n");
for(nCutoff in arrCutoffs) {
  nPassingCount <- length(dat$completeness[dat$completeness>=nCutoff]);
  cat(" >= " , nCutoff, " : " , nPassingCount , "\n");
  if (nCutoff == 0.9) {
    n90 <- nPassingCount;
    
  }
}

cat("<3AA insertion sum bp: " , sum(dat$smallinsertionlen, na.rm = T) , " count : " , sum(dat$smallinsertioncount, na.rm = T) , " per gene avg AA: " , sum(dat$smallinsertionlen, na.rm = T) / n90 , "\n" );
cat("<3AA deletion sum bp: " , sum(dat$smalldeletionlen, na.rm = T) , " count : " , sum(dat$smalldeletioncount, na.rm = T) , " per gene avg AA: " , sum(dat$smalldeletionlen, na.rm = T) / n90 ,"\n" );

cat("3-10AA insertion sum bp: " , sum(dat$biginsertionlen, na.rm = T) , " count : " , sum(dat$biginsertioncount, na.rm = T) , " per gene avg AA: " , sum(dat$biginsertionlen, na.rm = T) / n90 ,"\n" );
cat("3-10AA deletion sum bp: " , sum(dat$bigdeletionlen, na.rm = T) , " count : " , sum(dat$bigdeletioncount, na.rm = T) ,  " per gene avg AA: " , sum(dat$bigdeletionlen, na.rm = T) / n90 ,"\n" );

cat(">10AA insertion sum bp: " , sum(dat$hugeinsertionlen, na.rm = T) , " count : " , sum(dat$hugeinsertioncount, na.rm = T) ,  " per gene avg AA: " , sum(dat$hugeinsertionlen, na.rm = T) / n90 ,"\n" );
cat(">10AA deletion sum bp: " , sum(dat$hugedeletionlen, na.rm = T) , " count : " , sum(dat$hugedeletioncount, na.rm = T) ,  " per gene avg AA: " , sum(dat$hugedeletionlen, na.rm = T) / n90 ,"\n" );

sink()

