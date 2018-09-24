dat <- read.table("./mapped/BUSCO.cov", header=F);
#hist(dat$V5)
cat(mean(dat$V5))
