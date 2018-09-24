f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
args = strsplit(line, '\t') [[1]]
if (length(args)!=2) next;
if (is.nan(as.numeric(args[1])) || is.nan(as.numeric(args[2])) ) next;
nProb <- ppois(as.numeric(args[1]) ,as.numeric( args[2]) );
if (nProb > 0.50) {
	nProb <- 1 - nProb;
}
cat(nProb*2,"\n");

}
