#
#
#		asymptoticMK: Asymptotic McDonald-Kreitman Test web service
#
#		By Benjamin C. Haller and Philipp W. Messer
#		Copyright (C) 2017 Philipp Messer.
#
#		This web service should be available at http://benhaller.com/messerlab/asymptoticMK.html
#		The Github repository for asymptoticMK is at https://github.com/MesserLab/asymptoticMK
#
#

#  This file is part of asymptoticMK.
#
#  asymptoticMK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#  asymptoticMK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with asymptoticMK.  If not, see <http://www.gnu.org/licenses/>.


#######################################################################################
#
#	To run asymptoticMK locally, first source all the functions between this comment
#	and the next comment, below, that is in the same box style.  Then, read on from
#	that second comment.  Everything between the two box comments may be considered a
#	black box; you should not need to understand it or modify it to use asymptoticMK.
#
#######################################################################################

# Get a CI using Monte Carlo simulation based upon a fitted model.  This is necessary because
# getting confidence intervals for non-linear models is a complicated business, apparently.
# Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
# See: https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# Or, if that link goes stale: http://stats.stackexchange.com/a/251501/141766
predictNLS <- function(
  object, 
  newdata,
  level = 0.95, 
  nsim = 10000,
  ...
)
{
  require(MASS, quietly = TRUE)
  
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## coefficients
  COEF <- coef(object)
  
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  
  ## check that 'newdata' has same name as predVAR
  if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
  
  ## get parameter coefficients
  COEF <- coef(object)
  
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
  
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
  
  ## define counter function
  counter <- function (i) 
  {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  
  outMAT <- NULL 
  
  for (i in 1:NR) {
    #counter(i)		# show a counter for lengthy fits; commented out to reduce noise here...
    
    ## get predictor values and optional errors
    predVAL <- newdata[i, 1]
    if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
    
    ## create mean vector for 'mvrnorm'
    MU <- c(COEF, predVAL)
    
    ## create variance-covariance matrix for 'mvrnorm'
    ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
    
    ## create MC simulation matrix
    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
    
    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
    
    ## collect statistics
    PRED <- data.frame(predVAL)
    colnames(PRED) <- predNAME   
    FITTED <- predict(object, newdata = data.frame(PRED))
    MEAN.sim <- mean(EVAL, na.rm = TRUE)
    SD.sim <- sd(EVAL, na.rm = TRUE)
    MEDIAN.sim <- median(EVAL, na.rm = TRUE)
    MAD.sim <- mad(EVAL, na.rm = TRUE)
    QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
    RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
    outMAT <- rbind(outMAT, RES)
  }
  
  colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
  rownames(outMAT) <- NULL
  
  #cat("\n")	# commented out along with the call to counter() above
  
  return(outMAT)  
}

# core code for two-step nls2() model fit at a given level of precision (res)
fitMKmodel <- function(alpha_trimmed, f_trimmed, res)
{
  require(nls2)
  
  mod <- tryCatch({
    st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
  },
  error=function(cond) {})
  
  if (length(mod) == 0)
    return(NULL)
  
  mod2 <- tryCatch({
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
  },
  error=function(cond) {})
  
  if (length(mod2) == 0)
    return(NULL)
  
  return(mod2)
}

# You can customize this function to support additional output devices.  This is probably
# not a very good design; I'm not sure what a better design would be, however, since I've
# never really understood R's device model very well...
createOutputDevice <- function(output, plotNum, width=3.5, height=3.5, bg="white", ...)
{
  titleString <- paste0("Figure ", plotNum)
  filenameString <- paste0(titleString, ".pdf")
  if (output == 'default') {
    return();
  }
  
  if (output == "quartz")
  {
    quartz(title=titleString, width=width, height=height, bg=bg, ...)
  }
  else if (output == "pdf")
  {
    pdf(file=filenameString, width=width, height=height, bg=bg, ...)
  }
  else
  {
    stop(paste0("Unsupported output device name: ", output))
  }
}

closeOutputDevice <- function(output, plotNum)
{
  if (output == 'default') {
    return();
  }
  if (output == "quartz")
  {
    # leave quartz windows open, since they are interactive
  }
  else if (output == "pdf")
  {
    # close PDF devices, since they connect to the filesystem
    dev.off()
  }
  else
  {
    stop(paste0("Unsupported output device name: ", output))
  }
}


#
#
#	asymptoticMK() function for use in R on a local machine
#
#
asymptoticMK <- function(d0, d, xlow, xhigh, df, true_alpha=NA, output="table", arrYLims=NULL,  ...)
{
  require(nls2)
  require(MASS)
  
  if (is.na(d0) || is.null(d0))
    stop("Malformed d0 (must be numeric).")
  if (is.na(d) || is.null(d))
    stop("Malformed d (must be numeric).")
  if (is.na(xlow) || is.null(xlow))
    stop("Malformed xlow (must be numeric).")
  if (is.na(xhigh) || is.null(xhigh))
    stop("Malformed xhigh (must be numeric).")
  
  #	Bounds-check response variables
  #
  if (d0 <= 0)
    stop("d0 must greater than zero.")
  if (d <= 0)
    stop("d must greater than zero.")
  if ((xlow < 0.0) || (xlow > 1.0))
    stop("xlow must be in the interval [0,1].")
  if ((xhigh < 0.0) || (xhigh > 1.0))
    stop("xhigh must be in the interval [0,1].")
  if (xlow >= xhigh)
    stop("xlow must be less than xhigh.")
  
  #	Read in the file and check its format
  #
  if (NCOL(df) != 3)
    stop("Dataframe df does not contain exactly three tab-separated columns.")
  if (NROW(df) <= 0)
    stop("Dataframe df contains no data rows.")
  
  cols <- names(df)
  
  suppressWarnings(	# the goal is to generate NAs here, so we don't want to see the warnings...
    if (!is.na(as.numeric(cols[1])) || !is.na(as.numeric(cols[2])) || !is.na(as.numeric(cols[3])))
      stop("Dataframe df has a numeric column name; probably the required header row is missing.")
  )
  
  f <- df[[1]]
  p <- df[[2]]
  p0 <- df[[3]]
  
  if (!is.numeric(f))
    stop("The first column of the dataframe df, frequency, is not numeric.")
  if (!is.numeric(p))
    stop("The second column of the dataframe df, p, is not numeric.")
  if (!is.numeric(p0))
    stop("The third column of the dataframe df, p0, is not numeric.")
  if (any(is.na(f)))
    stop("The first column of the dataframe df, frequency, contains NA values (not allowed).")
  if (any(is.na(p)))
    stop("The second column of the dataframe df, p, contains NA values (not allowed).")
  if (any(is.na(p0)))
    stop("The third column of the dataframe df, p0, contains NA values (not allowed).")
  if (any(f < 0.0) || any(f > 1.0))
    stop("The first column of the dataframe df, frequency, contains values out of the required range [0,1].")
  if (any(p < 0))		# note that zero is allowed, although not recommended
    stop("The second column of the dataframe df, p, contains values < 0 (not allowed).")
  if (all(p == 0))		# not all can be zero, however
    stop("The second column of the dataframe df, p, contains all values == 0 (not allowed).")
  if (any(p0 <= 0))
    stop("The third column of the dataframe df, p0, contains values <= 0 (not allowed).")
  
  if (NROW(df) < 3)
    stop("At least three data rows are required, to constrain the fit.")
  
  #	Compute alpha values and trim
  #
  alpha <- 1 - (d0/d) * (p/p0)
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  
  if (sum(trim) < 3)
    stop("At least three data points are required after trimming the frequency range, to constrain the fit.")
  
  f_trimmed <- f[trim]
  alpha_trimmed <- alpha[trim]
  
  #	Compute the original McDonald-Kreitman alpha; we decided to use the trimmed data for this.
  #
  alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim]))			# from trimmed data
  #alpha_nonasymp <- 1 - (d0/d) * (sum(p)/sum(p0))						# from untrimmed data
  
  #	Fit models
  #
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
  
  if (length(mod1) == 0)
  {
    # try a deeper scan for a decent fit
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
  }
  
  tryCatch({
    mod2 <- lm(alpha_trimmed ~ f_trimmed)
  },
  error=function(cond) {})
  
  linear_better <- FALSE
  
  if ((length(mod1) == 0) || (AIC(mod2) < AIC(mod1)))
    linear_better <- TRUE
  
  if (!linear_better)
  {
    # if we're leaning toward the exponential model, check for ridiculously wide confidence intervals; sometimes
    # we should reject the exponential model for that reason, because it is basically just a linear model with
    # a "cheat" of a swing up or down to fit one additional data point perfectly, which is lame :->
    ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
    alpha_1_low <- ci_pred[6]
    alpha_1_high <- ci_pred[7]
    
    if ((alpha_1_low < -100) || (alpha_1_high > 100))
      linear_better <- TRUE
  }
  
  # Prepare for output and plotting
  full_seq <- seq(from=min(f), to=max(f), by=0.001)
  trimmed_seq <- seq(from=min(f_trimmed), to=max(f_trimmed), by=0.001)
  
  if (linear_better)
  {
    alpha_1_est <- predict(mod2, newdata=data.frame(f_trimmed=1.0))
    ci_pred <- predict(mod2, newdata=data.frame(f_trimmed=1.0), interval="confidence")	# we want confidence, not prediction
    alpha_1_low <- ci_pred[2]
    alpha_1_high <- ci_pred[3]
    const_a <- coef(mod2)["(Intercept)"]
    const_b <- coef(mod2)["f_trimmed"]
    const_c <- NA
    
    full_predicts <- predict(mod2, newdata=data.frame(f_trimmed=full_seq))
    trimmed_predicts <- predict(mod2, newdata=data.frame(f_trimmed=trimmed_seq))
    fit_color <- "red"
  }
  else
  {
    alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
    const_a <- coef(mod1)["const_a"]
    const_b <- coef(mod1)["const_b"]
    const_c <- coef(mod1)["const_c"]
    
    full_predicts <- predict(mod1, newdata=data.frame(f_trimmed=full_seq))
    trimmed_predicts <- predict(mod1, newdata=data.frame(f_trimmed=trimmed_seq))
    fit_color <- "red"
  }
  
  
  #	BEGIN OUTPUT
  #
  
  result_df <- data.frame(model=(if ((length(mod1) == 0) || linear_better) "linear" else "exponential"), a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=alpha_1_low, CI_high=alpha_1_high, alpha_original=alpha_nonasymp, row.names=NULL)
  
  if (output == "table")
  {
    # Just return the results dataframe; no plots or other output
    return(result_df)
  }
  else
  {
    cat("\nAsymptotic McDonald-Kreitman Evaluator: Results\n")
    cat("\n")
    cat("   from: http://benhaller.com/messerlab/asymptoticMK.html\n")
    cat("\n")
    cat("Analysis dataset:\n")
    cat("\n")
    cat(paste0("   d0 = ", format(d0, scientific=FALSE)), "\n")
    cat(paste0("   d = ", format(d, scientific=FALSE)), "\n")
    cat(paste0("   x interval = [", format(xlow, digits=3, nsmall=3, scientific=FALSE), ", ", format(xhigh, digits=3, nsmall=3, scientific=FALSE), "]"), "\n")
    cat(paste0("   f = (", paste0(format(f, scientific=FALSE), collapse=", "), ")\n"))
    cat(paste0("   p0 = (", paste0(format(p0, scientific=FALSE), collapse=", "), ")\n"))
    cat(paste0("   p = (", paste0(format(p, scientific=FALSE), collapse=", "), ")\n"))
    cat("\n")
    
    if ((length(mod1) == 0) || linear_better)
      cat("Fitted model: linear, alpha(x) = a + bx\n")
    else
      cat("Fitted model: exponential, alpha(x) = a + b * exp(−cx)\n")
    if (length(mod1) == 0)
      cat("   (The exponential fit failed to converge, usually because the data are not exponential in shape.)\n")
    
    cat("\n")
    
    cat(paste0("a\t", format(const_a, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("b\t", format(const_b, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("c\t", format(const_c, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("alpha_asymptotic\t", format(alpha_1_est, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("95% CI(lower)\t", format(alpha_1_low, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("95% CI(upper)\t", format(alpha_1_high, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("alpha_original\t", format(alpha_nonasymp, digits=5, nsmall=5, scientific=FALSE)), "\n")
    
    cat("\nIf you use this service, please cite our paper:\n\n   B.C. Haller, P.W. Messer. (2017). asymptoticMK: A web-based tool\n      for the asymptotic McDonald-Kreitman test. G3: Genes, Genomes,\n      Genetics 7(5), 1569-1575. doi:10.1534/g3.117.039693\n\nPlease let us know of any issues with this service at philipp {dot} messer <at> gmail [dot] com.  Thanks!\n\n")
    
    
    #	Output plots
    #
    
    #	PLOT 1: Frequency spectrum: p and p0 versus x
    #
    createOutputDevice(output=output, plotNum=1)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(c(p,p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="polymorphism counts")
    points(x=f, y=p0, col="black", pch=19, cex=0.7)
    points(x=f, y=p, col="red", pch=19, cex=0.7)
    legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
    
    closeOutputDevice(output=output, plotNum=1)
    
    
    #	PLOT 2: Frequency spectrum: p and p0 versus x
    #
    createOutputDevice(output=output, plotNum=2)
    
    normalized_p0 <- p0 / sum(p0)
    normalized_p <- p / sum(p)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(c(normalized_p,normalized_p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="normalized SFS")
    points(x=f, y=normalized_p0, col="black", pch=19, cex=0.7)
    points(x=f, y=normalized_p, col="red", pch=19, cex=0.7)
    legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
    
    closeOutputDevice(output=output, plotNum=2)
    
    
    #	PLOT 3: alpha(x) ~ x
    #
    createOutputDevice(output=output, plotNum=3)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(alpha), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
    points(x=f, y=alpha, col="black", pch=19, cex=0.7)
    
    closeOutputDevice(output=output, plotNum=3)
    
    
    #	PLOT 4: alpha(x) ~ x plus fit to that data
    #
    createOutputDevice(output=output, plotNum=4)
    
    yr <- na.omit(c(alpha, alpha_1_est, max(0.0, alpha_1_low), min(1.0, alpha_1_high), alpha_nonasymp))
    if (!is.na(true_alpha))
      yr <- c(yr, true_alpha)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(yr), cex.axis=0.8, cex.lab=1.0, ylim=arrYLims, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
    polygon(x=c(-1, -1, 2.0, 2.0), y=c(alpha_1_low, alpha_1_high, alpha_1_high, alpha_1_low), col="#DDDDDD", border=NA)
    
    if (!is.na(alpha_nonasymp))
      abline(h=alpha_nonasymp, lty=3, col="#777777")
    
    points(x=f, y=alpha, col=ifelse(trim, "black", "#999999"), pch=19, cex=0.7)
    abline(v=cutoff_f1, col="#5555FF")
    abline(v=cutoff_f2, col="#5555FF")
    
    if (!is.na(true_alpha))
      abline(h=true_alpha, col="#00AA00")
    
    lines(x=full_seq, full_predicts, col="#333333", lwd=1)
    lines(x=trimmed_seq, trimmed_predicts, col=fit_color, lwd=2)
    abline(h=alpha_1_est, lty=2, col=fit_color)
    
    box()
    closeOutputDevice(output=output, plotNum=4)
    
    # return the same table as for output="table", but using invisible()
    return(invisible(result_df))
  }
}


#######################################################################################
#
#	Now you have a function named asymptoticMK that you can call.  It looks like:
#
#		asymptoticMK(d0, d, xlow, xhigh, df, true_alpha=NA, output="table", ...)
#
#	The d0, d, xlow, xhigh, and df parameters are as defined in the web interface.
#	(For df, supply a dataframe with three numeric columns, for f, p, and p0; the
#	column names don't matter.)
#
#	You may supply true_alpha if you know it (i.e. from simulation data), but
#	typically you don't.
#
#	The output parameter tells the function how to produce output.  Supported
#	values at present are:
#
#		"table" : return a table of fit-related values, without plot generation
#		"quartz" : generate plot windows with the quartz() function, for OS X
#		"pdf" : save PDF files for the plots
#
#	It would be easy to add support for other devices if you wish; just modify the
#	createOutputDevice() function above to support your chosen device.
#
#
#	The asymptoticMK.Drosophila() function defined below shows how to use the
#	asymptoticMK() function to do the analysis for the Drosophila example dataset.
#
#######################################################################################


asymptoticMK.Drosophila <- function(xlow=0.1, xhigh=0.9, output="table", ...)
{
  # This would read the dataset from a file on disk:
  #
  # df <- read.csv("~/sample_polymorphism_levels.txt", sep="	")
  
  # Instead, to avoid having to deal with file locations and such, we just hard-code the dataset here
  f <- c(0.054, 0.062, 0.069, 0.077, 0.085, 0.092, 0.1, 0.108, 0.115, 0.123, 0.131, 0.138, 0.146, 0.154, 0.162, 0.169, 0.177, 0.185, 0.192, 0.2, 0.208, 0.215, 0.223, 0.231, 0.238, 0.246, 0.254, 0.262, 0.269, 0.277, 0.285, 0.292, 0.3, 0.308, 0.315, 0.323, 0.331, 0.338, 0.346, 0.354, 0.362, 0.369, 0.377, 0.385, 0.392, 0.4, 0.408, 0.415, 0.423, 0.431, 0.438, 0.446, 0.454, 0.462, 0.469, 0.477, 0.485, 0.492, 0.5, 0.508, 0.515, 0.523, 0.531, 0.538, 0.546, 0.554, 0.562, 0.569, 0.577, 0.585, 0.592, 0.6, 0.608, 0.615, 0.623, 0.631, 0.638, 0.646, 0.654, 0.662, 0.669, 0.677, 0.685, 0.692, 0.7, 0.708, 0.715, 0.723, 0.731, 0.738, 0.746, 0.754, 0.762, 0.769, 0.777, 0.785, 0.792, 0.8, 0.808, 0.815, 0.823, 0.831, 0.838, 0.846, 0.854, 0.862, 0.869, 0.877, 0.885, 0.892, 0.9, 0.908, 0.915, 0.923, 0.931, 0.938, 0.946)
  
  p <- c(1018, 865, 746, 620, 578, 466, 445, 412, 360, 338, 294, 299, 291, 235, 227, 262, 241, 211, 232, 214, 196, 183, 166, 126, 177, 159, 159, 139, 124, 132, 139, 126, 127, 111, 112, 132, 95, 130, 114, 111, 116, 85, 107, 96, 85, 89, 91, 77, 99, 94, 78, 94, 94, 101, 90, 88, 102, 67, 84, 67, 68, 84, 60, 73, 79, 65, 57, 73, 69, 84, 67, 69, 75, 78, 75, 64, 80, 70, 53, 56, 74, 67, 72, 59, 55, 58, 54, 74, 59, 57, 46, 54, 47, 62, 52, 68, 67, 55, 59, 63, 65, 68, 61, 58, 65, 61, 73, 45, 68, 72, 71, 66, 78, 73, 75, 97, 84)
  
  p0 <- c(2477, 2053, 1879, 1710, 1528, 1389, 1352, 1231, 1207, 1059, 1097, 987, 969, 953, 934, 887, 822, 845, 734, 755, 720, 730, 675, 642, 626, 657, 613, 562, 654, 560, 583, 616, 552, 559, 572, 524, 512, 504, 498, 505, 518, 485, 444, 471, 462, 448, 422, 435, 450, 476, 444, 457, 414, 394, 400, 433, 438, 379, 427, 420, 404, 379, 392, 403, 414, 377, 358, 400, 395, 374, 381, 361, 381, 376, 388, 384, 386, 394, 364, 350, 391, 360, 352, 359, 360, 371, 376, 343, 332, 354, 362, 366, 331, 377, 346, 372, 338, 354, 395, 405, 393, 363, 358, 384, 376, 369, 407, 420, 380, 436, 425, 423, 426, 451, 473, 508, 577)
  
  d <- 59570
  
  d0 <-  159058
  
  asymptoticMK(d0=d0, d=d, df=data.frame(f=f, p=p, p0=p0, row.names=NULL), xlow=xlow, xhigh=xhigh, output=output, ...)
}

# To just get a table of results, with no plots and no verbose analysis, do:
#asymptoticMK.Drosophila(output="table")

# To generate PDF files for the plots, just do:
#setwd("/Users/bhaller/Desktop")		# set the directory where you want the PDF files to appear
#asymptoticMK.Drosophila(output="pdf")

# To make plot windows on OS X using quartz(), do:
#asymptoticMK.Drosophila(output="quartz")


# Look at the asymptoticMK.Drosophila() function above to see how to run an analysis using your own dataset.
# Just substitute your own values for f, p, p0, d, and d0.




