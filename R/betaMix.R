
# The right-hand side of the support of the non-null beta component.
# Find x, such that the number of values from the null distribution to the
# right of x is at least (1-delta)*100% of the total number of values to the
# right of x. The result will be used as the maximum of the support of the 
# non-null distribution.
rightBetaThreshold <- function(x, z_jz, p0, eta, delta=1e-3) {
  ((1-delta)*mean(z_jz > x) - p0*(1-pbeta(x,eta,1/2)) )^2
}

# the Jacobian, to speed up the MLE calculation of a and b in the non-null
# component.
jacmle <- function(par, z_j0, m, xmax) {
  m <- m[which(z_j0 < xmax)]
  z_j0 <- z_j0[which(z_j0 < xmax)]
  z_j0 <- z_j0/xmax
  a <- par[1]
  b <- par[2]
  matrix(c(2*(digamma(a)-digamma(a+b) - 
                sum((1-m)*log(z_j0))/sum(1-m))*(trigamma(a)-trigamma(a+b)),
           2*(digamma(a)-digamma(a+b) - 
                sum((1-m)*log(z_j0))/sum(1-m))*(-trigamma(a+b)),
           2*(digamma(b)-digamma(a+b) - sum((1-m)*log(1-z_j0))/sum(1-m))*
             (-trigamma(a+b)),
           2*(digamma(b)-digamma(a+b) - sum((1-m)*log(1-z_j0))/sum(1-m))*
             (trigamma(a)-trigamma(a+b))),2,2)
  
}

# The maximum likelihood estimation of a,b of the nonnull component.
# nu (eta = (nu-1)/2) is estimated separately.
MLEfun <- function(par, z_j0, m, xmax) {
  m <- m[which(z_j0 < xmax)]
  z_j0 <- z_j0[which(z_j0 < xmax)]
  z_j0 <- z_j0/xmax
  a <- par[1]
  b <- par[2]
  y <- numeric(2)
  y[1] <- (digamma(a)-digamma(a+b) - sum((1-m)*log(z_j0))/sum(1-m))^2
  y[2] <- (digamma(b)-digamma(a+b) - sum((1-m)*log(1-z_j0))/sum(1-m))^2
  y  
}

# for estimation of the null parameter (if samples are dependent).
etafun <- function(eta,z_j0, m, xmax) {
  m <- m[which(z_j0 < xmax)]
  z_j0 <- z_j0[which(z_j0 < xmax)]
  z_j0 <- z_j0/xmax
  (digamma(eta)-digamma(eta+0.5) - sum(m*log(z_j0))/sum(m))
}


#' Fit a two-component beta mixture model to the matrix of all pairwise correlations obtained from an N*P matrix.
#'
#' From the pairwise correlations, the function calculates the statistics z_j=sin^2(arccos(cor(y_i,y_j))) and fits the two-component model using the EM algorithm.
#' @param M A matrix with N rows (samples) and P columns (variables).
#' @param tol The convergence threshold for the EM algorithm (default= the maximum of 1e-6 and 1/(P(P-1)/2)).
#' @param calcAcc The calculation accuracy threshold (to avoid values greater than 1 when calling asin) Default=1e-9.
#' @param delta The probability of Type I error (default=1e-3)
#' @param ppr The null posterior probability threshold (default=0.01)
#' @param mxcnt The maximum number of EM iterations (default=200)
#' @param ahat The initial value for the first parameter of the nonnull beta distribution (default=1).
#' @param bhat The initial value for the second parameter of the nonnull beta distribution (default=2).
#' @param nnmax The value of z_j above which it is not expected to find nonnull edges (default=0.999).
#' @param subsamplesize If greater than 20000, take a random sample of size subsamplesize to fit the model. Otherwise, use all the data (default=0, but for very large P, it's highly recommended to change this value.)
#' @param seed The random seed to use if selecting a subset with the subsamplesize parameter (default=912469).
#' @param ind Whether the N samples should be assumed to be independent (default=FALSE).
#' @param msg Whether to print intermediate output messages (default=TRUE).
#' @return A list with the following:
#' \itemize{
#' \item{angleMat} {A PxP matrix with angles between pairs of vectors.}
#' \item{z_j} {The statistics z_j=sin^2(angles).}
#' \item{m_0} {The posterior null probabilities.}
#' \item{p_0} {The estimated probability of the null component.}
#' \item{ahat} {The estimated first parameter of the nonnull beta component.}
#' \item{bhat} {The estimated second parameter of the nonnull beta component.}
#' \item{etahat} {If the samples are not assumed to be independent, this corresponds to the effective sample size, ESS=2*etahat+1}
#' \item {nonNullMax} {The estimated right-hand side of the support of the non-null component.}
#' \item {ppthr} {The estimated posterior probability threshold, under which all the z_j correspond to nonnull edges.}
#' }
#' @export
#' @import edgefinder stats nleqslv
#' @importFrom stats cor dbeta pbeta qbeta optimize
#' @examples
#' \donttest{
#'    data(SIM) # from the edgefinder package. Note that for betaMix is has 
#'              # to be transposed (variables correspond to columns, samples to rows)
#'    res <- betaMix(t(SIM), delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    plotFittedBetaMix(t(SIM), res)
#'    
#'    # can use some plotting functions from edgefinder to visualize the network
#'    adjMat <- sin(res$angleMat)^2 < res$ppthr
#'    diag(adjMat) <- FALSE
#'    graphComp <- graphComponents(adjMat)
#'    head(summarizeClusters(graphComp))
#'    plotCluster(adjMat, 1, graphComp, labels=TRUE, nodecol = "red")
#' }

betaMix <- function(M, tol=1e-6, calcAcc=1e-9, delta=1e-3, ppr=0.01, mxcnt=200,
                    ahat=1, bhat=2, nnmax=0.999, subsamplesize=0, seed=912469,
                    ind=FALSE, msg=TRUE) {
  if(msg) { cat("Fitting the model...\n") }
  N <- nrow(M)
  P <- ncol(M)
  etahat <- (N-1)/2
  angleMat <- acos(cor(M))
  z_j <- pmin(1-calcAcc, pmax(calcAcc, (sin(angleMat[which(lower.tri(angleMat))]))^2))
  p0 <- length(which(z_j > qbeta(0.1, etahat, 0.5)))/(0.9*length(z_j))
  z_jall <- c()
  if ((length(z_j) > subsamplesize) & (subsamplesize >= 20000)) {
    z_jall <- z_j
    set.seed(seed)
    z_j <- z_j[sample(length(z_j), subsamplesize)]
  }
  tol <- max(tol, 1/length(z_j))
  inNonNullSupport <- which(z_j < nnmax)
  p0f0 <- p0*dbeta(z_j, etahat, 0.5)
  p1f1 <- rep(0,length(z_j))
  p1f1[inNonNullSupport] <- (1-p0)*dbeta(z_j[inNonNullSupport]/nnmax, ahat, bhat)
  m0 <- pmax(0, pmin(1, p0f0/(p0f0+p1f1)))
  p0new <- mean(m0)
  cnt <- 0
  nonNullMax <- nnmax
  while (abs(p0-p0new) > tol & (cnt <- cnt+1) < mxcnt) {
    p0 <- p0new
    ests <- try(nleqslv(c(ahat,bhat), MLEfun, jac=jacmle,
                    z_j0=z_j, m=m0, xmax=nonNullMax)$x, silent=T)
    if (class(ests) == "try-error") {
      cat("betaMix error when estimating a and b:", ests,"\n")
    }
    
    ahat <- ests[1]
    bhat <- ests[2]
    if(ind) {
      etahat <- (N-1)/2
    } else {
      etahat <- try(uniroot(etafun,c(1,(N-1)/2), z_j0=z_j, m=m0,
                            xmax=nonNullMax)$root, silent=T)
      if (class(etahat) == "try-error") {
        cat("betaMix error when estimating eta:", etahat,"\n")
      }
    }
    p0f0 <- p0*dbeta(z_j,etahat,1/2)
    p1f1 <- rep(0, length(z_j))
    p1f1[inNonNullSupport] <- (1-p0)*dbeta(z_j[inNonNullSupport]/nonNullMax, ahat, bhat)
    m0 <- pmax(0, pmin(1, p0f0/(p0f0+p1f1)))
    nonNullMax <- optimize(rightBetaThreshold, interval=c(qbeta(0.99,etahat,0.5),1), z_jz=z_j, 
                           eta=etahat,p0=p0, delta=delta)$minimum
    inNonNullSupport <- which(z_j < nonNullMax)
    p0new <- mean(m0)
  }
  if (length(z_jall) > subsamplesize) {
    z_j <- z_jall
  }
  inNonNullSupport <- which(z_j < nonNullMax)
  p0f0 <- p0*dbeta(z_j,etahat,1/2)
  p1f1 <- rep(0,length(z_j))
  p1f1[inNonNullSupport] <- (1-p0)*dbeta(z_j[inNonNullSupport]/nonNullMax,ahat,bhat)
  m0 <- p0f0/(p0f0+p1f1) # the posterior null probability
  ppthr = max(min(z_j[which(m0 > ppr)]), qbeta(delta,etahat,1/2))
  p0 <- mean(m0)
  if(msg) { cat("Done.\n") }
  list(angleMat=angleMat, z_j=z_j, m0=m0,p0=p0, ahat=ahat, bhat=bhat, etahat=etahat,
       nonNullMax=nonNullMax, ppthr=ppthr)
}

#' Plot the histogram of the z_j and the fitted mixture distribution.
#' @param M A matrix with N rows (samples) and P columns (variables).
#' @param betaMixObj An object returned from betaMix()
#' @param yLim The maximum value on the y-axis (default=5)
#' @export
#' @examples
#' \donttest{
#'    data(SIM) # from the edgefinder package. Note that for betaMix is has to be transposed
#'    res <- betaMix(t(SIM), delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    plotFittedBetaMix(t(SIM), res)
#' }
plotFittedBetaMix <- function(M, betaMixObj, yLim=5) {
  with(betaMixObj, {
    N <- nrow(M)
    ccc <- seq(0.001,0.999,length=1000)
    hist(z_j,freq=F,breaks=300, border="grey",main="",ylim=c(0,yLim),
         xlab=expression(sin^{2} ~ (theta)))
    lines(ccc,(p0)*dbeta(ccc,(N-1)/2,0.5),col=3, lwd=4)
    lines(ccc,(p0)*dbeta(ccc,etahat,0.5),col=5, lwd=4,lty=2)
    lines(ccc,(1-p0)*dbeta(ccc/nonNullMax,ahat,bhat),lwd=1,col=2)
    lines(ccc, (1-p0)*dbeta(ccc/nonNullMax,ahat,bhat)+(p0)*dbeta(ccc,etahat,0.5), col=4, lwd=2)
    cccsig <- ccc[which(ccc<ppthr)]
    rect(0,0, ppthr, yLim, col='#FF7F5020', border = "orange")
    cat("ahat=",ahat, "\nbhat=",bhat, "\netahat=",etahat,
               "\nPost. Pr. threshold=", ppthr,"\n")
  })
}

#' Metabolite Expression data for the the dry seed group
#'
#' DrySeeds is a matrix with normalized metabolite data containing with 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name DrySeeds
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL

#' Metabolite Expression data for the the 6-hour imbibed seed group
#'
#' SixHourImbibed is a matrix with normalized metabolite data containing with 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name SixHourImbibed
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL
