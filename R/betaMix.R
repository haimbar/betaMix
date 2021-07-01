#' @useDynLib betaMix
#' @importFrom Rcpp sourceCpp evalCpp
#' @export calcCorr
NULL

#' Fit a two-component beta mixture model to the matrix of all pairwise correlations.
#'
#' From the pairwise correlations, the function calculates the statistics z_j=sin^2(arccos(cor(y_i,y_j))) and fits the two-component model using the EM algorithm.
#' @param M A matrix with N rows (samples) and P columns (variables).
#' @param dbname The sqlite database, if one is used to store the pairwise correlation data instead of using the cor function and storing the cor(M) matrix in memory (for situations in which P is very large).
#' @param tol The convergence threshold for the EM algorithm (default= the maximum of 1e-6 and 1/(P(P-1)/2)).
#' @param calcAcc The calculation accuracy threshold (to avoid values greater than 1 when calling asin) Default=1e-9.
#' @param delta The probability of Type I error (default=1e-3).
#' @param ppr The null posterior probability threshold (default=0.01).
#' @param mxcnt The maximum number of EM iterations (default=200).
#' @param ahat The initial value for the first parameter of the nonnull beta distribution (default=1).
#' @param bhat The initial value for the second parameter of the nonnull beta distribution (default=2).
#' @param nnmax The value of z_j above which it is not expected to find nonnull edges (default=0.999).
#' @param subsamplesize If greater than 20000, take a random sample of size subsamplesize to fit the model. Otherwise, use all the data (default=0, but for very large P, it's highly recommended to change this value.)
#' @param seed The random seed to use if selecting a subset with the subsamplesize parameter (default=912469).
#' @param ind Whether the N samples should be assumed to be independent (default=FALSE).
#' @param msg Whether to print intermediate output messages (default=TRUE).
#' @return A list with the following:
#' \itemize{
#' \item{angleMat} {A PxP matrix with angles between pairs of vectors. If the correlation data is stored in SQLite, then the returned value is database name.}
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
#' @import DBI stats nleqslv
#' @importFrom stats cor dbeta pbeta qbeta optimize
#' @examples
#' \donttest{
#'    data(SIM) # variables correspond to columns, samples to rows
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    plotFittedBetaMix(res)
#' }

betaMix <- function(M, dbname=NULL, tol=1e-6, calcAcc=1e-9, delta=1e-3,
                    ppr=0.01, mxcnt=200, ahat=1, bhat=2, nnmax=0.999,
                    subsamplesize=0, seed=912469, ind=FALSE, msg=TRUE) {
  if(msg) { cat("Generating the z_ij statistics...\n") }
  if (!is.null(dbname)) {
    if(!file.exists(dbname)) {
      message(paste("SQLite database",dbname,"not found!\n"))
      return(NULL)
    }
    con <- dbConnect(RSQLite::SQLite(), dbname)
    res <- dbSendQuery(con, "SELECT * FROM metadata")
    metadata <- dbFetch(res)
    dbClearResult(res)
    P <- metadata$p
    N <- metadata$n
    etahat <- (N-1)/2
    if (subsamplesize < 20000)
      subsamplesize <- 20000
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations ORDER BY random()  LIMIT %d",subsamplesize))
    subtable <- dbFetch(res)
    dbClearResult(res)
    z_j <- pmin(1 - calcAcc, pmax(calcAcc,subtable$zij))
    angleMat <- dbname
  } else {
    N <- nrow(M)
    P <- ncol(M)
    etahat <- (N-1)/2
    corM <- cor(M)
    if(any(is.na(corM)))
      corM[which(is.na(corM))] <- 0
    angleMat <- acos(corM)
    z_j <- pmin(1-calcAcc, pmax(calcAcc, (sin(angleMat[which(lower.tri(angleMat))]))^2))
    z_jall <- c()
    if ((length(z_j) > subsamplesize) & (subsamplesize >= 20000)) {
      z_jall <- z_j
      set.seed(seed)
      z_j <- z_j[sample(length(z_j), subsamplesize)]
    }
  }
  p0 <- length(which(z_j > qbeta(0.1, etahat, 0.5)))/(0.9*length(z_j))
  if(msg) { cat("Fitting the model...\n") }
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
      message("betaMix error when estimating a and b:", ests,"\n")
    }
    
    ahat <- ests[1]
    bhat <- ests[2]
    if(ind) {
      etahat <- (N-1)/2
    } else {
      etahat <- try(uniroot(etafun,c(1,(N-1)/2), z_j0=z_j, m=m0,
                            xmax=nonNullMax)$root, silent=T)
      if (class(etahat) == "try-error") {
        message("betaMix error when estimating eta:", etahat,"\n")
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
  if (!is.null(dbname)) {
    ppthr <- max(min(z_j[which(m0 > ppr)]), qbeta(delta, etahat, 1/2))
    p0 <- mean(m0)
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations  WHERE zij < %f",ppthr))
    selected <- dbFetch(res)
    edges <- nrow(selected)
    dbClearResult(res)
    dbDisconnect(con)
  } else{ 
    if (length(z_jall) > subsamplesize) {
      z_j <- z_jall
    }
    inNonNullSupport <- which(z_j < nonNullMax)
    p0f0 <- p0*dbeta(z_j,etahat,1/2)
    p1f1 <- rep(0,length(z_j))
    p1f1[inNonNullSupport] <- (1-p0)*dbeta(z_j[inNonNullSupport]/nonNullMax,ahat,bhat)
    m0 <- p0f0/(p0f0+p1f1) # the posterior null probability of all pairs
    ppthr = max(min(z_j[which(m0 > ppr)]), qbeta(delta,etahat,1/2))
    p0 <- mean(m0)
    edges <- (sum(sin(angleMat)^2 < ppthr) - P)/2
    
  }
  if(msg) { cat("Done.\n") }
  list(angleMat=angleMat, z_j=z_j, m0=m0,p0=p0, ahat=ahat, bhat=bhat, etahat=etahat,
       nonNullMax=nonNullMax, ppthr=ppthr, nodes=P, edges=edges)
}


#' Get the adjacency matrix if the pairwise correlation data is stored in a SQL database.
#' 
#' SQL storage should be used if P is large and the calculating and storing the correlation matrix will require too much memory.
#' @param res The returned value from the betaMix function if the correlation matrix is strored in memory.
#' @param dbname The sqlite database with the pairwise correlations, if the data is not stored in memory (for very large P). Default=NULL.
#' @param ppthr A threshold to select edges. Default=NULL, in which case the returned value from betaMix() is used (available in res).
#' @param signed If TRUE, the returned matrix will contain -1 for negatively correlated pairs. Otherwise, all correlated pairs will have 1 in the returned matrix (Default).
#' @param nodes An optional parameter which allows to create a sparse adjacency matrix containing only the neighbors of the selected nodes (default=NULL).
#' @export
#' @examples
#' \donttest{
#'    res <- betaMix(SIM, delta = 1e-5,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    image(adjMat[1:80,1:80])
#'    # See online documentation for usage when SQLite is used.
#' }
getAdjMat <- function(res, dbname=NULL, ppthr=NULL, signed=FALSE, nodes=NULL) {
  if (any(nodes > res$nodes)) {
    cat("Node number out of range. Max=",res$nodes,"\n")
    return(NULL)
  }
  if (is.null(ppthr))
    ppthr <- res$ppthr
  if(is.null(dbname)) {
    if (!is.null(nodes)){
      Atmp <- Matrix(res$angleMat[nodes,])
      nbrs <- which(rowSums(sin(Atmp)^2 < ppthr) > 0)
      if (length(setdiff(nbrs,nodes)) == 0) {
        cat("No neighbors found.\n")
        return(NULL)
      }
      selected <- sort(union(nodes,nbrs))
      angMat <- res$angleMat[selected, selected]
    } else {
      angMat <- res$angleMat
    }
    A <- sin(angMat)^2 < ppthr
    if(signed){
      negedges <- intersect(which(angMat > 0), which(angMat > pi/2))
      if (length(negedges) > 0)
        A[negedges] <- -1
    }
    diag(A) <- FALSE
    return(A)
  }
  if(!file.exists(dbname)) {
    message(paste("SQLite database",dbname,"not found!\n"))
    return(NULL)
  }
  con <- dbConnect(RSQLite::SQLite(), dbname)
  res <- dbSendQuery(con, "SELECT * FROM metadata")
  metadata <- dbFetch(res)
  dbClearResult(res)
  P <- metadata$p
  res <- dbSendQuery(con, "select name from predictors")
  pnames <- unlist(dbFetch(res))
  dbClearResult(res)
  A <- Matrix::Matrix(FALSE,P,P)
  colnames(A) <- pnames
  if (is.null(nodes)) {
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations  WHERE zij < %f",ppthr))
  } else {
    sset <- paste(nodes, collapse=",")
    res <- dbSendQuery(con, sprintf("SELECT * FROM correlations  WHERE (node1 IN (%s)  or node2 IN (%s)) AND zij < %f", sset, sset,ppthr))
  }
  selected <- dbFetch(res)
  dbClearResult(res)
  dbDisconnect(con)
  if(nrow(selected) == 0) {
    message("Zero rows selected.\n")
    return(NULL)
  }
  selectedCells <- P*(selected$node1-1) + selected$node2
  A[selectedCells]  <- TRUE
  if(signed) {
    A[selectedCells[which(selected$corr < 0)]] <- -1
  }
  return(A + Matrix::t(A))
}


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

# The maximum likelihood estimation of the null parameter (if samples are dependent).
etafun <- function(eta,z_j0, m, xmax) {
  m <- m[which(z_j0 < xmax)]
  z_j0 <- z_j0[which(z_j0 < xmax)]
  z_j0 <- z_j0/xmax
  (digamma(eta)-digamma(eta+0.5) - sum(m*log(z_j0))/sum(m))
}


#' Plot the histogram of the z_j and the fitted mixture distribution.
#' 
#' @param betaMixObj An object returned from betaMix()
#' @param yLim The maximum value on the y-axis (default=5)
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-5,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    plotFittedBetaMix(res)
#' }
plotFittedBetaMix <- function(betaMixObj, yLim=5) {
  with(betaMixObj, {
    ccc <- seq(0.001,0.999,length=1000)
    hist(z_j,freq=F,breaks=300, border="grey",main="",ylim=c(0,yLim),
         xlim=c(0,1),xlab=expression(sin^{2} ~ (theta)))
    lines(ccc,(p0)*dbeta(ccc,(betaMixObj$nodes-1)/2,0.5),col=3, lwd=4)
    lines(ccc,(p0)*dbeta(ccc,etahat,0.5),col=5, lwd=4,lty=2)
    lines(ccc,(1-p0)*dbeta(ccc/nonNullMax,ahat,bhat),lwd=1,col=2)
    lines(ccc, (1-p0)*dbeta(ccc/nonNullMax,ahat,bhat)+(p0)*dbeta(ccc,etahat,0.5), col=4, lwd=2)
    cccsig <- ccc[which(ccc<ppthr)]
    rect(0,0, ppthr, yLim, col='#FF7F5020', border = "orange")
  })
}

#' Print a short summary of the fitted mixture model.
#'
#' Show the number of nodes, the number of possible and detected edges, the estimated proportion of uncorrelated pairs.
#' @param betamixobj The object (list) returned from the betaMix function.
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    shortSummary(res)
#' }
shortSummary <- function(betamixobj) {
  with(betamixobj,{
    cat("ahat=",ahat, "\nbhat=",bhat, "\netahat=",etahat,
        "\nPost. Pr. threshold=", ppthr,"\n")
    cat("No. nodes =", prettyNum(nodes,big.mark = ","),"\n")
    cat("Max no. edges =", prettyNum(choose(nodes, 2),big.mark = ","),"\n")
    cat("No. edges detected =", prettyNum(edges,big.mark = ","),"\n")
    cat("p0 =",format(p0,digits=3),"\n")
  })
}



#' Find spherical caps with more than one node.
#'
#' Takes a PxP adjacency Matrix as input and find apherical caps with more than one node in R^n. Cap centers can only appear in one cap, but other nodes are allowed to appear in multiple caps.
#' @param A An adjacency Matrix(0/1), zeros on the main diagonal.
#' @return A data frame with the following columns:
#' \itemize{
#' \item{node} {Node number (row/column in A.)}
#' \item{capNum} {The cap number.}
#' \item{isCtr} {1 if the node is the center of the cap, 0 otherwise.}
#' \item{deg} {Node degree.}
#' \item{cc} {Clustering coefficient.}
#' }
#' @importFrom Matrix Matrix rowSums
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    caps <- sphericalCaps(adjMat)
#'    head(caps)
#' }
sphericalCaps <- function(A) {
  stopifnot(grep("Matrix", class(A)) > 0)
  A <- abs(A)
  diag(A) <- FALSE
  deg <- Matrix::rowSums(A)
  possibleCtrs <- which(deg > 0)
  if (length(possibleCtrs) == 0) {
    cat("No edges!\n")
    return(NULL)
  }
  CC <- clusteringCoef(A)
  retdf <- data.frame(matrix(vector(), 0, 5))
  orddeg <- order(deg, decreasing = TRUE)
  capNum <- 1
  while (max(deg[orddeg]) > 0) {
    capCtr <- orddeg[1]
    if(deg[capCtr] == 0)
      break
    nbrs <- setdiff(which(A[capCtr,] != 0), capCtr)
    orddeg <- setdiff(orddeg, union(capCtr, nbrs))
    retdf <- rbind(retdf, cbind(c(capCtr,nbrs), 
                                rep(capNum,length(c(capCtr,nbrs))),
                                c(1,rep(0,length(nbrs))),
                                deg[c(capCtr,nbrs)],
                                CC[c(capCtr,nbrs)]))
    capNum <- capNum + 1
    if(length(orddeg) == 0)
      break
  }
  colnames(retdf) <- c("node", "capNum", "isCtr", "deg", "cc")
  return(retdf)
}


#' Find clusters, and return node characteristics.
#'
#' Take an adjacency Matrix as input and find clusters. For each node, find the degree and clustering coefficient (CC). Then, calculate a centrality measure (type\*CC+1)\*deg. For type=0, it's just the degree. Note that setting type=1 means that we assign a higher value to nodes that not only have many neighbors, but the neighbors are highly interconnected. For example, suppose we have two components with k nodes, one has a star shape, and the other is a complete graph. With type=0 both graphs will get the same value, but with type=1 the complete graph will be picked by the algorithm first. Setting type to a negative value gives CC\*deg as the centrality measure.
#' @param A An adjacency Matrix(0/1).
#' @param minCtr The minimum centrality value to be considered for a cluster center (default=5).
#' @param type Determines how the centrality measure is computed.
#' @return A data frame with the following columns:
#' \itemize{
#'  \item{labels} {Node label (e.g. gene names).}
#' \item{degree} {Node degree.}
#' \item{cc} {Node clustering coefficient.}
#' \item{ctr} {Node centrality measure: (type\*CC+1)\*deg, or CC\*deg if type is negative.}
#' \item{clustNo} {Cluster number.}
#' \item {iscenter} {1 for the node was chosen as the cluster's center, 0 otherwise.}
#' \item {intEdges} {Number of edges from the node to nodes in the same cluster.}
#' \item {extEdges} {Number of edges from the node to nodes NOT in the same cluster.}
#' \item {distCenter} {Standardized Manhattan distance to the central node.}
#' }
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    SimComp <- graphComponents(adjMat)
#'    head(SimComp)
#' }
graphComponents <- function(A, minCtr=5, type=1) {
  stopifnot(grep("Matrix", class(A)) > 0)
  A <- abs(A)
  Vn <- ncol(A)
  ctrs <- rep(2*Vn, Vn)
  labels <- 1:Vn
  if(!is.null(rownames(A)))
    labels <- rownames(A)
  deg <- Matrix::rowSums(A)
  CC <- clusteringCoef(A)
  ctrs <- (type*CC+1)*deg
  if (type < 0)
    ctrs <- CC*deg
  clustersInfo <- data.frame(labels=labels, degree=deg, cc=CC, ctr=ctrs,
                             clustNo=rep(0,Vn), iscenter=rep(0,Vn),
                             intEdges=rep(0,Vn), extEdges=rep(0,Vn),
                             distCenter=rep(0,Vn))
  clustNo <- 1
  clustered <- which(deg < 1)
  while(length(clustered) < Vn) {
    notInCluster <- setdiff(1:Vn, clustered)
    if (max(ctrs[notInCluster]) < minCtr)
      return(clustersInfo)
    ctrnode <- notInCluster[which.max(ctrs[notInCluster])]
    # candidate cluster neighbors
    nbrs <- setdiff(sort(c(ctrnode, which(A[ctrnode,] != 0))), clustered)
    if(length(nbrs) > 1) {
      if (length(nbrs) > minCtr) {
        clustersInfo$iscenter[ctrnode] <- 1
        clustersInfo$clustNo[union(ctrnode,nbrs)] <- clustNo
        clustersInfo$intEdges[nbrs] <- Matrix::rowSums(A[nbrs,nbrs])
        if (length(nbrs) < ncol(A)) {
          clustersInfo$extEdges[nbrs] <- Matrix::rowSums(A[nbrs,-nbrs])
        } else {
          clustersInfo$extEdges[nbrs] <- 0
        }
        for (i in 1:length(nbrs)) {
          clustersInfo$distCenter[nbrs[i]] <- mean(xor(A[ctrnode,], A[nbrs[i],]))
        }
        clustNo <- clustNo + 1
      }  else {
        nbrs <- c()
      }
    } else {
      nbrs <- c()
    }
    clustered <- union(clustered, c(nbrs, ctrnode))
  }
  return(clustersInfo)
}


#' Show cluster characteristics.
#'
#' Takes an object obtained from graphComponents and prints and returns summary statistics.
#' @param clustersInfo Obtained from graphComponents.
#' @return A matrix with cluster number, number of nodes, and fivenum summaries for the degrees of nodes in the cluster, and the percentage of edges that are within the cluster.
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    SimComp <- graphComponents(adjMat)
#'    (summtab <- summarizeClusters(SimComp))
#' }
summarizeClusters <- function(clustersInfo) {
  cat("Num of nodes:", nrow(clustersInfo),"\n")
  cat("Num of edges:", sum(clustersInfo$degree)/2,"\n")
  cat("Num of clusters:", max(clustersInfo$clustNo),"\n")
  cat("Num of unclustered nodes:", length(which(clustersInfo$clustNo == 0)),"\n")
  percentInCluster <- clustersInfo$intEdges/clustersInfo$degree
  percentInCluster[which(clustersInfo$degree == 0)] <- 0
  if (max(clustersInfo$clustNo) == 0)
    return(NULL)
  tab <- matrix(0,nrow=max(clustersInfo$clustNo),ncol=12)
  for (cnum in 1:max(clustersInfo$clustNo)) {
    tmpclusterInfo <- clustersInfo[which(clustersInfo$clustNo == cnum),]
    tab[cnum,] <- c(cnum,nrow(tmpclusterInfo), fivenum(tmpclusterInfo$degree),
                    fivenum(percentInCluster[which(clustersInfo$clustNo == cnum)]))
  }
  colnames(tab) <- c("Cluster","Nodes","degreeMin","degreeQ25","degreeMedian",
                     "degreeQ75","degreeMax","pctInClstMin","pctInClstQ25",
                     "pctInClstMedian", "pctInClstQ75","pctInClstMax")
  tab
}


#' Reduce an adjacency matrix to include only central nodes and nodes not assigned to any cluster.
#' 
#' Return an adjacency matrix after collapsing clusters into their central nodes.
#' @param A An adjacency Matrix.
#' @param clustersInfo Obtained from graphComponents
#' @return A weighted adjacency matrix between clusters and unclustered nodes.
#' @export
#' @examples
#' \donttest{
#'    data(DrySeeds)
#'    res <- betaMix(DrySeeds, delta = 1e-4,ppr = 0.01, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    rownames(adjMat) = colnames(DrySeeds)
#'    SimComp <- graphComponents(adjMat)
#'    Adj1 <- collapsedGraph(adjMat, SimComp) > 0
#'    # with the igraph package, we can do:
#'    # plot(graph.adjacency(Adj1, mode="undirected"), vertex.label=rownames(Adj1),
#'    # vertex.label.cex=0.7, vertex.size=0, vertex.frame.color="white",
#'    # vertex.label.color='blue', edge.color='grey80',  asp=1)
#' }
collapsedGraph <- function(A, clustersInfo) {
  collDim <- length(which(clustersInfo$clustNo == 0)) + max(clustersInfo$clustNo)
  collA <- Matrix::Matrix(0, ncol=collDim, nrow=collDim)
  inCluster <- which(clustersInfo$clustNo > 0)
  notInCluster <- which(clustersInfo$clustNo == 0)
  if (length(notInCluster) > 0) {
    collA[1:length(notInCluster), 1:length(notInCluster)] <- A[notInCluster, notInCluster]>0
  }
  if (length(rownames(A)) != nrow(A)) {
    rownames(A) <- 1:nrow(A)
  }
  rownames(collA) <- c(rownames(A)[notInCluster],
                       paste0("CLS",1:max(clustersInfo$clustNo)))
  for (i in 1:max(clustersInfo$clustNo)) {
    Ci <- which(clustersInfo$clustNo == i)
    if (length(notInCluster) > 0) {
      Atmp <- matrix(A[notInCluster,which(clustersInfo$clustNo==i)],
                     nrow=length(notInCluster), ncol=length(which(clustersInfo$clustNo==i)))
      collA[i+length(notInCluster),1:length(notInCluster)] <- Matrix::rowSums(Atmp)
    }
    if (i < max(clustersInfo$clustNo)) {
      for (j in (i+1):max(clustersInfo$clustNo)) {
        Cj <- which(clustersInfo$clustNo == j)
        collA[i+length(notInCluster),j+length(notInCluster)] <- sum(A[Ci,Cj])
      }
    }
  }
  collA + Matrix::t(collA)
}


#' Calculate the clustering coefficient of each node.
#'
#' @param A an adjacency Matrix (0/1).
#' @return A vector with the clustering coefficient of each node.
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    clusteringCoef(adjMat)
#' }
#'
clusteringCoef <- function(A) {
  rsum <- Matrix::rowSums(A)
  cc <- rep(0,nrow(A))
  for (i in 1:nrow(A)) {
    if (rsum[i] <= 1)
      cc[i] <- 0
    else {
      nbrs <- which(A[i,] == 1)
      At <- A[nbrs, nbrs]
      cc[i] <- 0.5*sum(At)/choose(rsum[i],2)
    }
  }
  cc
}


#' Plot the degree of nodes versus the degree times the clustering coefficient.
#'
#' The x-axis represents the number of neighbors of each node, and the y-axis represents the proportion of neighbors which are connected to each other.
#' @param betamixobj The object (list) returned by betaMix
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param highlightNodes A vector of node-numbers which will be shown in red. Default is NULL.
#' @export
#' @import stats graphics
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    SimComp <- graphComponents(adjMat)
#'    plotDegCC(res,SimComp)
#' }
plotDegCC <- function(betamixobj, clusterInfo=NULL, highlightNodes=NULL) {
  if (is.null(clusterInfo))
    clusterInfo <-  graphComponents(betamixobj$AdjMat)
  cc0 <- clusterInfo$cc
  deg0 <- clusterInfo$degree
  plot(deg0, deg0*cc0,axes=F,xlim=c(0,max(deg0)),
       ylim=c(0,1.1*max(deg0*cc0)),main="",
       xlab=bquote("degree"),ylab=bquote("CC*degree"),
       col="thistle",pch=24,cex=0.5); axis(1); axis(2)
  grid(); abline(0,1,col="seagreen1", lwd=2)
  if (!is.null(highlightNodes))
    points(deg0[highlightNodes],(deg0*cc0)[highlightNodes],col=2,pch=24,cex=0.5)
}


#' Edge-indicator bitmap plot.
#'
#' Plot a bitmap in which a black dot corresponds to a pair of highly correlated genes (an edge in the graph).
#' The default is to show the nodes according to their order in the input.
#' By setting orderByDegree=T as below, it is possible to change the order and cluster them, and show them in increasing degree order (from left to right.)
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param orderByCluster If false, show the bitmap is the original node order. If TRUE, show nodes by clusters, and sort by distance from the center of the cluster.
#' @param showMinDegree Non-negative integer indicating the minimum degree of nodes that should be displayed. Default=0 (all nodes).
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    SimComp <- graphComponents(adjMat)
#'    plotBitmapCC(adjMat)
#'    plotBitmapCC(adjMat, SimComp, orderByCluster=TRUE)
#' }
plotBitmapCC <- function(AdjMat, clusterInfo=NULL, orderByCluster=FALSE, showMinDegree=0) {
  if(!is.null(clusterInfo))
    orderByCluster <- TRUE
  if (orderByCluster) {
    if (is.null(clusterInfo))
      clusterInfo <- graphComponents(AdjMat)
    nodeOrder <- order(clusterInfo$clustNo,clusterInfo$distCenter)
    AdjMat <- AdjMat[nodeOrder, nodeOrder]
  }
  showNodes <- which(Matrix::rowSums(AdjMat) >= showMinDegree)
  Matrix::image(AdjMat[showNodes, showNodes])
}


#' Plot cluster network
#'
#' Plot a cluster network with all the nodes and edges - the central node is marked by a black circle. The radius of each point corresponds to its degree. The opacity corresponds to the percentage of edges from the node that is in the cluster (the darker it is, the larger the percentage of edges is within the cluster.) The distance from the center corresponds to the relative dissimilarity with the central node. This is computed as the number of neighbors the node and the central node do not have in common.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clustNo The chosen cluster.
#' @param clusterInfo Obtained from graphComponents.
#' @param labels If set to TRUE, show node names (default=FALSE).
#' @param nodecol The color(s) of the nodes. Can be a single value or a vector of length equal to the number of rows in AdjMat
#' @param labelsize Text size of node labels.
#' @param figtitle The title of the plot (default=NULL).
#' @param edgecols The colors to be used for edges. Default="grey88". If one value is given, all edges will be drawn using this color. If edgecol contains two valid colors, the first is used for positive correlations, and the second for negative ones.
#' @export
#' @importFrom grDevices col2rgb colours rgb
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-5,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res, signed=TRUE)
#'    SimComp <- graphComponents(adjMat)
#'    plotCluster(adjMat, 2, SimComp)
#'    plotCluster(adjMat, 2, SimComp, edgecols = c("blue","orange"), labels=TRUE)
#' }
plotCluster <- function(AdjMat, clustNo, clusterInfo=NULL, labels=FALSE, nodecol="blue", labelsize=1, figtitle=NULL, edgecols="grey88") {
  if(is.null(clusterInfo))
    clusterInfo <- graphComponents(AdjMat)
  if(length(nodecol) < nrow(AdjMat))
    nodecol <- rep(nodecol[1],length=nrow(AdjMat))
  ids <- which(clusterInfo$clustNo == clustNo)
  if (length(ids) > 0) {
    tmpA <- AdjMat[ids,ids]
    tmpclusterInfo <- clusterInfo[ids,]
    rads <- round(10*tmpclusterInfo$distCenter/
                    max(0.1,max(tmpclusterInfo$distCenter)))
    thetas <- rep(0,length(rads))
    intvls <- findInterval(rads,seq(1,10))
    for (intvl in unique(sort(intvls))) {
      pts <- which(intvls == intvl)
      thetas[pts] <- 3*intvl*pi/max(0.01,max(intvls))+seq(0,1.9*pi,length=length(pts))
    }
    sizes <- pmax(0.3,tmpclusterInfo$degree/max(tmpclusterInfo$degree))
    opacity <- 0.25+tmpclusterInfo$intEdges/tmpclusterInfo$degree
    opacity <- opacity/max(opacity)
    nodecol <- rgb(t(col2rgb(nodecol)/255),alpha=opacity)[ids]
    plot(rads*cos(thetas), rads*sin(thetas),cex=sizes*3, pch=19,axes=F,
         xlab="",ylab="",col=nodecol, main=figtitle,
         ylim=c(min(rads*sin(thetas)), 1.1*max(rads*sin(thetas))))
    for (i in 1:ncol(tmpA)) {
      nbrs <- setdiff(which(abs(tmpA[i,]) == 1), 1:i)
      if(length(nbrs) > 0) {
        edgecol <- rep(edgecols[1], ncol(tmpA))
        if (edgecols[2] %in% colours()) {
          edgecol[which(tmpA[i,nbrs] == -1)] <- edgecols[2]
        }
        for (j in nbrs) {
          lines(c(rads[i]*cos(thetas[i]), rads[j]*cos(thetas[j])),
                c(rads[i]*sin(thetas[i]), rads[j]*sin(thetas[j])),
                col=edgecol[j], lwd=0.5)
        }
      }
    }
    points(rads*cos(thetas), rads*sin(thetas),cex=sizes*3, pch=19, col=nodecol)
    if (labels)
      text(rads*cos(thetas), rads*sin(thetas), tmpclusterInfo$labels, pos=3, cex=labelsize)
    ctr <- which(tmpclusterInfo$iscenter==1)
    points(rads[ctr]*cos(thetas[ctr]), rads[ctr]*sin(thetas[ctr]),pch=21,
           cex=sizes[ctr]*3, col="black",lwd=2)
  } else {
    cat("Invalid cluster number\n")
  }
}


#' Return a Matrix with the shortest path distance between nodes (check up to numSteps.)
#'
#' return the adjacency matrix of expMat connecting neighbors up to numSteps away.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param numSteps The maximum number of edges between pairs of nodes. If numSteps=0, returns the input matrix. numSteps=1 adds neighbors of direct neighbors, etc.
#' @return A Matrix containing the shortset paths between nodes i and j
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    res <- betaMix(SIM, delta = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
#'    adjMat <- getAdjMat(res)
#'    AdjMat <- shortestPathDistance(adjMat, numSteps=2)
#'    Matrix::image( (AdjMat>0)[1:200,1:200])
#'    adjMat1 <- AdjMat>0
#'    SimComp <- graphComponents(adjMat1)
#'    head(summarizeClusters(SimComp))
#' }
shortestPathDistance <- function(AdjMat, numSteps=0) {
  degs <- 1:ncol(AdjMat)
  if (numSteps == 0)
    return(AdjMat)
  An <- Ap <- minDist <- AdjMat
  for (i in 1:numSteps) {
    An <- Ap%*%AdjMat
    if (sum((An | Ap) - (An & Ap)) == 0)
      break
    minDist[(An > 0) & (Ap == 0) & (minDist == 0)] <- i
    Ap <- An
  }
  rownames(minDist) <- colnames(minDist) <- rownames(AdjMat)
  minDist
}



#' Metabolite Expression data for the the dry seed group.
#'
#' DrySeeds is a matrix with normalized metabolite data containing with 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name DrySeeds
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL

#' Metabolite Expression data for the the 6-hour imbibed seed group.
#'
#' SixHourImbibed is a matrix with normalized metabolite data containing with 68 named metabolites, and 50 samples.
#'
#' @docType data
#' @keywords datasets
#' @name SixHourImbibed
#' @format A matrix with 50 rows and 68 columns
#' @references \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000561&StudyType=MS&ResultType=1}
NULL


#' Simulated gene Expression data using the huge package
#'
#' SIM is a simulated dataset with a hub structure, consisting of 1000 nodes and 50 hubs.
#'
#' @docType data
#' @keywords datasets
#' @name SIM
#' @usage data(SIM)
#' @format A 1000 by 200 matrix, representing 50 hubs
NULL
