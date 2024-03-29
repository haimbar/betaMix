% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaMix.R
\name{sphericalCaps}
\alias{sphericalCaps}
\title{Find spherical caps with more than one node.}
\usage{
sphericalCaps(A)
}
\arguments{
\item{A}{An adjacency Matrix(0/1), zeros on the main diagonal.}
}
\value{
A data frame with the following columns:
\itemize{
\item{node} {Node number (row/column in A.)}
\item{capNum} {The cap number.}
\item{isCtr} {1 if the node is the center of the cap, 0 otherwise.}
\item{deg} {Node degree.}
\item{cc} {Clustering coefficient.}
}
}
\description{
Takes a PxP adjacency Matrix as input and find apherical caps with more than one node in R^n. Cap centers can only appear in one cap, but other nodes are allowed to appear in multiple caps.
}
\examples{
\dontrun{
   data(SIM,package = "betaMix")
   res <- betaMix(betaMix::SIM, maxalpha = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
   adjMat <- getAdjMat(res)
   caps <- sphericalCaps(adjMat)
   head(caps)
   am <- getAdjMat(res, signed = TRUE, nodes=caps$node[which(caps$capNum==1)])
   plotCluster(am,1,edgecols = c("blue","red"), labels=TRUE)
}
}
