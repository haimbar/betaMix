% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaMix.R
\name{shortestPathDistance}
\alias{shortestPathDistance}
\title{Return a Matrix with the shortest path distance between nodes (check up to numSteps.)}
\usage{
shortestPathDistance(AdjMat, numSteps = 0)
}
\arguments{
\item{AdjMat}{An adjacency Matrix (0/1).}

\item{numSteps}{The maximum number of edges between pairs of nodes. If numSteps=0, returns the input matrix. numSteps=1 adds neighbors of direct neighbors, etc.}
}
\value{
A Matrix containing the shortset paths between nodes i and j
}
\description{
return the adjacency matrix of expMat connecting neighbors up to numSteps away.
}
\examples{
\dontrun{
   data(SIM,package = "betaMix")
   res <- betaMix(betaMix::SIM, maxalpha = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
   adjMat <- getAdjMat(res)
   AdjMat <- shortestPathDistance(adjMat, numSteps=2)
   Matrix::image( (AdjMat>0)[1:200,1:200])
   adjMat1 <- AdjMat>0
   SimComp <- graphComponents(adjMat1)
   head(summarizeClusters(SimComp))
}
}
