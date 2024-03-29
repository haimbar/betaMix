% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betaMix.R
\name{summarizeClusters}
\alias{summarizeClusters}
\title{Show cluster characteristics.}
\usage{
summarizeClusters(clustersInfo)
}
\arguments{
\item{clustersInfo}{Obtained from graphComponents.}
}
\value{
A matrix with cluster number, number of nodes, and fivenum summaries for the degrees of nodes in the cluster, and the percentage of edges that are within the cluster.
}
\description{
Takes an object obtained from graphComponents and prints and returns summary statistics.
}
\examples{
\dontrun{
   data(SIM,package = "betaMix")
   res <- betaMix(betaMix::SIM, maxalpha = 1e-6,ppr = 0.01,subsamplesize = 30000, ind=TRUE)
   adjMat <- getAdjMat(res)
   SimComp <- graphComponents(adjMat)
   (summtab <- summarizeClusters(SimComp))
}
}
