% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response_to_pathway.R
\name{response_to_pathway}
\alias{response_to_pathway}
\title{Change from (A1,R,A2) to Treatment Pathway}
\usage{
response_to_pathway(A1, R, A2)
}
\arguments{
\item{A1}{Length n vector with value in {-1,1}, represents the result of the cluster-level randomization for the first-stage treatment.}

\item{R}{Length n vector with value in {0,1}, represents the cluster-level response status to the first-stage treatment.}

\item{A2}{Length n vector with values in {-1,1,NA}, represents the result of the cluster-level randomization for the second-stage treatment.}
}
\value{
Length n vector with values in {1,2,...,6}, represents the label of treatment pathway.
}
\description{
This is a function for merging the randomization result and response status (A1,R,A2) to treatment pathway label under prototypical SMART.
}
\examples{
library(clusterSMART)
A1=c(1,1,1,-1,-1,-1)
R=c(1,0,0,1,0,0)
A2=c(NA,1,-1,NA,1,-1)
response_to_pathway(A1,R,A2)

}
