#'@title Change from Treatment Pathway to (A1,R,A2)
#'
#'@description
#'This is a function for extracting the randomization result and response status (A1,R,A2) from treatment pathway label under prototypical SMART.
#'
#'@param L Length n vector with values in {1,2,...,6}, represents the label of treatment pathway.
#'
#'@return A n*3 matrix of the extracted results:\tabular{ll}{
#'  \code{A1}\tab The result of the cluster-level randomization for the first-stage treatment.\cr
#'  \tab \cr
#'  \code{R}\tab The cluster-level response status to the first-stage treatment. \cr
#'  \tab \cr
#'  \code{A2}\tab The result of the cluster-level randomization for the second-stage treatment.\cr
#'  }
#'
#'@examples
#'library(clusterSMART)
#'L=c(1:6)
#'pathway_to_response(L)
#'
#'@export
#'

pathway_to_response=function(L){
  if (FALSE %in% (unique(L) %in% c(1,2,3,4,5,6))) stop("Error: value of L should be in 1,2,...,6")
  A1=ifelse(L>=4,-1,1)
  R=ifelse(L %in% c(1,4),1,0)
  A2=ifelse(L %in% c(1,4),NA,ifelse(L %in% c(2,5),1,-1))
  result=cbind(A1,R,A2)
  return(result)
}
