#'@title Change from (A1,R,A2) to Treatment Pathway
#'
#'@description
#'This is a function for merging the randomization result and response status (A1,R,A2) to treatment pathway label under prototypical SMART.
#'
#'@param A1 Length n vector with value in {-1,1}, represents the result of the cluster-level randomization for the first-stage treatment.
#'@param R Length n vector with value in {0,1}, represents the cluster-level response status to the first-stage treatment.
#'@param A2 Length n vector with values in {-1,1,NA}, represents the result of the cluster-level randomization for the second-stage treatment.
#'
#'@return Length n vector with values in {1,2,...,6}, represents the label of treatment pathway.
#'
#'@examples
#'library(clusterSMART)
#'A1=c(1,1,1,-1,-1,-1)
#'R=c(1,0,0,1,0,0)
#'A2=c(NA,1,-1,NA,1,-1)
#'response_to_pathway(A1,R,A2)
#'
#'@export
#'

response_to_pathway=function(A1,R,A2){
  n=length(A1)
  if (length(R)!=n) stop("Error: Variables array not of same length")
  if (length(A2)!=n) stop("Error: Variables array not of same length")
  A2=ifelse(is.na(A2),0,A2)
  if (!(all.equal(unique(A2[which(R==1)]),0)==TRUE)) stop("Error: Randomizing responders of first-stage treatment")
  label=1+(1-A1)*3/2+(3*A2-1)*A2/2
  return(label)
}
