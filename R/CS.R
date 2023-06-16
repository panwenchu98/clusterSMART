#'@title Compound Symmetry (Exchangeable) Matrix
#'
#'@description
#'This is a function for generating the compound symmetry (exchangeable) matrix
#'
#'@param rho The correlation parameter
#'@param n The size of the matrix
#'
#'@return An n*n matrix with diagonal term 1 and off-diagonal term rho
#'
#'@examples
#'library(clusterSMART)
#'CS(0.5,3)
#'
#'@export
#'

CS=function(rho,n){
  return((matrix(rho,n,n)+diag(1-rho,n)))
}
