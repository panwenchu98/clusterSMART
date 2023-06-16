#'@title Autocorrelation Matrix
#'
#'@description
#'This is a function for generating the autocorrelated matrix
#'
#'@param rho The correlation parameter
#'@param n The size of the matrix
#'
#'@return An n*n matrix with the (i,j) term having value rho^abs(i-j)
#'
#'@examples
#'library(clusterSMART)
#'AR1(0.5,3)
#'
#'@importFrom stats toeplitz
#'@export
#'


AR1=function(rho,n){
  result=toeplitz(rho^(0:(n-1)))
  return(result)
}
