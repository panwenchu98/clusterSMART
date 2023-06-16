#'@title Changing Conditional Parameters to Marginal Parameters
#'
#'@description
#'This is a function for changing the conditional (treatment-pathway level) parameters to marginal (Adaptive-Intervention level) parameters in a prototypical SMART.
#'The parameter is defined without considering the baseline covariates
#'
#'@param p1,p2 The response rate of the first-stage treatment in each randomization branch
#'@param conditional_paras A matrix with 6 rows, the first two entries of each row represent the expected mean and variance of each treatment pathway
#'
#'@return A 4*2 matrix, each row contains the marginal expected mean and variance of each Adaptive Intervention
#'
#'@examples
#'library(clusterSMART)
#'conditional_paras_original=matrix(c(6,8,4,4,4,4,100,81,64,49,36,36,0.20,0.18,0.16,0.14,0.12,0.10),6,3)
#'p1=0.25; p2=0.55
#'conditional_to_marginal(p1,p2,conditional_paras_original)
#'@export
#'

#change from conditional parameters (treatment pathway level mean and variance) to marginal parameters (AI level)
conditional_to_marginal=function(p1,p2,conditional_paras){
  ans=matrix(0,4,2)
  colnames(ans)=c("mu","sigma2")
  ans[1,1]=conditional_paras[1,1]*p1+conditional_paras[2,1]*(1-p1)
  ans[2,1]=conditional_paras[1,1]*p1+conditional_paras[3,1]*(1-p1)
  ans[3,1]=conditional_paras[4,1]*p2+conditional_paras[5,1]*(1-p2)
  ans[4,1]=conditional_paras[4,1]*p2+conditional_paras[6,1]*(1-p2)

  ans[1,2]=conditional_paras[1,2]*p1+conditional_paras[2,2]*(1-p1)+p1*(1-p1)*(conditional_paras[1,1]-conditional_paras[2,1])^2
  ans[2,2]=conditional_paras[1,2]*p1+conditional_paras[3,2]*(1-p1)+p1*(1-p1)*(conditional_paras[1,1]-conditional_paras[3,1])^2
  ans[3,2]=conditional_paras[4,2]*p2+conditional_paras[5,2]*(1-p2)+p2*(1-p2)*(conditional_paras[4,1]-conditional_paras[5,1])^2
  ans[4,2]=conditional_paras[4,2]*p2+conditional_paras[6,2]*(1-p2)+p2*(1-p2)*(conditional_paras[4,1]-conditional_paras[6,1])^2
  return(ans)
}
