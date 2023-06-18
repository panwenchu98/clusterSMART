#'@title Adjust parameters to satisfy certain requirements
#'
#'@description
#'This is a function for changing the given initial conditional (treatment-pathway level) parameters based on the required effect size of aimed comparison and effect size of covariates.
#'
#'@param conditional_paras_original A 6*3 matrix, each row contains of the given initial conditional (treatment-pathway level) parameters of mean, variance and intra-cluster correlation.
#'@param p1,p2 The response rate of the first-stage treatment given the first randomization
#'@param aimed_comparison A length 4 vector containing two 0s and two 1s, representing which two Adaptive Interventions we want to fix the effect size (default is c(1,0,0,1) which represents comparing AI(1,1) vs AI(-1,-1))
#'@param effectsize_dtr The desired effect size of the aimed comparison
#'@param effectsize_eta The desired effect size of each baseline covariates (can be NULL, value or vectors in [0,1]) (default is NULL)
#'@param rhoadjust The scale modification parameter that is multiplied onto the baseline correlation
#'@param sigma2_x The variance of desired covaraites. (can be NULL, positive value or vectors of the same length as effectsize_eta) (default is NULL)
#'
#'@return A list of adjusted paramters containing:\tabular{ll}{
#'  \code{beta_margin}\tab A 4*1 numeric vector of the parameters \eqn{\beta_0,\beta_1,\beta_2,\beta_3}\cr
#'  \tab \cr
#'  \code{conditional_paras}\tab A 6*3 numeric matrix of the mean, variance and ICC of each treatment trajectory \cr
#'  \tab \cr
#'  \code{marginal_paras}\tab A 4*2 numeric matrix of the marginal mean and variance of each Adaptive Intervention\cr
#'  \tab \cr
#'  \code{eta}\tab NULL, or a value or a numeric vector of the coefficients of each baseline covariate\cr
#'  }
#'
#'@examples
#'library(clusterSMART)
#'conditional_paras_original=matrix(c(6,8,4,4,4,4,100,81,64,49,36,36,0.20,0.18,0.16,0.14,0.12,0.10),6,3)
#'p1=0.25; p2=0.55
#'effectsize_dtr=0.5; rhoadjust=1
#'effectsize_eta=c(0.2,0.5); sigma2_x=c(1,2)
#'paras=parameter_adjust(conditional_paras_original=conditional_paras_original,p1=p1,p2=p2,aimed_comparison=c(1,0,0,1),effectsize_dtr=effectsize_dtr,effectsize_eta=effectsize_eta,rhoadjust[rhoadjust],sigma2_x=sigma2_x)
#'
#'@export
#'



#Given original conditional parameters and desired effect size, adjust the parameters to achieve the aim
#  the adjustment is made by changing the coefficients of covariates, as well as modifying the variance and correlation
#conditional_paras_original is a 6*3 matrix containing the starting mean, variance and ICC for each pathway
#  the conditional mean effect, and the proportion of variance and ICC between pathways will be kept with the adjustment
#If for no covariates, then just not inputing effectsize_eta
#aimed_comparison is for the aimed comparison of Adaptive Interventions (default is comparing AI(1,1) vs AI(-1,-1))
#  format of aimed_comparison: should be a vector length in 4, containing 2 0's and 2 1's
parameter_adjust=function(conditional_paras_original,p1,p2,aimed_comparison=c(1,0,0,1),
                          effectsize_dtr,effectsize_eta=NULL,rhoadjust,sigma2_x=NULL){
  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)  #the matrix that changes coefficients to marginal mean
  q=4
  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras_original)
  beta_margin=solve(coeff) %*% as.matrix(marginal_paras[,1],4,1)

  if (!(all.equal(aimed_comparison[order(aimed_comparison)],c(0,0,1,1))==TRUE)){
    stop("Error: aimed_comparison should be a length 4 vector with two 0's and two 1's")
  }
  num_1=which(aimed_comparison==1)[1]; num_2=which(aimed_comparison==1)[2]  #which two AIs are being compared

  weight_beta=coeff[num_1,]-coeff[num_2,]
  sum_beta=sum(weight_beta*beta_margin)

  if (is.null(effectsize_eta)){
    p=0
    eta=NULL
  }else{
    if (sum(effectsize_eta^2)>=1) stop("Error: total effect size of covariates too large, their sum of squares exceeding 1")
    p=length(effectsize_eta)
    if (is.null(sigma2_x)) sigma2_x=rep(1,p)
    if (length(sigma2_x)!=p) stop("Error: different length of effect size and variance of covariates")
    eta=effectsize_eta/effectsize_dtr*sum_beta/sqrt(sigma2_x)
  }
  sum_varx=sum(eta^2*sigma2_x)

  coef1=2*((sum_beta/effectsize_dtr)^2-sum_varx)

  coeff_new=matrix(0,4,2)
  #the within-pathway part that can be duplicated by variance_adjust
  coeff_new[1,1]=conditional_paras_original[1,2]*p1+conditional_paras_original[2,2]*(1-p1)
  coeff_new[2,1]=conditional_paras_original[1,2]*p1+conditional_paras_original[3,2]*(1-p1)
  coeff_new[3,1]=conditional_paras_original[4,2]*p2+conditional_paras_original[5,2]*(1-p2)
  coeff_new[4,1]=conditional_paras_original[4,2]*p2+conditional_paras_original[6,2]*(1-p2)
  #the between-pathway part that cannot be duplicated
  coeff_new[1,2]=p1*(1-p1)*(conditional_paras_original[1,1]-conditional_paras_original[2,1])^2
  coeff_new[2,2]=p1*(1-p1)*(conditional_paras_original[1,1]-conditional_paras_original[3,1])^2
  coeff_new[3,2]=p2*(1-p2)*(conditional_paras_original[4,1]-conditional_paras_original[5,1])^2
  coeff_new[4,2]=p2*(1-p2)*(conditional_paras_original[4,1]-conditional_paras_original[6,1])^2

  a=t(as.matrix(aimed_comparison)) %*% coeff_new
  varianceadjust=(coef1-a[2])/a[1]
  if (varianceadjust<0) stop("Error: total effect size too large, consider reducing effect of covariates or difference between pathway mean effects")

  conditional_paras=conditional_paras_original
  conditional_paras[,2]=conditional_paras[,2]*varianceadjust
  conditional_paras[,3]=conditional_paras[,3]*rhoadjust
  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras)

  result=list(beta_margin=beta_margin,conditional_paras=conditional_paras,marginal_paras=marginal_paras,eta=eta)
  return(result)
}
