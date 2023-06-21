#'@title Generate clustered SMART data based on given parameters
#'
#'@description
#'This is a function for generating clustered SMART data based on given conditional (treatment-pathway level) mean, variance and intra-cluster correlation, as well as possible coefficient and variance of mean-zero normal baseline covariates.\cr
#'Default parameters gives a set of data with 100 clusters, each of same size 10, and the effect size of comparing AI(1,1) to AI(-1,-1) is 0.5, and has no baseline covaraites.
#'
#'@param conditional_paras A 6*3 matrix, each row contains the conditional (treatment-pathway level) parameters of mean, variance and intra-cluster correlation.
#'@param eta_x NULL, or a number or numeric vector of coefficients of baseline covariates.
#'@param sigma2_x NULL, or a positive number of numeric vector of the variance of baseline covariates.
#'@param detail_x NULL, or a number or vector with value in {0,1}, specifying each baseline covariate being cluster-level(0) or individual_level(1), default is cluster level.
#'@param p1,p2 The response rate of the first-stage treatment given the first randomization.
#'@param N Number of clusters.
#'@param Mmax,Mmin The largest and smallest size of each cluster.
#'@param seed The random seed.
#'
#'
#'@return A list of generated clustered SMART data containing:\tabular{ll}{
#'  \code{data}\tab A list of the data, each row represents an individual.
#'  \itemize{
#'    \item{Y}{ - The final outcome}
#'    \item{id}{ - The cluster id that individual belongs to}
#'    \item{L}{ - The treatment pathway of the cluster that individual belongs to}
#'    \item{A1}{ - The result of first cluster-level randomization}
#'    \item{A2}{ - The result of second cluster-level randomization}
#'    \item{R}{ - The cluster-level response status to the first-stage treatment}
#'    \item{X}{ - NULL, or a number or numeric vector of the baseline covaraites}
#'  }
#'  \cr
#'  \tab \cr
#'  \code{marginal_paras}\tab A 4*2 numeric matrix of the marginal mean and variance of each Adaptive Intervention.\cr
#'  \tab \cr
#'  \code{beta}\tab The coefficients of \eqn{\beta_0,\beta_1,\beta_2,\beta_3} and the baseline covaraites.\cr
#'  }
#'
#'@examples
#'library(clusterSMART)
#'conditional_paras_original=matrix(c(6,8,4,4,4,4,
#'                                    100,81,64,49,36,36,
#'                                    0.20,0.18,0.16,0.14,0.12,0.10),6,3)
#'p1=0.25; p2=0.55
#'effectsize_dtr=0.5; rhoadjust=1
#'effectsize_eta=c(0.2,0.5); sigma2_x=c(1,2)
#'paras=parameter_adjust(conditional_paras_original=conditional_paras_original,
#'                       p1=p1,p2=p2,aimed_comparison=c(1,0,0,1),
#'                       effectsize_dtr=effectsize_dtr,effectsize_eta=effectsize_eta,
#'                       rhoadjust=rhoadjust,sigma2_x=sigma2_x)
#'temp=generate_SMART(conditional_paras=paras$conditional_paras,
#'                    eta_x=paras$eta,sigma2_x=c(1,2),detail_x=c(1,1),
#'                    p1=0.25,p2=0.55,N=500,Mmax=10,Mmin=1,seed=1234)
#'
#'@importFrom stats rbinom
#'@importFrom stats rnorm
#'@importFrom stats runif
#'@importFrom MASS mvrnorm
#'@export
#'




#Default setting is N=100, cluster size constant=10, effectsize_DTR=0.5, effectsize_eta=0 (no cluster-level covariates)
#eta_x: an array of coefficients of covariates (if not imputed, assuming there is no cluster-level covariates)
#sigma2_x: an array of variances of each covariate, default is all 1
#detail_x: an array specifying each covariate being cluster-level(0) or individual_level(1), default is cluster level
#p1 and p2 are response rate corresponding to different first-stage treatments
#N is cluster number, Mmin and Mmax are lower and upper bound of cluster size (and the actual size is uniformly drawn between)
generate_SMART=function(conditional_paras=matrix(c(6,8,4,4,4,4,75.4461,61.1113,48.2855,36.9686,27.1606,27.1606,0.20,0.18,0.16,0.14,0.12,0.10),6,3),
                        eta_x=NULL,sigma2_x=NULL,detail_x=NULL,
                        p1=0.25,p2=0.55,N=100,Mmax=10,Mmin=10,seed=1234){
  if (Mmax<Mmin) stop("The maximum cluster size should not be smaller than the minimum")

  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)

  marginal_paras=conditional_to_marginal(p1,p2,conditional_paras)
  beta_margin=rep(0,4)
  beta_margin[1]=(marginal_paras[1,1]+marginal_paras[2,1]+marginal_paras[3,1]+marginal_paras[4,1])/4
  beta_margin[2]=(marginal_paras[1,1]+marginal_paras[2,1]-marginal_paras[3,1]-marginal_paras[4,1])/4
  beta_margin[3]=(marginal_paras[1,1]-marginal_paras[2,1]+marginal_paras[3,1]-marginal_paras[4,1])/4
  beta_margin[4]=(marginal_paras[1,1]-marginal_paras[2,1]-marginal_paras[3,1]+marginal_paras[4,1])/4
  q=4

  if (is.null(eta_x)){
    p=0
    beta=beta_margin
  }else{
    p=length(eta_x)
    beta=c(beta_margin,eta_x)
    if (is.null(sigma2_x)) sigma2_x=rep(1,p)
    if (is.null(detail_x)) detail_x=rep(0,p)
  }

  #generating the belonging pathway and corresponding AI for each cluster
  set.seed(seed)
  correct=0
  while (correct==0){
    a1=rbinom(N,size=1,prob=0.5)
    r=rep(0,N); r[a1==1]=rbinom(sum(a1),size=1,prob=p1); r[a1==0]=rbinom(N-sum(a1),size=1,prob=p2)
    a2=rep(0,N); a2[r==0]=2*rbinom(N-sum(r),size=1,prob=0.5)-1
    a1=2*a1-1
    label=1+(1-a1)*3/2+(3*a2-1)*a2/2
    correct=as.numeric(length(table(label))==6)
  }
  within=matrix(0,N,4)
  within[,1]=as.numeric((a1==1)&(a2>=0)); within[,2]=as.numeric((a1==1)&(a2<=0))
  within[,3]=as.numeric((a1==-1)&(a2>=0)); within[,4]=as.numeric((a1==-1)&(a2<=0))


  M=as.integer(runif(N,Mmin,Mmax+1-1e-8))

  #Building covariates
  if (p>0){
    x=list()
    for (i in 1:p){
      if (detail_x[i]==0){
        x[[i]]=matrix(rep(rnorm(N,mean=0,sd=sqrt(sigma2_x[i])),Mmax),N,Mmax)
      }else{
        x[[i]]=matrix(rnorm(N*Mmax,mean=0,sd=sqrt(sigma2_x[i])),N,Mmax)
      }
    }
  }

  #Building residuals
  epsilon=matrix(0,N,Mmax); mu=matrix(0,N,Mmax)
  for (i in 1:N){
    epsilon[i,]=mvrnorm(1,mu=rep(0,Mmax),Sigma=conditional_paras[label[i],2]*CS(conditional_paras[label[i],3],Mmax))
    mu[i,]=rep(conditional_paras[label[i],1],Mmax)
  }

  #merging all data to generate the outcome
  y=mu+epsilon
  if (p>0){
    for (i in 1:p){
      y=y+eta_x[i]*x[[i]]
    }
  }

  #changing into long format data
  result=matrix(0,sum(M),6+p)
  colnames(result)=c(1:(6+p))
  colnames(result)[1:6]=c("Y","id","L","A1","A2","R")
  if (p>1) colnames(result)[7:(6+p)]=paste("X",c(1:p),sep="")
  if (p==1) colnames(result)[7]="X"
  result=data.frame(result)
  count=0
  for (i in 1:N){
    for (j in 1:M[i]){
      count=count+1
      result[count,1]=y[i,j]
      result[count,2]=i
      result[count,3]=label[i]
      result[count,4]=ifelse(label[i] %in% c(1,2,3),1,-1)
      result[count,5]=ifelse(label[i] %in% c(1,4),NA,ifelse(label[i] %in% c(2,5),1,-1))
      result[count,6]=as.numeric(is.na(result[count,5]))
      if (p>0){
        for (k in 1:p) result[count,6+k]=x[[k]][i,j]
      }
    }
  }
  if (p>0){
    data=list(Y=result$Y,id=result$id,L=result$L,A1=result$A1,A2=result$A2,R=result$R,X=result[,7:(6+p)])
  }else{
    data=list(Y=result$Y,id=result$id,L=result$L,A1=result$A1,A2=result$A2,R=result$R,X=NULL)
  }

  return(list(data=data,marginal_paras=marginal_paras,beta=beta))
}
