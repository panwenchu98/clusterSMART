#'@title General Hypothesis Testing
#'
#'@description
#'This is a function for doing general hypothesis testing with given solve_SMART output.
#'
#'@param result List output of solve_SMART function, containing details of parameter estimation and variance estiamtor.
#'@param aimed_test r*(q+p+1) matrix, each row corresponds to a linear combination of coefficients that is checked in hypothesis. q=4 is number of Adaptive Interventions, p is number of baseline covariates, the last row is constant term.
#'@param alpha The significance level (in inference we construct the (1-alpha) confidence interval).
#'@param use_t Logical. If TRUE uses t-distribution instead of normal when constructing the confidence interval and calculating the p-value.
#'
#'@return A matrix of test result.'
#'
#'
#'@examples
#'library(clusterSMART)
#'library(geepack)
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
#'                    p1=0.25,p2=0.55,N=500,Mmax=5,Mmin=2,seed=1234)
#'res=solve_SMART(Y=temp$data$Y,X=temp$data$X,cluster_id=temp$data$id,
#'                A1=temp$data$A1,R=temp$data$R,A2=temp$data$A2,
#'                aimed_comparison=matrix(c(1,0,0,1,1,0,1,0),2,4,byrow=TRUE),
#'                variance_structure=0,correlation_structure=1,ICC_lower_thresh=0,
#'                max_iter=10,convergence_thresh=1e-5,alpha=0.05,
#'                estimate_weight=FALSE,dof_adjustment=FALSE,use_t=FALSE,
#'                bias_correction=FALSE,verbose=2)
#'hypothesis_testing(res,aimed_test=matrix(c(1,-1,0,0,1,0,-5,
#'                                           1,0,0,0,0,-1,1),2,7,byrow=TRUE),
#'                   alpha=0.05,use_t=TRUE)
#'
#'@importFrom stats qnorm
#'@importFrom stats pnorm
#'@importFrom stats qt
#'@importFrom stats pt
#'
#'@export
#'
#'


hypothesis_testing=function(result,aimed_test,alpha=0.05,use_t=TRUE){
  q=4; p=nrow(result$summary_paras)-q
  if (is.null(dim(aimed_test))){
    aimed_test=matrix(aimed_test,nrow=1)
    r=1
  }else{
    r=nrow(aimed_test)
    if (!is.matrix(aimed_test)) aimed_test=as.matrix(aimed_test)
  }
  if (ncol(aimed_test)!=(p+q+1)) stop("Error: aimed_test should have q+p+1 columns")

  summary_test=matrix(0,r,7)
  colnames(summary_test)=c("Hypothesis_Estimand","Estimate","Std.Err","CI Lower","CI Higher","Test Score","p-value")
  summary_test=data.frame(summary_test)
  variable_names=result$summary_paras[,1]
  for (rr in 1:r){
    weight_beta=aimed_test[rr,1:(p+q)]

    H0=""
    t=weight_beta[1]
    if (t>0) {
      if (t==1){
        H0=variable_names[1]
      }else{
        H0=paste(round(t,digits=2)," * ",variable_names[1],sep="")
      }
    }
    if (t<0){
      if (t==-1){
        H0=paste("- ",variable_names[1],sep="")
      }else{
        H0=paste("- ",round(-t,digits=2)," * ",variable_names[1],sep="")
      }
    }

    for (i in 2:(p+q)){
      t=weight_beta[i]
      if (t>0) {
        if (t==1){
          H0=paste(H0," + ",variable_names[i],sep="")
        }else{
          H0=paste(H0," + ",round(t,digits=2)," * ",variable_names[i],sep="")
        }
      }
      if (t<0){
        if (t==-1){
          H0=paste(H0," - ",variable_names[i],sep="")
        }else{
          H0=paste(H0," - ",round(-t,digits=2)," * ",variable_names[i],sep="")
        }
      }
    }

    t=aimed_test[rr,p+q+1]
    if (t>0){H0=paste(H0," + ",round(t,digits=2),sep="")}
    if (t<0){H0=paste(H0," - ",round(-t,digits=2),sep="")}

    #H0=paste(H0," = 0",sep="")
    summary_test[rr,1]=H0

    summary_test[rr,2]=sum(weight_beta*as.numeric(result$summary_paras[,2]))+aimed_test[rr,p+q+1]
    summary_test[rr,3]=sqrt(t(matrix(weight_beta))%*%result$var_estimator%*%matrix(weight_beta))
    summary_test[rr,6]=summary_test[rr,2]/summary_test[rr,3]
    if (!use_t){
      summary_test[rr,4]=summary_test[rr,2]-qnorm(1-alpha/2)*summary_test[rr,3]
      summary_test[rr,5]=summary_test[rr,2]+qnorm(1-alpha/2)*summary_test[rr,3]
      summary_test[rr,7]=2*(1-pnorm(abs(summary_test[rr,4])))
    }else{
      summary_test[rr,4]=summary_test[rr,2]-qt(1-alpha/2,df=result$N-p-q)*summary_test[rr,3]
      summary_test[rr,5]=summary_test[rr,2]+qt(1-alpha/2,df=result$N-p-q)*summary_test[rr,3]
      summary_test[rr,7]=2*(1-pt(abs(summary_test[rr,4]),df=result$N-p-q))
    }
  }
  return(summary_test)
}
