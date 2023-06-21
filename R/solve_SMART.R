#'@title Solve clustered SMART
#'
#'@description
#'This is a function for solving clustered SMART data and do inference.\cr
#'This function uses iterative algorithm to get the coefficient estimates for main effects and covaraites, and uses sandwich estimator to get proper inference.\cr
#'Also, this function is flexible to use different types of working variance model and implement different small-sample modification of variance estiamtors.
#'
#'@param Y Length n numeric vector or n*1 matrix of the individual outcome.
#'@param X NULL, length n vector or n*p matrix for the baseline.
#'@param cluster_id n*1 numeric vector of the cluster id which the individual belongs to.
#'@param A1 Length n vector with value in {-1,1}, represents the result of the cluster-level randomization for the first-stage treatment.
#'@param R Length n vector with value in {0,1}, represents the cluster-level response status to the first-stage treatment.
#'@param A2 Length n vector with values in {-1,1,NA}, represents the result of the cluster-level randomization for the second-stage treatment.
#'@param aimed_comparison NULL, r*4 matrix, each row has 2 0's and 2 1's, corresponding to the proposed comparison between AIs (default is all pairs of comparison).
#'@param variance_structure The marginal variance in the working model, 0 for same variance across Adaptive Intervention, 1 for different variance (default = 1).
#'@param correlation_structure The marginal intra-cluster correlation in the working model, 0 for independence, 1 for same ICC across AI, 2 for different ICC (default = 2).
#'@param ICC_lower_thresh A value in [-1,1] representing the set lower cutting threshold for ICC in working model, default is 0, setting to -1 means no cutting.
#'@param max_iter A natural integer representing the max allowed iteration in solving the effect estimates.
#'@param convergence_thresh A positive value representing the convergence threshold in the iterative algorithm of solving the effect estiamtes.
#'@param alpha The significance level (in inference we construct the (1-alpha) confidence interval).
#'@param estimate_weight Logical. If TRUE uses empirical estiamted inverse probability weight and adjust for the uncertainty introduced in the estimation in the variance estimator.
#'@param dof_adjustment Logical. If TRUE uses degree-of-freedom adjustment in the variance estiamtor.
#'@param use_t Logical. If TRUE uses t-distribution instead of normal when constructing the confidence interval and calculating the p-value.
#'@param bias_correction Logical. If TRUE uses the bias-corrected variance estimator for small-sample inference.
#'@param verbose Numeric value in {0,1,2,3} representing the printed output format. 0 for printing no output, 1 for coefficient estimation, 2 for coefficient estimation and working variance parameters, 3 for most detailed output including aimed comparison and variance estiamtor.
#'
#'
#'@return A list of result of analysis:\tabular{ll}{
#'  \code{summary_paras}\tab A list of the result of coefficient estimate.\cr
#'  \tab \cr
#'  \code{marginal_paras}\tab A 4*3 numeric matrix of the marginal mean,variance and ICC of each Adaptive Intervention.\cr
#'  \tab \cr
#'  \code{var_estimator}\tab The variance estimator of the coefficients.\cr
#'  \tab \cr
#'  \code{formula}\tab The formula of marginal mean regression model.\cr
#'  \tab \cr
#'  \code{comparison}\tab NULL or a list of the result of the aimed comparison between Adaptive Interventions.\cr
#'  \tab \cr
#'  \code{N}\tab Number of clusters.\cr
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
#'res=solve_SMART(Y=temp$data$Y,X=temp$data$X,cluster_id=temp$data$id,
#'                A1=temp$data$A1,R=temp$data$R,A2=temp$data$A2,
#'                aimed_comparison=matrix(c(1,0,0,1,1,0,1,0),2,4,byrow=TRUE),
#'                variance_structure=0,correlation_structure=1,ICC_lower_thresh=0,
#'                max_iter=10,convergence_thresh=1e-5,alpha=0.05,
#'                estimate_weight=FALSE,dof_adjustment=FALSE,use_t=FALSE,
#'                bias_correction=FALSE,verbose=2)
#'
#'@importFrom stats as.formula
#'@importFrom stats qnorm
#'@importFrom stats pnorm
#'@importFrom stats qt
#'@importFrom stats pt
#'@export
#'

solve_SMART=function(Y,X=NULL,cluster_id,A1,R,A2,
                     aimed_comparison=matrix(c(1,1,1,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,0,1,1),6,4),
                     variance_structure=1,correlation_structure=2,ICC_lower_thresh=0,
                     max_iter=10,convergence_thresh=1e-5,alpha=0.05,
                     estimate_weight=FALSE,dof_adjustment=TRUE,use_t=TRUE,bias_correction=TRUE,verbose=2){
  coeff=matrix(c(1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1),4,4)
  AI_names=c("AI(1,1)","AI(1,-1)","AI(-1,1)","AI(-1,-1)")

  #need to check imput validity
  if ((ICC_lower_thresh<(-1))|(ICC_lower_thresh>1)) stop("ICC_lower_thresh should be between -1 and 1")

  id=cluster_id
  N=length(unique(id))
  id_name=names(table(id))
  Mmax=max(table(id))
  M=as.numeric(table(id))

  if (all.equal(sort(unique(A1)),c(-1,1))!=TRUE) stop("A1 must be in -1 and 1")
  if (all.equal(sort(unique(R)),c(0,1))!=TRUE) stop("R must be in 0 and 1")
  if (length(unique(A2))!=3) stop("A2 must be in -1, 1 and NA")
  if (length(sort(unique(A2)))!=2) stop("A2 must be in -1, 1 and NA")
  if (all.equal(sort(unique(A2)),c(-1,1))!=TRUE) stop("A2 must be in -1, 1 and NA")
  if (length(unique(A2[which(R==1)]))>1) stop("The responders cannot be randomized in prototypical SMART")
  if (!is.na(unique(A2[which(R==1)]))) stop("The responders cannot be randomized in prototypical SMART")
  if (length(unique(A2[which(R==0)]))!=2) stop("The non-responders' A2 should be in -1 and 1")
  if (all.equal(sort(unique(A2[which(R==0)])),c(-1,1))!=TRUE) stop("The non-responders' A2 should be in -1 and 1")
  L=ifelse(A1==1,ifelse(R==1,1,ifelse(A2==1,2,3)),ifelse(R==1,4,ifelse(A2==1,5,6)))

  inside_list=list()
  for (i in 1:N) inside_list[[i]]=which(id==id_name[i])
  label=rep(0,N)
  for (i in 1:N){
    ll=L[inside_list[[i]]]
    if (length(unique(ll))!=1){
      stop("id is not nested in treatment pathway")
    }else{
      label[i]=unique(ll)
    }
  }
  if (length(unique(label))<6) stop("not every treatment pathway has observations, postitivity assumption violated")
  a1=ifelse(label>=4,-1,1)
  a2=ifelse(label %in% c(1,4),0,ifelse(label %in% c(2,5),1,-1))
  within=matrix(0,N,4)
  within[,1]=as.numeric((a1==1)&(a2>=0)); within[,2]=as.numeric((a1==1)&(a2<=0))
  within[,3]=as.numeric((a1==-1)&(a2>=0)); within[,4]=as.numeric((a1==-1)&(a2<=0))

  q=4
  if (is.null(X)){
    p=0
  }else{
    if (is.null(dim(X))){
      X=matrix(X,nrow=length(X),ncol=1)
      colnames(X)="X"
      p=1
    }else{
      p=ncol(X)
      if (!is.matrix(X)) X=as.matrix(X)
    }
  }

  if (estimate_weight==FALSE){
    count=as.numeric(table(label))
    p1_hat=count[1]/sum(count[1:3]); p2_hat=count[4]/sum(count[4:6])
    pa1_hat=1/2; p1a2_hat=1/2; p2a2_hat=1/2
    w=c(2,4,4,2,4,4)
  }else{
    count=as.numeric(table(label))
    pa1_hat=sum(count[1:3])/N
    p1_hat=count[1]/sum(count[1:3]); p2_hat=count[4]/sum(count[4:6])
    p1a2_hat=count[2]/(count[2]+count[3]); p2a2_hat=count[5]/(count[5]+count[6])
    w=c(1/pa1_hat,1/pa1_hat/p1a2_hat,1/pa1_hat/(1-p1a2_hat),1/(1-pa1_hat),1/(1-pa1_hat)/p2a2_hat,1/(1-pa1_hat)/(1-p2a2_hat))
  }


  #Use iterative method to calculate beta_hat and working covariance matrix (model based V_naive)
  sigma2_hat=rep(1,q); rho_hat=rep(0,q); beta_hat=rep(0,q+p); diff=1
  count_iter=0
  while ((diff>convergence_thresh)&(count_iter<max_iter)){
    count_iter=count_iter+1
    V_naive=matrix(0,q+p,q+p)
    weighty=matrix(0,q+p,1)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          V_naive=V_naive+w[label[i]]*t(Di)%*%solve(Vi)%*%Di
          weighty=weighty+w[label[i]]*t(Di)%*%solve(Vi)%*%yi
        }
      }
    }
    V_naive=solve(V_naive)
    beta_new=V_naive%*%weighty
    diff=sum(abs(beta_new-beta_hat))
    beta_hat=beta_new

    sigma_group=rep(0,q); weight_sigma=rep(0,q)
    rho_group=rep(0,q); weight_rho=rep(0,q)

    for (group in 1:q){
      for (i in 1:N){
        if (within[i,group]==1){
          mu_hat=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1) mu_hat=mu_hat+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mu_hat=mu_hat+X[inside_list[[i]],]*beta_hat[q+1]
          epsilon_hat=as.matrix(Y[inside_list[[i]]],M[i],1)-mu_hat
          sigma_group[group]=sigma_group[group]+w[label[i]]*sum(epsilon_hat^2)
          weight_sigma[group]=weight_sigma[group]+w[label[i]]*M[i]
          if (correlation_structure>=1){
            rho_group[group]=rho_group[group]+w[label[i]]*(sum(epsilon_hat)^2-sum(epsilon_hat^2))
            weight_rho[group]=weight_rho[group]+w[label[i]]*M[i]*(M[i]-1)}
        }
      }
    }
    if (variance_structure==0){
      sigma2_hat=rep(sum(sigma_group)/sum(weight_sigma),q)}
    else{
      sigma2_hat=sigma_group/weight_sigma}
    weight_rho=weight_rho*sigma2_hat
    if (correlation_structure==0){rho_hat=rep(0,q)}

    if (correlation_structure==1){rho_hat=rep(max(ICC_lower_thresh,sum(rho_group)/sum(weight_rho)),q)}
    if (correlation_structure==2){
      rho_hat=rho_group/weight_rho
      for (group in 1:q) rho_hat[group]=max(rho_hat[group],ICC_lower_thresh)
    }
  }


  if(estimate_weight){
    S_mat=matrix(c(1/p1_hat,-1/(1-p1_hat),-1/(1-p1_hat),0,0,0,
                   0,0,0,1/p2_hat,-1/(1-p2_hat),-1/(1-p2_hat),
                   1/pa1_hat,1/pa1_hat,1/pa1_hat,-1/(1-pa1_hat),-1/(1-pa1_hat),-1/(1-pa1_hat),
                   0,1/p1a2_hat,-1/(1-p1a2_hat),0,0,0,
                   0,0,0,0,1/p2a2_hat,-1/(1-p2a2_hat))
                 ,5,6,byrow=T)
    W_der_mat=matrix(c(0,0,0,0,0,0,
                       0,0,0,0,0,0,
                       -1/(pa1_hat^2),-1/(pa1_hat^2*p1a2_hat),-1/(pa1_hat^2*(1-p1a2_hat)),1/((1-pa1_hat)^2),1/((1-pa1_hat)^2*p2a2_hat),1/((1-pa1_hat)^2*(1-p2a2_hat)),
                       0,-1/(pa1_hat*p1a2_hat^2),1/(pa1_hat*(1-p1a2_hat)^2),0,0,0,
                       0,0,0,0,-1/((1-pa1_hat)*p2a2_hat^2),1/((1-pa1_hat)*(1-p2a2_hat)^2))
                     ,5,6,byrow=T)

    B_mat=matrix(0,5,5)
    for (i in 1:N){
      B_mat=B_mat+S_mat[,label[i]]%*%t(S_mat[,label[i]])
    }
  }


  if (!bias_correction){
    #Calculating Liang&Zeger Robust Sandwich Estimator
    meat_lz=matrix(0,p+q,p+q)
    if (estimate_weight) C_mat=matrix(0,p+q,5)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      Ui=matrix(0,p+q,1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          mui=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1) mui=mui+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mui=mui+X[inside_list[[i]],]*beta_hat[q+1]
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          Ui=Ui+w[label[i]]*t(Di)%*%solve(Vi)%*%(yi-mui)
          if (estimate_weight) C_mat=C_mat+t(Di)%*%solve(Vi)%*%(yi-mui)%*%t(W_der_mat[,label[i]])
        }
      }
      meat_lz=meat_lz+Ui%*%t(Ui)
    }

    if (estimate_weight){
      V_lz=V_naive%*%(meat_lz+C_mat%*%solve(B_mat)%*%t(C_mat))%*%V_naive}
    else{
      V_lz=V_naive%*%meat_lz%*%V_naive}
    V_result=V_lz
  }else{
    #Calculating Mancl&Derhoen Bias Corrected Estimator
    meat_md=matrix(0,p+q,p+q)
    if (estimate_weight) C_mat=matrix(0,p+q,5)
    for (i in 1:N){
      yi=as.matrix(Y[inside_list[[i]]],M[i],1)
      Ui=matrix(0,p+q,1)
      for (group in 1:q){
        if (within[i,group]==1){
          Di=matrix(rep(c(coeff[group,]),each=M[i]),M[i],q)
          if (p>0) Di=cbind(Di,matrix(X[inside_list[[i]],],M[i],p))
          mui=beta_hat[1]*coeff[group,1]+beta_hat[2]*coeff[group,2]+beta_hat[3]*coeff[group,3]+beta_hat[4]*coeff[group,4]
          if (p>1)mui=mui+X[inside_list[[i]],]%*%beta_hat[(q+1):(q+p),1]
          if (p==1) mui=mui+X[inside_list[[i]],]*beta_hat[q+1]
          Vi=sigma2_hat[group]*CS(rho_hat[group],M[i])
          Hi=Di%*%V_naive%*%t(Di)%*%solve(Vi)
          Ui=Ui+w[label[i]]*t(Di)%*%solve(Vi)%*%solve(diag(1,M[i],M[i])-Hi)%*%(yi-mui)
          if (estimate_weight) C_mat=C_mat+t(Di)%*%solve(Vi)%*%(yi-mui)%*%t(W_der_mat[,label[i]])
        }
      }
      meat_md=meat_md+Ui%*%t(Ui)
    }
    if (estimate_weight){
      V_md=V_naive%*%(meat_md+C_mat%*%solve(B_mat)%*%t(C_mat))%*%V_naive}
    else{
      V_md=V_naive%*%meat_md%*%V_naive}
    V_result=V_md
  }

  if (dof_adjustment) V_result=V_result*N/(N-p-q)





  #--------------------------------output session------------------------
  #output model structure
  formu="Y ~ a1 + a2 + I(a1 * a2)"
  if (p>0){
    for (i in 1:p) formu=paste(formu," + ",colnames(X)[i],sep="")
  }
  if (verbose>=2){
    cat("Marginal Mean Model: ",formu,"\n")
    cat(N,"Number of Clusters, Minimum Cluster Size =",min(M),", Maximum Cluster Size =",max(M),"\n")
    cat("Algorithm stops after ",count_iter," iterations.\n")
  }
  formu=as.formula(formu)

  #construct parameter output
  summary_paras=matrix(0,nrow=p+q,ncol=7)
  colnames(summary_paras)=c("Parameter","Estimate","Std.Err","CI Lower","CI Higher","Z Score","p-value")
  summary_paras=data.frame(summary_paras)
  summary_paras[1:4,1]=c("(Intercept)","a1","a2","I(a1*a2)")
  if (p>0) summary_paras[5:(4+p),1]=colnames(X)
  summary_paras[,2]=beta_hat
  summary_paras[,3]=sqrt(diag(V_result))
  summary_paras[,6]=summary_paras[,2]/summary_paras[,3]
  if (!use_t){
    summary_paras[,4]=beta_hat-qnorm(1-alpha/2)*sqrt(diag(V_result))
    summary_paras[,5]=beta_hat+qnorm(1-alpha/2)*sqrt(diag(V_result))
    summary_paras[,7]=2*(1-pnorm(abs(summary_paras[,6])))
  }else{
    colnames(summary_paras)[6]=c("T Score")
    summary_paras[,4]=beta_hat-qt(1-alpha/2,df=N-p-q)*sqrt(diag(V_result))
    summary_paras[,5]=beta_hat+qt(1-alpha/2,df=N-p-q)*sqrt(diag(V_result))
    summary_paras[,7]=2*(1-pt(abs(summary_paras[,6]),df=N-p-q))
  }

  #construct aimed comparison
  if (is.null(aimed_comparison)){
    r=0
  }else{
    if (is.null(dim(aimed_comparison))){
      aimed_comparison=matrix(aimed_comparison,nrow=1,ncol=4)
      r=1
    }else{
      r=nrow(aimed_comparison)
      if (!is.matrix(aimed_comparison)) aimed_comparison=as.matrix(aimed_comparison)
    }

    for (rr in 1:r){
      if (all.equal(aimed_comparison[order(aimed_comparison[rr,])],c(0,0,1,1))==FALSE){
        stop("Formula Error: aimed_comparison should be a length 4 vector with two 0's and two 1's")
      }
      num_1=which(aimed_comparison[rr,]==1)[1]; num_2=which(aimed_comparison[rr,]==1)[2]  #which two AIs are being compared
      weight_beta=coeff[num_1,]-coeff[num_2,]

      sum_test=sum(weight_beta*beta_hat[1:4])
      var_test=t(matrix(weight_beta))%*%V_result[1:4,1:4]%*%matrix(weight_beta)
      sd_test=sqrt(var_test)
      value_test=sum_test/sd_test
      if (!use_t){
        ci_low=sum_test-qnorm(1-alpha/2)*sd_test
        ci_high=sum_test+qnorm(1-alpha/2)*sd_test
        p_test=2*(1-pnorm(abs(value_test)))
      }else{
        ci_low=sum_test-qt(1-alpha/2,df=N-p-q)*sd_test
        ci_high=sum_test+qt(1-alpha/2,df=N-p-q)*sd_test
        p_test=2*(1-pt(abs(value_test),df=N-p-q))
      }
      addon=c(paste(AI_names[num_1],"-",AI_names[num_2],sep=""),sum_test,sd_test,ci_low,ci_high,value_test,p_test)
      summary_paras=rbind(summary_paras,addon)
    }
  }

  #formatting the parameter output
  output_paras=data.frame(cbind(summary_paras,rep(0,nrow(summary_paras))))
  output_paras[,2]=round(as.numeric(output_paras[,2]),digits=5)
  output_paras[,3]=round(as.numeric(output_paras[,3]),digits=5)
  output_paras[,4]=round(as.numeric(output_paras[,4]),digits=5)
  output_paras[,5]=round(as.numeric(output_paras[,5]),digits=5)
  output_paras[,6]=round(as.numeric(output_paras[,6]),digits=2)
  for (i in 1:nrow(summary_paras)){
    t=as.numeric(output_paras[i,7])
    output_paras[i,8]=ifelse(t>0.1," ",ifelse(t>0.05,".",ifelse(t>0.01,"*",ifelse(t>0.001,"**","***"))))
    output_paras[i,7]=ifelse(t<2e-16,"<2e-16",ifelse(t>0.001,as.character(round(t,digits=4)),formatC(t,format="e",digits=2)))
  }
  output_paras=rbind(colnames(output_paras),output_paras)
  output_paras[1,8]=""
  output_paras[1,7]=ifelse(!use_t,"Pr(>|z|)","Pr(>|t|)")
  colnames(output_paras)=NULL
  rownames(output_paras)=NULL

  #outputting the model coefficients
  if (verbose>=1){
    cat("Summary of model Coefficients:\n")
    print(output_paras[1:(p+q+1),],row.names=FALSE)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }



  marginal_paras=matrix(0,4,3)
  colnames(marginal_paras)=c("Mean","Variance","ICC")
  rownames(marginal_paras)=AI_names
  marginal_paras[,1]=coeff%*%matrix(as.numeric(summary_paras[1:4,2]),4,1)
  marginal_paras[,2]=sigma2_hat
  marginal_paras[,3]=rho_hat

  if(verbose>=2){
    cat("\n")
    cat("Working Variance Structure: ")
    if (correlation_structure==0){cat("Independence")}else{cat("Exchangeable")}
    cat("\nDetailed Structure: ")
    if (variance_structure==1){cat("sigma2_{a1,a2} * ")} else{cat("sigma2 * ")}
    if (correlation_structure==0) cat("I_{M[i]}")
    if (correlation_structure==1) cat("EXCH(rho,M[i])")
    if (correlation_structure==2) cat("EXCH(rho_{a1,a2},M[i])")
    cat("\n")
  }

  if (verbose==2){
    cat("Marginal Parameters:\n")
    if ((variance_structure==0)&(correlation_structure<=1)){
      if (correlation_structure==0) cat("sigma2 =",marginal_paras[1,2])
      if (correlation_structure==1) cat("sigma2 =",marginal_paras[1,2]," , rho =",marginal_paras[1,3])
    }else{
      if (correlation_structure==0) print(marginal_paras[,2])
      if (correlation_structure>0) print(marginal_paras[,-1])
    }
    cat("\n")
  }


  if ((verbose>=3)&(p+r>0)){
    cat("\nWarning:\n")
    cat("The following marginal mean and inference are valid only when:\n")
    cat("    The model only contains baseline covariates and their interaction with A1,A2\n")
    cat("    The baseline covaraites should have mean zero\n")
  }

  if (verbose>=3){
    cat("\nMarginal Parameters:\n")
    if (correlation_structure==0) print(marginal_paras[,-3]) else print(marginal_paras)
    cat("\n")
  }

  if (r>0){
    summary_comparison=summary_paras[(p+q+1):(p+q+r),]
    colnames(summary_comparison)[1]="Comparison"
    rownames(summary_comparison)=c()
    summary_paras=summary_paras[1:(p+q),]

    if (verbose>=3){
      output_comparison=output_paras[c(1,(p+q+2):(p+q+r+1)),]
      output_comparison[1,1]="Comparison"
      if (r==1){
        cat("Result of Proposed Comparison between ",AI_names[num_1]," and ",AI_names[num_2],"\n")
      }else{
        cat("Result of Proposed Comparison between Adaptive Interventions:\n")
      }
      print(output_comparison,row.names=FALSE)
      cat("---\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      cat("\n")
    }
  }




  rownames(V_result)=summary_paras[1:(p+q),1]
  colnames(V_result)=summary_paras[1:(p+q),1]

  if (verbose>=3){
    cat("\n")
    cat("Variance-Covariance matrix of the estimates\n")
    print(V_result)
  }

  for (i in 2:7){
    summary_paras[,i]=as.numeric(summary_paras[,i])
  }

  if (!is.null(aimed_comparison)){
    for (i in 2:7) {summary_comparison[,i]=as.numeric(summary_comparison[,i])}
    return(list(summary_paras=summary_paras,marginal_paras=marginal_paras,var_estimator=V_result,
                formula=formu,comparison=summary_comparison,N=N))
  }else{
    return(list(summary_paras=summary_paras,marginal_paras=marginal_paras,var_estimator=V_result,
                formula=formu,comparison=NULL,N=N))
  }
}
