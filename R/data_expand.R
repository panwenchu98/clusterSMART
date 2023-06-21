#'@title Expand Data to use in geeglm
#'
#'@description
#'This is a function for doing data expansion trick to clustered SMART data to make it available for using off-the-shelf geeglm function to do proper inference.\cr
#'Such analysis is only valid when using independence working variance structure, or when the cluster sizes are equal and use the customed variance structure for exchangeable.
#'
#'@param Y Length n numeric vector or n*1 matrix of the individual outcome.
#'@param X NULL, length n vector or n*p matrix for the baseline.
#'@param cluster_id n*1 numeric vector of the cluster id which the individual belongs to.
#'@param A1 Length n vector with value in {-1,1}, represents the result of the cluster-level randomization for the first-stage treatment.
#'@param R Length n vector with value in {0,1}, represents the cluster-level response status to the first-stage treatment.
#'@param A2 Length n vector with values in {-1,1,NA}, represents the result of the cluster-level randomization for the second-stage treatment.
#'@param seperate Logical. If TRUE then each response cluster is seperated into two clusters, if FALSE the response cluster is expanded to be twice large. Default is FALSE.
#'
#'@return A list of following two objects:\tabular{ll}{
#'  \code{data}\tab The expanded dataset.\cr
#'  \tab \cr
#'  \code{structure}\tab If cluster sizes are not equal or seperate is TRUE, NULL; if cluster sizes are equal, a userdefined correlation structure that can be used in the 'zcor' entry in geeglm. \cr
#'  }
#'
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
#'effectsize_eta=NULL; sigma2_x=NULL
#'paras=parameter_adjust(conditional_paras_original=conditional_paras_original,
#'                       p1=p1,p2=p2,aimed_comparison=c(1,0,0,1),
#'                       effectsize_dtr=effectsize_dtr,effectsize_eta=effectsize_eta,
#'                       rhoadjust=rhoadjust,sigma2_x=sigma2_x)
#'temp=generate_SMART(conditional_paras=paras$conditional_paras,
#'                   p1=0.25,p2=0.55,N=50,Mmax=3,Mmin=3,seed=1234)
#'
#'expanded=data_expand(Y=temp$data$Y,cluster_id=temp$data$id,
#'                     A1=temp$data$A1,R=temp$data$R,A2=temp$data$A2,seperate=FALSE)
#'model1=geeglm(Y~A1+A2+I(A1*A2),id=id,family=gaussian(),weight=W,data=expanded$data,
#'              corstr="userdefined",zcor=expanded$structure)
#'summary(model1)
#'
#'@importFrom geepack genZcor
#'
#'@export
#'
#'

data_expand=function(Y,X=NULL,cluster_id,A1,R,A2,seperate=FALSE){
  id=cluster_id
  N=length(unique(id))
  M=as.numeric(table(id))
  id_name=names(table(id))

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
  print("Warning:  The analysis using geeglm is only valid when all clusters have same size or when using indepdendence working variance structure")

  if (is.null(X)){
    data=cbind(Y,id,L,A1,R,A2)
  }else{
    data=cbind(Y,id,L,A1,R,A2,X)
  }
  data=data.frame(data)

  ndata=nrow(data)

  if (seperate){
    maxid=max(data$id)
    maxid=10^trunc(log(maxid,base=10)+1)
    for (i in 1:ndata){
      if (is.na(data$A2[i])){
        new_row=data[i,]
        data$A2[i]=1
        new_row$A2=-1
        new_row$id=new_row$id+maxid
        data=rbind(data,new_row)
      }
    }
    data$W=ifelse(data$L %in% c(1,4),2,4)
    return(list(data=data,structure=NULL))
  }

  for (i in 1:ndata){
    if (is.na(data$A2[i])){
      new_row=data[i,]
      data$A2[i]=1
      new_row$A2=-1
      new_row$id=new_row$id
      data=rbind(data,new_row)
    }
  }
  data$W=ifelse(data$L %in% c(1,4),2,4)
  data=data[order(data$id),]

  if (length(unique(M))==1){
    m=nrow(data)
    wave=rep(0,m)
    wave[1]=1
    for (i in 2:m){
      if (data$id[i]==data$id[i-1]){
        wave[i]=wave[i-1]+1
      }else{
        wave[i]=1
      }
    }

    M=M[1]
    zcor1=genZcor(clusz=table(data$id),waves=wave,corstrv=2*M)
    zcor.toep1=matrix(0,nrow(zcor1),1)
    count=0
    for (i in 1:(2*M-1)){
      for (j in (i+1):(2*M)){
        count=count+1
        if ((j<=M)|(i>M)) zcor.toep1[,1]=zcor.toep1[,1]+zcor1[,count]
      }
    }
    return(list(data=data,structure=zcor.toep1))
  }else{
    return(list(data=data,structure=NULL))
  }
}
