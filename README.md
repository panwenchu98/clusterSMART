# clusterSMART
<!-- badges: start -->
[![R-CMD-check](https://github.com/panwenchu98/clusterSMART/actions/workflows/check-release.yaml/badge.svg)](https://github.com/panwenchu98/clusterSMART/actions/workflows/check-release.yaml)
[![test-coverage](https://github.com/panwenchu98/clusterSMART/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/panwenchu98/clusterSMART/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Overview  
This package includes a series of functions, listed below are the main functioning ones: 
   * parameter_adjust helps to adjust initial parameters to achieve certain effect size (signal-to-noise ratio) of main effect and covaraite effect.
   * generate_SMART generates clustered SMART data based on given parameters for simulation testing.
   * solve_SMART does all analysis for clustered SMART data, including effect estimate, working variance model estimate, as well as variance estimator.
   * hypothesis_testing does more generalized hypothesis testing using solve_SMART output.
   * data_expand duplicates data to make it capable of being analysed in off-the-shelf geepack.
----

## installation  
```{r}
  install.packages("devtools")
  library(devtools)
```
```{r}
  install_github("panwenchu98/clusterSMART",build_vignettes = T)
  library(clusterSMART)
```

----
## Generate cluster SMART data 
```{r}
library(stats)
library(MASS)
library(geepack)

conditional_paras_original=matrix(c(6,8,4,4,4,4,100,81,64,49,36,36,0.20,0.18,0.16,0.14,0.12,0.10),6,3)
p1=0.25; p2=0.55
effectsize_dtr=0.5; rhoadjust=1

effectsize_eta=c(0.2,0.5); sigma2_x=c(1,2)
#Adjusting Parameters
paras=parameter_adjust(conditional_paras_original=conditional_paras_original,p1=p1,p2=p2,aimed_comparison=c(1,0,0,1),
             effectsize_dtr=effectsize_dtr,effectsize_eta=effectsize_eta,rhoadjust[rhoadjust],sigma2_x=sigma2_x)
             
#Building SMART data
temp=generate_SMART(conditional_paras=paras$conditional_paras,
                    eta_x=paras$eta,sigma2_x=c(1,2),detail_x=c(1,1),
                    p1=0.25,p2=0.55,N=500,Mmax=5,Mmin=2,seed=1234)
```

## Analysing cluster SMART
```{r}
res=solve_SMART(Y=temp$data$Y,X=temp$data$X,cluster_id=temp$data$id,A1=temp$data$A1,R=temp$data$R,A2=temp$data$A2,
                aimed_comparison=matrix(c(1,0,0,1,1,0,1,0),2,4,byrow=T),
                variance_structure=0,correlation_structure=1,ICC_lower_thresh=0,
                max_iter=10,convergence_thresh=1e-5,alpha=0.05,
                estimate_weight=FALSE,dof_adjustment=FALSE,use_t=FALSE,bias_correction=FALSE,verbose=3)
```

## General Hypothesis Testing
```{r}
hypothesis_testing(res,aimed_test=matrix(c(1,-1,0,0,1,0,-5,1,0,0,0,0,-1,1),2,7,byrow=T),alpha=0.05,use_t=TRUE)
```

For further analysis and comparison with geepack, see attached vignettes.
