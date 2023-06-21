# clusterSMART
<!-- badges: start -->
[![R-CMD-check](https://github.com/panwenchu98/clusterSMART/actions/workflows/check-release.yaml/badge.svg)](https://github.com/panwenchu98/clusterSMART/actions/workflows/check-release.yaml)
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

```{r}
##Marginal Mean Model:  Y ~ a1 + a2 + I(a1 * a2) + X1 + X2 
##500 Number of Clusters, Minimum Cluster Size = 2 , Maximum Cluster Size = 5 
##Algorithm stops after  4  iterations.
##Summary of model Coefficients:
##                                                                     
##   Parameter Estimate Std.Err CI.Lower CI.Higher Z.Score Pr(>|z|)    
## (Intercept)  4.84318 0.15774  4.53401   5.15235    30.7   <2e-16 ***
##          a1  0.87258 0.15853  0.56186    1.1833     5.5 3.71e-08 ***
##          a2  0.48733 0.12033  0.25148   0.72318    4.05 5.12e-05 ***
##    I(a1*a2)  0.47307 0.12031  0.23726   0.70887    3.93 8.42e-05 ***
##          X1  1.45899 0.14119  1.18226   1.73571   10.33   <2e-16 ***
##          X2  2.52112 0.10144   2.3223   2.71994   24.85   <2e-16 ***
##---
##Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
##
##Working Variance Structure: Exchangeable
##Detailed Structure: sigma2 * EXCH(rho,M[i])
##
##Warning:
##The following marginal mean and inference are valid only when:
##    The model only contains baseline covariates and their interaction with A1,A2
##    The baseline covaraites should have mean zero
##
##Marginal Parameters:
##              Mean Variance       ICC
##AI(1,1)   6.676156 32.82146 0.1242524
##AI(1,-1)  4.755363 32.82146 0.1242524
##AI(-1,1)  3.984866 32.82146 0.1242524
##AI(-1,-1) 3.956334 32.82146 0.1242524
##
##Result of Proposed Comparison between Adaptive Interventions:
##                                                                           
##        Comparison Estimate Std.Err CI.Lower CI.Higher Z.Score Pr(>|z|)    
## AI(1,1)-AI(-1,-1)  2.71982 0.42699  1.88295    3.5567    6.37 1.89e-10 ***
##  AI(1,1)-AI(-1,1)  2.69129 0.43439   1.8399   3.54268     6.2 5.81e-10 ***
##---
##Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
##
##
##Variance-Covariance matrix of the estimates
##              (Intercept)           a1            a2     I(a1*a2)            X1            X2
##(Intercept)  0.0248832615 0.0089065023  0.0037041157 0.0029525102  0.0011169075 -0.0006440195
##a1           0.0089065023 0.0251332935  0.0029829652 0.0037828408  0.0007065177  0.0020113721
##a2           0.0037041157 0.0029829652  0.0144800623 0.0089576884  0.0015132126 -0.0001609414
##I(a1*a2)     0.0029525102 0.0037828408  0.0089576884 0.0144749962  0.0011831942  0.0002554717
##X1           0.0011169075 0.0007065177  0.0015132126 0.0011831942  0.0199344244 -0.0003183243
##X2          -0.0006440195 0.0020113721 -0.0001609414 0.0002554717 -0.0003183243  0.0102905667
```

## General Hypothesis Testing
```{r}
hypothesis_testing(res,aimed_test=matrix(c(1,-1,0,0,1,0,-5,1,0,0,0,0,-1,1),2,7,byrow=T),alpha=0.05,use_t=TRUE)
```

```{r}
##        Hypothesis_Estimand  Estimate   Std.Err    CI.Lower CI.Higher Test.Score     p.value
##1 (Intercept) - a1 + X1 - 5 0.4295864 0.2301277 -0.02256334 0.8817362   1.866731 0.982007700
##2      (Intercept) - X2 + 1 3.3220586 0.1909499  2.94688444 3.6972327  17.397540 0.003361658
```

#For further analysis and comparison with geepack, see attached vignettes.
