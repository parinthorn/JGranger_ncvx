###############################################
# BootSE.R:  Functions for residual bootstrap #
###############################################

# (!!The functions should be in a folder named "Functions" since they are sources in this R-script)
# Ensure that the following lines work:
# source('Functions/MultiClass_VAR.R')
# source('Functions/BootSE.R')
# source('Functions/FunctionsVisualTools.R')


bootSE<-function(FITOBJECT,CFIT=NULL,J=NULL,K=NULL,P=NULL,n=NULL,maxit.both=10,maxit.beta=10,maxit.omega=10,mu=1e-4,
                 tol.both=10^-5,tol.beta=0.01,tol.omega=0.01,Nboot=1000,calculate.C=F,parallell=F,type="AdLasso",
                 groupindex=NULL, lambda_weights=F, Clusterinfo=NULL){
  
  #### Function for residual bootstrap with option for parallel computing ####  
  
  #########
  # INPUT #
  #########
  # FITOBJECT : output of the function  "MultiClass_VAR"
  # CFIT : output of the function "buildCmatrix" 
  # J : number of time series in each class
  # K : number of classes
  # P : order of the VAR
  # n : time series length
  # maxit.both : maximum iterations of the Multiclass VAR 
  # maxit.beta : maximum iterations  SPG algorithm
  # maxit.omega : maximum iterations JGL algorithm
  # mu : smoothness parameter
  # tol.both : tolerance  Multiclass VAR 
  # tol.beta : tolerance SPG algorithm
  # tol.omega : tolerance JGL algorithm
  # NBoot : number of bootstrap replicates
  # calculate.C : logical for calculating the C matrix. 
  # parallell : logical for residual bootstrap procedure carried out using parallel computing. Default  is NO parallel.
  # type: "Lasso", "AdLasso", "GrLasso" for Lasso, Adapative Lasso and Group Lasso respectively
  # groupindex : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda_weights : logical. If TRUE, then do the weighting of fusion parameter based on Clusterinfo. Default is FALSE.
  # Clusterinfo : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.
  # seed: fix seed of pseudo random number generator
  
  ##########
  # OUTPUT #
  ##########
  # BetaSign : significant coefficients based on standard errors computed with the residual bootstrap
  # Bootsd: bootstrap standard errors
  # BetaBoot: bootstrap replicates of the coefficients
  # Beta: estimated coefficients
  
  #### Start Code ####
  
  # Preliminaries: Elements from the MultiClass_VAR object
  Betahat_array<-FITOBJECT$Beta_array
  Betahat_vec<-FITOBJECT$Beta_new
  X<-FITOBJECT$X_array
  U<-FITOBJECT$E_array
  lam1.opt<-FITOBJECT$lambda1_opt
  lam2.opt<-FITOBJECT$lambda2_opt
  gam1.opt<-FITOBJECT$gamma1_opt
  gam2.opt<-FITOBJECT$gamma2_opt
  lamridge<-FITOBJECT$lambda_RIDGE
  J<-FITOBJECT$J
  K<-FITOBJECT$K
  P<-FITOBJECT$P
  n<-FITOBJECT$n
  
  if (is.null(CFIT)){
    CFIT<-FITOBJECT$CFIT
  }
  
  
  # 1. Center residuals
  Ucent<-apply(U,c(2,3),scale,scale=F,center=T) #U needs to be an array of dimension (n,J,K)
  
  BetaBoot<-matrix(NA,ncol=length(Betahat_vec),nrow=Nboot)
  
  if(parallell==T){# Parallel Computing
    
    checkpackage<-function(U){
      if((U %in% rownames(installed.packages()))==F){
        install.packages(U)
      }
    }
    packagelist<-list("foreach","doSNOW","doRNG")
    lapply(packagelist,checkpackage)
    # Load packages
    suppressMessages(suppressWarnings(packages <- lapply(packagelist, FUN = function(x) {
      library(x, character.only = TRUE)
    })))
    
    # Simulation with parallel computing
    ncl<-parallel:::detectCores()-2 #check number of cores, use two less than total number of cores
    R<-Nboot # number of simulations
    simestcl <- makeCluster(ncl) 
    registerDoSNOW(simestcl)
    
    
    BOOTFITS<- foreach(isim = 1:R, .packages=c('MASS', 'magic', 'rARPACK', 'Matrix', 'JGL', 'plyr')) %dorng% {
      
      # Additional functions to be sourced (!!The functions should be in a folder named "Functions")
      source('Functions/MultiClass_VAR.R')
      source('Functions/BootSE.R')
      source('Functions/FunctionsVisualTools.R')
      
      # if(!is.null(seed)){
      #   set.seed(seed + isim)
      # }

      BOOTFIT<-BootAUX(nboot=isim,Betahat_array=Betahat_array,Betahat_vec=Betahat_vec,X.data=X,lam1.opt=lam1.opt,lam2.opt=lam2.opt,gam1.opt=gam1.opt,gam2.opt=gam2.opt,
                       lamridge=lamridge,CFIT=CFIT,J=J,K=K,P=P,n=n,maxit.both=maxit.both,maxit.beta=maxit.beta,maxit.omega=maxit.omega,mu=mu,
                       tol.both=tol.both,tol.beta=tol.beta,tol.omega=tol.omega,Ucent=Ucent,calculate.C=calculate.C,type=type,groupindex=groupindex,
                       lambda_weights=lambda_weights, Clusterinfo=Clusterinfo)
      
      list("BOOTFIT"=BOOTFIT)
    }
    stopCluster(simestcl)
  }else{
    bootseq<-1:Nboot
    BOOTFITS<-lapply(bootseq,FUN=BootAUX,Betahat_array=Betahat_array,Betahat_vec=Betahat_vec,X.data=X,lam1.opt=lam1.opt,lam2.opt=lam2.opt,gam1.opt=gam1.opt,gam2.opt=gam2.opt,
                     lamridge=lamridge,CFIT=CFIT,J=J,K=K,P=P,n=n,maxit.both=maxit.both,maxit.beta=maxit.beta,maxit.omega=maxit.omega,mu=mu,
                     tol.both=tol.both,tol.beta=tol.beta,tol.omega=tol.omega,Ucent=Ucent,calculate.C=calculate.C,type=type,groupindex=groupindex,
                     lambda_weights=lambda_weights, Clusterinfo=Clusterinfo)
  }
  
  
  BetaBoot<-matrix(unlist(BOOTFITS),byrow=T,nrow=Nboot)
  
  # Compute Standard errors
  Bootsd<-apply(BetaBoot,2,sd)
  coef.index<-1:length(Bootsd)
  FLAGnorm<-unlist(lapply(coef.index,tstatnormal,estimate=Betahat_vec,sd=Bootsd))
  BetaSign<-rep(0,length(Betahat_vec))
  BetaSign[FLAGnorm]<-Betahat_vec[FLAGnorm]
  
  out<-list("BetaSign"=BetaSign,"Bootsd"=Bootsd,"BetaBoot"=BetaBoot, "Beta"=FITOBJECT$Beta_new)
}

BootAUX<-function(nboot, Betahat_array, Betahat_vec, X.data,
                  lam1.opt, lam2.opt, gam1.opt, gam2.opt, lamridge, CFIT,
                  J, K, P, n, maxit.both, maxit.beta, maxit.omega, mu,
                  tol.both, tol.beta, tol.omega, Ucent, calculate.C, 
                  criterion="BIC", type=NULL, groupindex=NULL, 
                  lambda_weights=NULL, Clusterinfo=NULL, max.it.test=100){
  #### Auxiliary Function: Perform residual bootstrap ####
  
  #########
  # INPUT #
  #########
  # nboot : bootstrap run
  # Betahat_array : array of estimated beta coefficients 
  # Betahat_vec : vector of estimated beta coefficients 
  # X.data : array of predictors 
  # lam1.opt : optimal lambda 1 (lasso)
  # lam2.opt : optimal lambda 2 (fused lasso)
  # gam1.opt : optimal gamma 1 (lasso)
  # gam2.opt : optimal gamma 2 (fused lasso)
  # lamridge : array (1,1,K) containing the optimal regularization parameter for the ridge, for each class
  # CFIT : output the function buildCmatrix.R
  # J : number of time series in each class
  # K : number of classes (K>1)
  # P : order of VAR
  # n : time series length
  # maxit.both : maximum iterations of the Multiclass VAR 
  # maxit.beta : maximum iterations  SPG algorithm
  # maxit.omega : maximum iterations JGL algorithm
  # mu : smoothness parameter
  # tol.both : tolerance  Multiclass VAR 
  # tol.beta : tolerance SPG algorithm
  # tol.omega : tolerance JGL algorithm
  # Ucent : centered residuals
  # calculate.C : logical for calculating the C matrix. 
  # criterion : criterion for the selection of the regularization parameters. Default is "BIC", otherwise set "AICc".
  # type: "Lasso", "AdLasso", "GrLasso" for Lasso, Adapative Lasso and Group Lasso
  # groupindex : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda_weights : logical. If TRUE, then do the weighting of fusion parameter based on Clusterinfo. Default is FALSE.
  # Clusterinfo : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.
  # max.it.test : maximum iteration for the multivariate white noise test. Default is 100.
  
  ##########
  # OUTPUT #
  ##########
  # BetaBoot : bootstrap beta coefficients 
  
  #### Start Code ####
  
  # 1. Draw residuals with replacement 
  flag<-1
  it<-1
  while(length(flag)!=0 & (it<max.it.test)){ #(check for multivariate white noise)
    
    n.index<-sample(1:n,n,replace=T)
    Uboot<-Ucent[n.index,,]
    
    pvalueMWN<-matrix(NA,ncol=dim(Uboot)[3],nrow=1)
    for(i.class in 1:dim(Uboot)[3]){
      MWNtest<-MultivariateWhiteNoise(Uboot[,,i.class], lag = floor(sqrt(dim(Uboot)[1])), adj = 0)
      pvalueMWN[1,i.class]<-MWNtest[nrow(MWNtest),4]
    }
    flag<-which(pvalueMWN<0.05)
    it<-it+1
  }
  
  
  # 2. Compute bootstrap time series
  BootDGP<-MultiClass_DG(Beta_DG=Betahat_array, Beta_vec_DG=Betahat_vec, X_DG=X.data, E_DG=Uboot,calculate.C=calculate.C)
  
  # 3. Multi-class fit with bootstrap time series, 
  n<-length(BootDGP$Y_input)/(J*K)
  BOOTFIT<-MultiClass_VAR(Data=NULL, calculate.C=F,X0=BootDGP$X0_input,XX_prod=BootDGP$XX_input,XY_prod=BootDGP$XY_input,Y=BootDGP$Y_input, X=BootDGP$X_input, C=CFIT$C, J=J, K=K, P=P, n=n, 
                          CNorm=CFIT$CNorm, maxit.both=maxit.both, maxit.beta=maxit.beta,maxit.omega=maxit.omega,b_init= array(0,c(ncol(BootDGP$X_input),1)), mu=mu, 
                          tol.both=tol.both, tol.beta=tol.beta,tol.omega=tol.omega, lambda1_OPT=lam1.opt, lambda2_OPT=lam2.opt, gamma1_OPT=gam1.opt, gamma2_OPT=gam2.opt,
                          lambdaRIDGE_OPT=lamridge,
                          criterion=criterion,type=type,group=groupindex,
                          lambda_weights=lambda_weights, Clusterinfo=Clusterinfo,
                          lambda1_min=NULL, lambda1_max=NULL, lambda1_steps=NULL, 
                          lambda2_min=NULL, lambda2_max=NULL, lambda2_steps=NULL,
                          gamma1_min=NULL, gamma1_max=NULL, gamma1_steps=NULL, 
                          gamma2_min=NULL, gamma2_max=NULL, gamma2_steps=NULL,
                          lambdaRIDGE_min=NULL, lambdaRIDGE_max=NULL, lambdaRIDGE_steps=NULL)
  
  
  BetaBoot<-BOOTFIT$Beta_new
  return(BetaBoot)
  
}

MultiClass_DG<-function(Beta_DG, Beta_vec_DG, X_DG, E_DG,calculate.C=FALSE){
  #### Function to generate data given the output of the function MultiClass_VAR ####
  
  #########
  # INPUT #
  #########
  # Beta_DG : array of autoregressive coefficients from the Multiclass_VAR
  # Beta_vec_DG : vector of autoregressive coefficients from the Multiclass_VAR
  # X_DG : array of responses from the Multiclass_VAR
  # E_DG : array of error from the  Multiclass_VAR
  # calculate.C : logical for calculating the C matrix. Default is FALSE.
  
  ##########
  # OUTPUT #
  ##########
  # Y_DG : NxJXK array generated responses
  # Y_input : NJKx1 vector of responses
  # X_input : NJKxP(J^2) matrix of inputs
  # X0_input : NxJPxK array of inputs
  # XX_input : product X'X
  # XY_input : product X'Y
  # C_input : C matrix for SPG algorithm
  # CNorm_input : norm of C for the SPG algorithm
  
  #### Start Code ####
  # Dimensions
  n_DG<-dim(E_DG)[1]
  J_DG<-dim(E_DG)[2]
  K_DG<-dim(E_DG)[3]
  P_DG<-dim(X_DG)[2]/J_DG
  
  # Data generation
  Y_DG<-array(NA, c(n_DG,J_DG,K_DG))
  
  for (i.class in 1:K_DG){
    BetaP_DG<-array(Beta_DG[,,i.class], c(J_DG,J_DG,P_DG)) # Array (J,J,P) of Beta coefficients at different lags p=1,..P
    
    # First time observation  
    Y_DG[1,,i.class]<-E_DG[1,,i.class]  # Y_1=E_1
    
    # Generate data until P+1 (included)
    for (i.n in 1:P_DG){
      X_P<-array(NA,c(J_DG,1,i.n)) # expanding inputs depending on the lag p=1,..,P
      arg_P<-array(NA, c(J_DG,1,i.n))
      for (ip in 1:i.n){
        X_P[,,ip]<-Y_DG[i.n-ip+1,,i.class]
        arg_P[,,ip]<-BetaP_DG[,,ip]%*%X_P[,,ip]
      }
      arg_P_LIST<-alply(arg_P, 3, .dims=TRUE)          # Sum across lags p=1,..P
      sum_DG<-Reduce('+', arg_P_LIST)
      Y_DG[(i.n+1),,i.class]<-sum_DG + E_DG[(i.n+1),,i.class]
    }
    
    # Generate data from P+2 (included)
    for (i.n in (P_DG+1):(n_DG-1)){
      X_P<-array(NA,c(J_DG,1,P_DG))
      arg_P<-array(NA, c(J_DG,1,P_DG))
      for (ip in 1:P_DG){
        X_P[,,ip]<-Y_DG[i.n-ip+1,,i.class]
        arg_P[,,ip]<-BetaP_DG[,,ip]%*%X_P[,,ip]
      }
      arg_P_LIST<-alply(arg_P, 3, .dims=TRUE)          # Sum across lags p=1,..P
      sum_DG<-Reduce('+', arg_P_LIST)
      Y_DG[(i.n+1),,i.class]<-sum_DG + E_DG[(i.n+1),,i.class]
    }
  }  
  
  # Generate the new input for the MultiClass_VAR based on the generated Y_DG
  inputs_multiclass<-Inputs_MultiClass_VAR(Data=Y_DG, P=P_DG,calculate.C=calculate.C)
  
  Y_input<-inputs_multiclass$Y
  X_input<-inputs_multiclass$X
  X0_input<-inputs_multiclass$X0
  XX_input<-inputs_multiclass$XX
  XY_input<-inputs_multiclass$XY
  C_input<-inputs_multiclass$C
  CNorm_input<-inputs_multiclass$CNorm
  
  
  MultiClass_DG<-list(Y_DG=Y_DG, Y_input=Y_input, 
                      X_input=X_input, X0_input=X0_input,
                      XX_input=XX_input, XY_input=XY_input,
                      C_input=C_input, CNorm_input=CNorm_input)  
  
  
}

tstatnormal<-function(index,estimate,sd){
  #### Auxiliary Function: Check the significance of the beta coefficients ####
  if(sd[index]==0){
    tstat<-0
  }else{
    tstat<-estimate[index]/sd[index]
  }
  
  flag<-abs(tstat)>qnorm(0.975) # return significant ones
}

MultivariateWhiteNoise<-function (x, lag = 24, adj = 0) {
  #### MTS PACKAGE: Multivariate Portmanteau from package MTS as in Tsay (2014), p72. ####
  # Input and Output: see documentation of function "mq" in package MTS
  
  if (!is.matrix(x)) 
    x = as.matrix(x)
  nr = nrow(x)
  nc = ncol(x)
  g0 = var(x)
  ginv = solve(g0)
  qm = 0
  QM = NULL
  df = 0
  for (i in 1:lag) {
    x1 = x[(i + 1):nr, ]
    x2 = x[1:(nr - i), ]
    g = cov(x1, x2)
    g = g * (nr - i - 1)/(nr - 1)
    h = t(g) %*% ginv %*% g %*% ginv
    qm = qm + nr * nr * sum(diag(h))/(nr - i)
    df = df + nc * nc
    dff = df - adj
    mindeg = nc^2 - 1
    pv = 1
    if (dff > mindeg) 
      pv = 1 - pchisq(qm, dff)
    QM = rbind(QM, c(i, qm, dff, pv))
  }
  pvs = QM[, 4]
  dimnames(QM) = list(names(pvs), c("  m  ", "    Q(m) ", "   df  ", 
                                    " p-value"))
  
  return(QM)
}

bootSE_GeneralOmega<-function(FITOBJECT,CFIT=NULL,J=NULL,K=NULL,P=NULL,n=NULL,maxit.both=10,maxit.beta=10,maxit.omega=10,mu=1e-4,
                              tol.both=10^-5,tol.beta=0.01,tol.omega=0.01,Nboot=1000,calculate.C=F,parallell=F,type="AdLasso",
                              groupindex=NULL, lambda_weights=F, Clusterinfo=NULL){
  
  #### Function for residual bootstrap when allowing for cross-error correlations ####  
  
  #########
  # INPUT #
  #########
  # FITOBJECT : output of the function  MultiClass_VAR_generalOmega.R
  # CFIT : output the function buildCmatrix.R
  # J : number of variables in each class
  # K : number of classes (K>1)
  # P : lag order of VAR
  # n : time series length
  # maxit.both : maximum iteration of the Multiclass VAR for joint estimation of autoregressive coeffcients and inverse error covariance matrix
  # maxit.beta : maximum iteration for the estimation of autoregressive coeffcients using the SPG algorithm
  # maxit.omega : maximum iteration for the estimation of inverse error covariance matrix using JGL 
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol.both : tolerance for the Multiclass VAR for joint estimation of autoregressive coeffcients and inverse error covariance matrix (if no default, take 0.00001)
  # tol.beta : tolerance for the estimation of autoregressive coeffcients using the SPG algorithm
  # tol.omega : tolerance for the estimation of inverse error covariance matrix using JGLautoregressive coeffcients
  # NBoot : number of bootstrap replicates
  # calculate.C : logical for calculating the C matrix. Default is FALSE.
  # parallell : logical for residual bootstrap procedure carried out using parallel computing. Default is FALSE (i.e. NO parallel).
  # type: "Lasso", "AdLasso", "GrLasso" for Lasso, Adapative Lasso and Group Lasso
  # groupindex : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda2_weights : logical. If TRUE, then do the weighting based on X_weight. Default is FALSE.
  # X_weight : K x w matrix of additional weights for lambda2. K is the number of classes and w the number of weights. Default is NULL.
  
  
  ##########
  # OUTPUT #
  ##########
  # BetaSign : significant coefficients based on standard errors built with the residual bootstrap
  
  # Elements from the MultiClass_VAR_generalOmega object
  Betahat_array<-FITOBJECT$Beta_array
  Betahat_vec<-FITOBJECT$Beta_new
  X<-FITOBJECT$X_array
  U<-FITOBJECT$E_array
  lam1.opt<-FITOBJECT$lambda1_opt
  lam2.opt<-FITOBJECT$lambda2_opt
  gam1.opt<-FITOBJECT$gamma1_opt
  gam2.opt<-FITOBJECT$gamma2_opt
  lamridge<-FITOBJECT$lambda_RIDGE
  J<-FITOBJECT$J
  K<-FITOBJECT$K
  P<-FITOBJECT$P
  n<-FITOBJECT$n
  
  if (is.null(CFIT)){
    CFIT<-FITOBJECT$CFIT
  }
  
  # 1. Center residuals
  Ucent<-apply(U,c(2,3),scale,scale=F,center=T) #U needs to be an array of dimension (n,J,K)
  
  BetaBoot<-matrix(NA,ncol=length(Betahat_vec),nrow=Nboot)
  
  if(parallell==T){# Parallell Computing
    
    checkpackage<-function(U){
      if((U %in% rownames(installed.packages()))==F){
        install.packages(U)
      }
    }
    packagelist<-list("foreach","doSNOW","doRNG")
    lapply(packagelist,checkpackage)
    # Load packages
    suppressMessages(suppressWarnings(packages <- lapply(packagelist, FUN = function(x) {
      library(x, character.only = TRUE)
    })))
    
    # Simulation with parallell computing
    ncl<-parallel:::detectCores()-2 #check number of cores, use two less than total number of cores
    R<-Nboot # number of simulations
    simestcl <- makeCluster(ncl) 
    registerDoSNOW(simestcl)
    
    BOOTFITS<- foreach(isim = 1:R, .packages=c('MASS', 'magic', 'rARPACK', 'Matrix', 'JGL', 'plyr')) %dopar% {
      
      # Additional functions to be sourced (The functions should be in a folder named "Functions")
      source('Functions/MultiClass_VAR.R')
      source('Functions/BootSE.R')
      source('Functions/FunctionsVisualTools.R')
      
      BOOTFIT<-BootAUX_GeneralOmega(nboot=isim,Betahat_array=Betahat_array,Betahat_vec=Betahat_vec,X.data=X,lam1.opt=lam1.opt,lam2.opt=lam2.opt,gam1.opt=gam1.opt,gam2.opt=gam2.opt,
                                    lamridge=lamridge,CFIT=CFIT,J=J,K=K,P=P,n=n,maxit.both=maxit.both,maxit.beta=maxit.beta,maxit.omega=maxit.omega,mu=mu,
                                    tol.both=tol.both,tol.beta=tol.beta,tol.omega=tol.omega,Ucent=Ucent,calculate.C=calculate.C,type=type,groupindex=groupindex,
                                    lambda_weights=lambda_weights, Clusterinfo=Clusterinfo)
      
      list("BOOTFIT"=BOOTFIT)
    }
    stopCluster(simestcl)
  }else{
    bootseq<-1:Nboot
    BOOTFITS<-lapply(bootseq,FUN=BootAUX_GeneralOmega,Betahat_array=Betahat_array,Betahat_vec=Betahat_vec,X.data=X,lam1.opt=lam1.opt,lam2.opt=lam2.opt,gam1.opt=gam1.opt,gam2.opt=gam2.opt,
                     lamridge=lamridge,CFIT=CFIT,J=J,K=K,P=P,n=n,maxit.both=maxit.both,maxit.beta=maxit.beta,maxit.omega=maxit.omega,mu=mu,
                     tol.both=tol.both,tol.beta=tol.beta,tol.omega=tol.omega,Ucent=Ucent,calculate.C=calculate.C,type=type,groupindex=groupindex,
                     lambda_weights=lambda_weights, Clusterinfo=Clusterinfo)
  }
  
  
  BetaBoot<-matrix(unlist(BOOTFITS),byrow=T,nrow=Nboot)
  
  # Standard interval
  Bootsd<-apply(BetaBoot,2,sd)
  coef.index<-1:length(Bootsd)
  FLAGnorm<-unlist(lapply(coef.index,tstatnormal,estimate=Betahat_vec,sd=Bootsd))
  BetaSign<-rep(0,length(Betahat_vec))
  BetaSign[FLAGnorm]<-Betahat_vec[FLAGnorm]
  
  out<-list("BetaSign"=BetaSign)
}

BootAUX_GeneralOmega<-function(nboot,Betahat_array,Betahat_vec,X.data,
                               lam1.opt,lam2.opt,gam1.opt,gam2.opt,lamridge,
                               CFIT,J,K,P,n,maxit.both,maxit.beta,maxit.omega,mu,
                               tol.both,tol.beta,tol.omega,Ucent,calculate.C,
                               criterion="BIC",type=NULL,groupindex=NULL,
                               lambda_weights=NULL, Clusterinfo=NULL,max.it.test=100){
  
  #### Auxiliary Function: Perform residual bootstrap when allowing for cross-error correlations in the Multi-class VAR ####
  
  #########
  # INPUT #
  #########
  # nboot : bootstrap run
  # Betahat_array : array of estimated beta coefficients (output of MultiClassGeneralOmega.R)
  # Betahat_vec : vector of estimated beta coefficients (output of MultiClassGeneralOmega.R)
  # X.data : array of inputs (output of MultiClassGeneralOmega.R)
  # lam1.opt : optimal lambda 1 
  # lam2.opt : optimal lambda 2 (fused lasso)
  # gam1.opt : optimal gamma 1
  # gam2.opt : optimal gamma 2 (fused lasso)
  # lamridge : array (1,1,K) containing the OPTIMAL regularization parameter on the autoregressive coeffcients in the RIDGE estimator for each CLASS
  # CFIT : output the function buildCmatrix.R
  # J : number of time series in each class
  # K : number of classes (K>1)
  # P : order of VAR
  # n : time series length
  # maxit.both : maximum iteration of the Multiclass VAR for joint estimation of autoregressive coeffcients and inverse error covariance matrix
  # maxit.beta : maximum iteration for the estimation of autoregressive coeffcients using the SPG algorithm
  # maxit.omega : maximum iteration for the estimation of inverse error covariance matrix using JGL 
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol.both : tolerance for the Multiclass VAR for joint estimation of autoregressive coeffcients and inverse error covariance matrix (if no default, take 0.00001)
  # tol.beta : tolerance for the estimation of autoregressive coeffcients using the SPG algorithm
  # tol.omega : tolerance for the estimation of inverse error covariance matrix using JGLautoregressive coeffcients
  # Ucent : centered residuals
  # calculate.C : logical for calculating the C matrix. 
  # criterion : criterion for the choice of the regularization parameters. Default is "BIC", otherwise set "AICc".
  # type: "Lasso", "AdLasso", "GrLasso" for Lasso, Adapative Lasso and Group Lasso
  # groupindex : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda2_weights : logical. If TRUE, then do the weighting based on X_weight. Default is FALSE.
  # X_weight : K x w matrix of additional weights for lambda2. K is the number of classes and w the number of weights. Default is NULL.
  # max.it.test : maximum iteration for the multivariate normality test. Default is 100.
  
  ##########
  # OUTPUT #
  ##########
  # BetaBoot : beta coefficients based on residual bootstrap 
  
  
  # Draw residuals with replacement (check for multivariate white noise)
  flag<-1
  it<-1
  while(length(flag)!=0 & (it<max.it.test)){
    n.index<-sample(1:n,n,replace=T)
    Uboot<-Ucent[n.index,,]
    
    pvalueMWN<-matrix(NA,ncol=dim(Uboot)[3],nrow=1)
    for(i.class in 1:dim(Uboot)[3]){
      MWNtest<-MultivariateWhiteNoise(Uboot[,,i.class], lag = floor(sqrt(dim(Uboot)[1])), adj = 0)
      pvalueMWN[1,i.class]<-MWNtest[nrow(MWNtest),4]
    }
    flag<-which(pvalueMWN<0.05)
    it<-it+1
  }
  
  
  # 2. Compute bootstrap time series
  BootDGP<-MultiClass_DG(Beta_DG=Betahat_array, Beta_vec_DG=Betahat_vec, X_DG=X.data, E_DG=Uboot,calculate.C=calculate.C)
  
  
  # 3. Multi-class fit with bootstrap time series, 
  n<-length(BootDGP$Y_input)/(J*K)
  BOOTFIT<-MultiClass_VAR_generalOmega(Data=NULL, calculate.C=F, X0=BootDGP$X0_input,XX_prod=BootDGP$XX_input,XY_prod=BootDGP$XY_input,Y=BootDGP$Y_input, X=BootDGP$X_input, C=CFIT$C, J=J, K=K, P=P, n=n, 
                                       CNorm=CFIT$CNorm, maxit.both=maxit.both, maxit.beta=maxit.beta,maxit.omega=maxit.omega,b_init= array(0,c(ncol(BootDGP$X_input),1)), mu=mu, 
                                       tol.both=tol.both, tol.beta=tol.beta,tol.omega=tol.omega, lambda1_OPT=lam1.opt, lambda2_OPT=lam2.opt, gamma1_OPT=gam1.opt, gamma2_OPT=gam2.opt,
                                       lambdaRIDGE_OPT=lamridge,
                                       criterion=criterion,type=type,group=groupindex,
                                       lambda_weights=lambda_weights, Clusterinfo=Clusterinfo,
                                       lambda1_min=NULL, lambda1_max=NULL, lambda1_steps=NULL, 
                                       lambda2_min=NULL, lambda2_max=NULL, lambda2_steps=NULL,
                                       gamma1_min=NULL, gamma1_max=NULL, gamma1_steps=NULL, 
                                       gamma2_min=NULL, gamma2_max=NULL, gamma2_steps=NULL,
                                       lambdaRIDGE_min=NULL, lambdaRIDGE_max=NULL, lambdaRIDGE_steps=NULL)
  
  BetaBoot<-BOOTFIT$Beta_new
  rm(BOOTFIT)
  return(BetaBoot)
  
}

